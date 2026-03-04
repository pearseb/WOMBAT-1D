from __future__ import annotations

import shutil
import subprocess
from pathlib import Path

import numpy as np
import xarray as xr

from .config import ColumnConfig, ModelConfig


class GotmDriver:
    """Adapter for running GOTM as the physical 1-D backbone."""

    KAPPA_CANDIDATES = ("nuh", "num", "Kz", "kz", "akt", "Av")

    def __init__(self, model_cfg: ModelConfig):
        self.cfg = model_cfg

    def available(self) -> bool:
        return shutil.which(self.cfg.gotm_executable) is not None

    def _stage_setup(self, run_dir: Path) -> None:
        run_dir.mkdir(parents=True, exist_ok=True)
        if self.cfg.gotm_setup_dir is None:
            return
        if not self.cfg.gotm_setup_dir.exists():
            raise FileNotFoundError(f"GOTM setup dir not found: {self.cfg.gotm_setup_dir}")
        for src in self.cfg.gotm_setup_dir.iterdir():
            dest = run_dir / src.name
            if src.is_dir():
                shutil.copytree(src, dest, dirs_exist_ok=True)
            else:
                shutil.copy2(src, dest)

    def run(self) -> xr.Dataset:
        run_dir = self.cfg.gotm_run_dir
        self._stage_setup(run_dir)
        subprocess.run([self.cfg.gotm_executable], cwd=run_dir, check=True)
        output = run_dir / self.cfg.gotm_output_path
        if not output.exists():
            raise FileNotFoundError(f"Expected GOTM output file not found: {output}")
        return xr.open_dataset(output)

    @classmethod
    def kappa_from_output(cls, gotm_output: xr.Dataset) -> xr.DataArray | None:
        for name in cls.KAPPA_CANDIDATES:
            if name in gotm_output:
                return gotm_output[name]
        return None


def _diffuse_single_step(prev: np.ndarray, kappa: np.ndarray, dz: float, dt: float) -> np.ndarray:
    nxt = prev.copy()
    interior = slice(1, -1)
    nxt[interior] = prev[interior] + dt * kappa[interior] * (
        prev[2:] - 2.0 * prev[1:-1] + prev[:-2]
    ) / (dz * dz)
    return nxt


def _estimate_dz(ds: xr.Dataset, column: ColumnConfig, nz: int) -> float:
    """Estimate vertical spacing using dataset depth coordinate when available."""

    if "depth" in ds.coords:
        depth = np.asarray(ds["depth"].values, dtype=float)
        if depth.ndim == 1 and depth.size > 1:
            dz = float(np.mean(np.diff(depth)))
            if np.isfinite(dz) and dz > 0.0:
                return dz
    # Depth grid in this repo is built with linspace(0, depth_m, nz), so spacing
    # is depth_m / (nz - 1), not depth_m / nz.
    return column.depth_m / max(nz - 1, 1)


def apply_gotm_mixing(ds: xr.Dataset, kappa_2d: np.ndarray, column: ColumnConfig) -> xr.Dataset:
    """Apply vertical diffusion using GOTM-provided diffusivity profile K(z,t)."""

    out = ds.copy(deep=True)
    if kappa_2d.ndim != 2:
        raise ValueError(f"Expected 2D kappa array [time, depth], got shape {kappa_2d.shape}")

    dt = column.dt_seconds

    for var in out.data_vars:
        arr = out[var].values
        if arr.ndim != 2:
            continue
        nt = min(arr.shape[0], kappa_2d.shape[0])
        nz = min(arr.shape[1], kappa_2d.shape[1])
        dz = _estimate_dz(out, column, nz)
        for t in range(1, nt):
            prev = arr[t - 1, :nz]
            kappat = np.maximum(kappa_2d[t - 1, :nz], 0.0)
            arr[t, :nz] = _diffuse_single_step(prev, kappat, dz, dt)
        out[var].values[:] = arr

    return out


def fallback_mix_advect(ds: xr.Dataset, column: ColumnConfig) -> xr.Dataset:
    """Fallback explicit diffusion when GOTM is not available."""

    out = ds.copy(deep=True)
    nz = int(ds.sizes.get("depth", column.nz))
    dz = _estimate_dz(ds, column, nz)
    kappa = 1e-4
    dt = column.dt_seconds

    for var in out.data_vars:
        arr = out[var].values
        if arr.ndim != 2:
            continue
        for t in range(1, arr.shape[0]):
            arr[t] = _diffuse_single_step(arr[t - 1], np.full(arr.shape[1], kappa), dz, dt)
        out[var].values[:] = arr
    return out
