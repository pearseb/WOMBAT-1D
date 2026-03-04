from __future__ import annotations

from pathlib import Path

import numpy as np
import xarray as xr

from .config import RuntimeConfig
from .gotm import GotmDriver, apply_gotm_mixing, fallback_mix_advect
from .io import ForcingProvider, InitialConditionProvider
from .wombat import WombatModel


def _column_inventory(ds: xr.Dataset) -> dict[str, float]:
    totals: dict[str, float] = {}
    for name, da in ds.data_vars.items():
        if da.ndim == 2 and "depth" in da.dims:
            totals[name] = float(np.asarray(da.isel(time=-1).sum().values))
    return totals


def _interp_profile_to_grid(var: xr.DataArray, target_z: np.ndarray) -> np.ndarray:
    """Interpolate (and linearly extrapolate) a 1-D profile onto target depths.

    Uses NumPy only to avoid optional SciPy dependency required by xarray.interp.
    """

    if "depth" not in var.coords:
        raise ValueError(f"Tracer '{var.name}' is missing required 'depth' coordinate in initialization dataset.")

    src_z = np.asarray(var["depth"].values, dtype=float)
    src_v = np.asarray(var.values, dtype=float)

    if src_z.ndim != 1 or src_v.ndim != 1:
        raise ValueError(f"Tracer '{var.name}' must be 1-D over depth for initialization.")

    order = np.argsort(src_z)
    src_z = src_z[order]
    src_v = src_v[order]

    interp = np.interp(target_z, src_z, src_v)

    if src_z.size >= 2:
        left_slope = (src_v[1] - src_v[0]) / (src_z[1] - src_z[0] + 1e-12)
        right_slope = (src_v[-1] - src_v[-2]) / (src_z[-1] - src_z[-2] + 1e-12)

        left_mask = target_z < src_z[0]
        right_mask = target_z > src_z[-1]

        interp[left_mask] = src_v[0] + left_slope * (target_z[left_mask] - src_z[0])
        interp[right_mask] = src_v[-1] + right_slope * (target_z[right_mask] - src_z[-1])

    return interp


def _limit_forcing_to_runtime_window(forcing: xr.Dataset, *, days: int) -> xr.Dataset:
    """Limit forcing to the configured runtime window if time is datetime-like."""

    if "time" not in forcing.coords:
        return forcing
    if forcing.sizes.get("time", 0) == 0:
        return forcing

    time_vals = np.asarray(forcing["time"].values)
    if np.issubdtype(time_vals.dtype, np.datetime64):
        start = np.datetime64(time_vals[0])
        stop = start + np.timedelta64(days, "D") - np.timedelta64(1, "ns")
        return forcing.sel(time=slice(start, stop))

    return forcing


class WombatColumnModel:
    def __init__(self, cfg: RuntimeConfig):
        self.cfg = cfg
        self.forcing_provider = ForcingProvider(cfg.data)
        self.init_provider = InitialConditionProvider(cfg.data)
        self.bio = WombatModel(cfg.model.scheme, cfg.model.parameters)
        self.gotm = GotmDriver(cfg.model)

    def _build_grid(self) -> np.ndarray:
        return np.linspace(0.0, self.cfg.column.depth_m, self.cfg.column.nz)

    def _run_biogeochemistry(self, forcing: xr.Dataset, init: xr.Dataset) -> xr.Dataset:
        time = forcing.time.values
        z = self._build_grid()
        tracers = self.bio.empty_state(self.cfg.column.nz)

        for name, arr in tracers.tracers.items():
            if name in init:
                source = _interp_profile_to_grid(init[name], z)
                tracers.tracers[name][:] = source

        light_name = "sw_down" if "sw_down" in forcing else list(forcing.data_vars)[0]

        history = {name: np.zeros((len(time), len(z)), dtype=float) for name in tracers.tracers}
        total_steps = len(time)
        for it, _ in enumerate(time):
            print(f"Running model at timestep {it + 1} of {total_steps} timesteps", flush=True)
            light_surface = float(forcing[light_name].isel(time=it).values)
            light_profile = np.exp(-0.04 * z) * max(light_surface, 0.0) / (abs(light_surface) + 1e-9)
            tracers = self.bio.step(tracers, dt_seconds=self.cfg.column.dt_seconds, light=light_profile)
            for name in history:
                history[name][it, :] = tracers.tracers[name]

        return xr.Dataset(
            {name: (("time", "depth"), values) for name, values in history.items()},
            coords={"time": time, "depth": z},
        )

    def _mix_and_advect(self, ds: xr.Dataset) -> xr.Dataset:
        pre_inventory = _column_inventory(ds)

        gotm_ok, gotm_status = self.gotm.availability_status()
        if self.cfg.model.use_gotm and gotm_ok:
            gotm_ds = self.gotm.run()
            kappa = self.gotm.kappa_from_output(gotm_ds)
            if kappa is not None and kappa.ndim >= 2:
                kappa_vals = np.asarray(kappa.values)
                if kappa_vals.ndim > 2:
                    kappa_vals = kappa_vals.reshape(kappa_vals.shape[0], -1)
                mixed = apply_gotm_mixing(ds, kappa_vals, self.cfg.column)
                mixed.attrs["mixing_scheme"] = "gotm"
                mixed.attrs["gotm_status"] = gotm_status
                mixed.attrs["gotm_kappa_var"] = kappa.name or "unknown"
                mixed.attrs["inventory_before_mixing"] = str(pre_inventory)
                mixed.attrs["inventory_after_mixing"] = str(_column_inventory(mixed))
                return mixed

        mixed = fallback_mix_advect(ds, self.cfg.column)
        mixed.attrs["mixing_scheme"] = "fallback_diffusion"
        if self.cfg.model.use_gotm:
            mixed.attrs["gotm_status"] = gotm_status
        mixed.attrs["inventory_before_mixing"] = str(pre_inventory)
        mixed.attrs["inventory_after_mixing"] = str(_column_inventory(mixed))
        return mixed

    def run(self) -> xr.Dataset:
        forcing = self.forcing_provider.load_jra55(
            year=self.cfg.year,
            latitude=self.cfg.latitude,
            longitude=self.cfg.longitude,
        )
        forcing = _limit_forcing_to_runtime_window(forcing, days=self.cfg.days)
        init = self.init_provider.load_profile(latitude=self.cfg.latitude, longitude=self.cfg.longitude)
        bio_ds = self._run_biogeochemistry(forcing, init)
        mixed = self._mix_and_advect(bio_ds)
        mixed.attrs["n_timesteps"] = int(mixed.sizes.get("time", 0))
        mixed.attrs["configured_days"] = int(self.cfg.days)
        if mixed.sizes.get("time", 0) > 0:
            mixed.attrs["time_start"] = str(np.asarray(mixed["time"].values[0]))
            mixed.attrs["time_end"] = str(np.asarray(mixed["time"].values[-1]))
        return mixed

    def run_to_file(self, out_path: Path) -> xr.Dataset:
        ds = self.run()
        out_path.parent.mkdir(parents=True, exist_ok=True)
        ds.to_netcdf(out_path)
        return ds
