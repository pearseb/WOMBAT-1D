import numpy as np
import xarray as xr

from womcol.config import ColumnConfig
from womcol.gotm import apply_gotm_mixing


def test_apply_gotm_mixing_changes_profile() -> None:
    nz = 10
    init = np.zeros(nz)
    init[: nz // 2] = 1.0
    ds = xr.Dataset(
        {"phy": (("time", "depth"), np.vstack([init, np.zeros(nz), np.zeros(nz)]))},
        coords={"time": np.arange(3), "depth": np.linspace(0.0, 100.0, nz)},
    )
    kappa = np.full((3, nz), 1e-4)
    out = apply_gotm_mixing(ds, kappa, ColumnConfig(depth_m=100, nz=nz, dt_seconds=1800))
    assert float(out["phy"].isel(time=2, depth=4)) != float(ds["phy"].isel(time=0, depth=4))


def test_apply_gotm_mixing_conserves_column_inventory_with_no_flux_boundaries() -> None:
    nz = 10
    profile = np.zeros(nz)
    profile[2] = 2.0
    profile[7] = 1.0
    ds = xr.Dataset(
        {"phy": (("time", "depth"), np.vstack([profile, np.zeros(nz), np.zeros(nz), np.zeros(nz)]))},
        coords={"time": np.arange(4), "depth": np.linspace(0.0, 90.0, nz)},
    )
    kappa = np.full((4, nz), 1e-3)
    out = apply_gotm_mixing(ds, kappa, ColumnConfig(depth_m=90, nz=nz, dt_seconds=60.0))

    initial_inventory = float(ds["phy"].isel(time=0).sum())
    final_inventory = float(out["phy"].isel(time=-1).sum())
    assert abs(final_inventory - initial_inventory) < 1e-5
