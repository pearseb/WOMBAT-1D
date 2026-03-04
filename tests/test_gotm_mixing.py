import numpy as np
import xarray as xr

from womcol.config import ColumnConfig
from womcol.gotm import apply_gotm_mixing


def test_apply_gotm_mixing_changes_profile() -> None:
    ds = xr.Dataset(
        {"phy": (("time", "depth"), np.array([[1.0, 0.5, 0.1], [1.0, 0.5, 0.1], [1.0, 0.5, 0.1]]))},
        coords={"time": np.arange(3), "depth": np.array([0.0, 50.0, 100.0])},
    )
    kappa = np.full((3, 3), 1e-4)
    out = apply_gotm_mixing(ds, kappa, ColumnConfig(depth_m=100, nz=3, dt_seconds=1800))
    assert float(out["phy"].isel(time=2, depth=1)) != 0.5
