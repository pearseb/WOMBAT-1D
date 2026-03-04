from pathlib import Path

import numpy as np
import xarray as xr

from womcol.config import RuntimeConfig
from womcol.simulation import WombatColumnModel


def _write_test_data(tmp_path: Path) -> tuple[Path, Path]:
    times = np.array(["2019-01-01", "2019-01-02"], dtype="datetime64[ns]")
    jra = xr.Dataset(
        {"sw_down": ("time", np.array([100.0, 120.0]))},
        coords={"time": times, "lat": -42.0, "lon": 147.0},
    )
    init = xr.Dataset(
        {
            "phy": (("depth",), np.linspace(0.1, 0.01, 8)),
            "no3": (("depth",), np.linspace(5.0, 20.0, 8)),
        },
        coords={"depth": np.linspace(0, 200, 8), "lat": -42.0, "lon": 147.0},
    )
    jra_path = tmp_path / "jra.nc"
    init_path = tmp_path / "init.nc"
    jra.to_netcdf(jra_path)
    init.to_netcdf(init_path)
    return jra_path, init_path


def test_smoke_run(tmp_path: Path) -> None:
    jra, init = _write_test_data(tmp_path)
    cfg = RuntimeConfig(
        year=2019,
        latitude=-42.0,
        longitude=147.0,
        data={"jra55_uri": str(jra), "init_uri": str(init)},
        model={"scheme": "wombat-lite", "use_gotm": False},
        column={"depth_m": 200, "nz": 8, "dt_seconds": 1800},
    )
    ds = WombatColumnModel(cfg).run()
    assert "phy" in ds
    assert ds.dims["depth"] == 8
