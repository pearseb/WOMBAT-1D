import numpy as np
import xarray as xr

from womcol.config import RuntimeConfig
from womcol.simulation import WombatColumnModel, _limit_forcing_to_runtime_window


def test_limit_forcing_to_runtime_window_uses_days() -> None:
    time = np.array(
        np.arange(
            np.datetime64("2019-01-01T00:00"),
            np.datetime64("2019-01-01T00:00") + 240 * np.timedelta64(1, "h"),
            np.timedelta64(1, "h"),
        ),
        dtype="datetime64[ns]",
    )
    forcing = xr.Dataset({"sw_down": ("time", np.ones(time.size))}, coords={"time": time})

    out = _limit_forcing_to_runtime_window(forcing, days=3)
    assert out.sizes["time"] == 72


class _DummyForcingProvider:
    def __init__(self, forcing: xr.Dataset):
        self._forcing = forcing

    def load_jra55(self, **_: object) -> xr.Dataset:
        return self._forcing


class _DummyInitProvider:
    def __init__(self, init: xr.Dataset):
        self._init = init

    def load_profile(self, **_: object) -> xr.Dataset:
        return self._init


def test_run_sets_timestep_attrs() -> None:
    cfg = RuntimeConfig(
        year=2019,
        latitude=-42.0,
        longitude=147.0,
        days=2,
        model={"scheme": "wombat-lite", "use_gotm": False},
        column={"depth_m": 100.0, "nz": 10, "dt_seconds": 3600.0},
        data={"jra55_uri": "unused", "init_uri": "unused"},
    )

    time = np.array(
        np.arange(
            np.datetime64("2019-01-01T00:00"),
            np.datetime64("2019-01-01T00:00") + 240 * np.timedelta64(1, "h"),
            np.timedelta64(1, "h"),
        ),
        dtype="datetime64[ns]",
    )
    forcing = xr.Dataset({"sw_down": ("time", np.ones(time.size) * 100.0)}, coords={"time": time})
    depth = np.linspace(0.0, 100.0, 10)
    init = xr.Dataset(
        {
            "phy": ("depth", np.linspace(0.2, 0.02, 10)),
            "no3": ("depth", np.linspace(1.0, 10.0, 10)),
            "dfe": ("depth", np.linspace(2e-4, 8e-4, 10)),
            "zoo": ("depth", np.linspace(0.1, 0.01, 10)),
            "det": ("depth", np.linspace(0.05, 0.2, 10)),
            "o2": ("depth", np.linspace(220.0, 180.0, 10)),
        },
        coords={"depth": depth},
    )

    model = WombatColumnModel(cfg)
    model.forcing_provider = _DummyForcingProvider(forcing)
    model.init_provider = _DummyInitProvider(init)

    out = model.run()
    assert out.attrs["n_timesteps"] == 48
    assert out.attrs["configured_days"] == 2
    assert out.attrs["mixing_scheme"] == "fallback_diffusion"
