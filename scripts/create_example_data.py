#!/usr/bin/env python3
"""Create minimal self-contained example NetCDF files for WOMBAT-1D runs."""

from __future__ import annotations

from pathlib import Path

import numpy as np
import xarray as xr


def build_init_dataset(lat: float, lon: float, depth: np.ndarray) -> xr.Dataset:
    nz = depth.size

    def lin(top: float, bottom: float) -> np.ndarray:
        return np.linspace(top, bottom, nz)

    ds = xr.Dataset(
        {
            # WOMBAT-lite + WOMBAT-mid union of tracers
            "o2": (("depth",), lin(220.0, 160.0)),
            "no3": (("depth",), lin(1.0, 18.0)),
            "po4": (("depth",), lin(0.1, 1.5)),
            "dfe": (("depth",), lin(0.2e-3, 1.2e-3)),
            "phy": (("depth",), lin(0.2, 0.02)),
            "zoo": (("depth",), lin(0.08, 0.01)),
            "det": (("depth",), lin(0.03, 0.2)),
            "doc": (("depth",), lin(0.5, 0.8)),
        },
        coords={"depth": depth, "lat": lat, "lon": lon},
    )
    return ds


def build_forcing_dataset(lat: float, lon: float, year: int, nsteps: int = 48) -> xr.Dataset:
    time = np.array(
        np.arange(np.datetime64(f"{year}-01-01T00:00"), np.datetime64(f"{year}-01-01T00:00") + nsteps * np.timedelta64(1, "h"), np.timedelta64(1, "h")),
        dtype="datetime64[ns]",
    )

    phase = np.linspace(0, 2 * np.pi, nsteps)
    sw = np.maximum(0.0, 300.0 * np.sin(phase))

    ds = xr.Dataset(
        {
            # Used directly by WOMBAT-1D light limitation
            "sw_down": (("time",), sw),
            # Typical GOTM-style atmospheric forcing variables
            "lw_down": (("time",), 320.0 + 20.0 * np.cos(phase)),
            "tair": (("time",), 289.0 + 2.5 * np.sin(phase)),
            "qair": (("time",), 0.010 + 0.002 * np.cos(phase)),
            "u10": (("time",), 6.0 + 2.0 * np.sin(phase + 0.2)),
            "v10": (("time",), 2.0 + 1.0 * np.cos(phase + 0.4)),
            "sp": (("time",), 101325.0 + 800.0 * np.sin(phase / 2.0)),
            "precip": (("time",), 1.0e-7 * (1.0 + np.cos(phase))),
            "heat_flux": (("time",), 30.0 * np.sin(phase)),
            "freshwater_flux": (("time",), 1.0e-8 * np.cos(phase)),
            "tau_x": (("time",), 0.07 * np.sin(phase + 0.3)),
            "tau_y": (("time",), 0.04 * np.cos(phase + 0.1)),
        },
        coords={"time": time, "lat": lat, "lon": lon},
    )
    return ds


def main() -> None:
    out_dir = Path("example_data")
    out_dir.mkdir(parents=True, exist_ok=True)

    lat = -42.0
    lon = 147.0
    year = 2019
    depth = np.linspace(0.0, 200.0, 80)

    init_ds = build_init_dataset(lat=lat, lon=lon, depth=depth)
    forcing_ds = build_forcing_dataset(lat=lat, lon=lon, year=year, nsteps=24 * 10)

    init_ds.to_netcdf(out_dir / "init_profiles.nc")
    forcing_ds.to_netcdf(out_dir / "jra55_point.nc")

    print(f"Wrote {out_dir / 'init_profiles.nc'}")
    print(f"Wrote {out_dir / 'jra55_point.nc'}")


if __name__ == "__main__":
    main()
