from __future__ import annotations

from pathlib import Path

import xarray as xr

from .config import DataConfig


class ForcingProvider:
    """Loads surface forcing from JRA55-like datasets."""

    def __init__(self, data_config: DataConfig):
        self.cfg = data_config
        self.cfg.cache_dir.mkdir(parents=True, exist_ok=True)

    def load_jra55(self, *, year: int, latitude: float, longitude: float) -> xr.Dataset:
        if self.cfg.jra55_uri is None:
            raise ValueError("DataConfig.jra55_uri must point to a local/remote JRA55-compatible dataset.")

        ds = xr.open_dataset(self.cfg.jra55_uri)
        lat_name = "lat" if "lat" in ds.coords else "latitude"
        lon_name = "lon" if "lon" in ds.coords else "longitude"

        sel = ds.sel(
            {
                lat_name: latitude,
                lon_name: longitude,
            },
            method="nearest",
        )

        if "time" in sel.coords:
            sel = sel.sel(time=slice(f"{year}-01-01", f"{year}-12-31T23:59:59"))

        return sel


class InitialConditionProvider:
    """Loads profile initialization from climatology/observations."""

    def __init__(self, data_config: DataConfig):
        self.cfg = data_config

    def load_profile(self, *, latitude: float, longitude: float) -> xr.Dataset:
        if self.cfg.init_uri is None:
            raise ValueError("DataConfig.init_uri must point to a local/remote profile climatology dataset.")

        ds = xr.open_dataset(self.cfg.init_uri)
        lat_name = "lat" if "lat" in ds.coords else "latitude"
        lon_name = "lon" if "lon" in ds.coords else "longitude"

        return ds.sel({lat_name: latitude, lon_name: longitude}, method="nearest")


def write_netcdf(ds: xr.Dataset, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    ds.to_netcdf(out_path)
