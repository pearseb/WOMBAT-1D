from __future__ import annotations

from pathlib import Path

import xarray as xr

from .config import DataConfig


def _open_dataset_with_backends(uri: str) -> xr.Dataset:
    try:
        return xr.open_dataset(uri)
    except ValueError as exc:
        msg = str(exc)
        if "xarray's IO backends" in msg:
            raise ValueError(
                f"Unable to open dataset '{uri}' because no compatible xarray NetCDF backend is installed. "
                "Install at least one backend, e.g. `python -m pip install netCDF4` "
                "(or `h5netcdf`/`scipy` depending on file format), then retry."
            ) from exc
        raise


def _select_nearest_point(ds: xr.Dataset, *, lat_name: str, lon_name: str, latitude: float, longitude: float) -> xr.Dataset:
    """Select nearest lat/lon when those coordinates are indexed dimensions.

    If lat/lon are scalar coordinates (common for single-point files), no selection
    is needed and this function returns the dataset unchanged.
    """

    indexers: dict[str, float] = {}
    if lat_name in ds.indexes:
        indexers[lat_name] = latitude
    if lon_name in ds.indexes:
        indexers[lon_name] = longitude

    if not indexers:
        return ds

    return ds.sel(indexers, method="nearest")


class ForcingProvider:
    """Loads surface forcing from JRA55-like datasets."""

    def __init__(self, data_config: DataConfig):
        self.cfg = data_config
        self.cfg.cache_dir.mkdir(parents=True, exist_ok=True)

    def load_jra55(self, *, year: int, latitude: float, longitude: float) -> xr.Dataset:
        if self.cfg.jra55_uri is None:
            raise ValueError("DataConfig.jra55_uri must point to a local/remote JRA55-compatible dataset.")

        ds = _open_dataset_with_backends(self.cfg.jra55_uri)
        lat_name = "lat" if "lat" in ds.coords else "latitude"
        lon_name = "lon" if "lon" in ds.coords else "longitude"

        sel = _select_nearest_point(ds, lat_name=lat_name, lon_name=lon_name, latitude=latitude, longitude=longitude)

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

        ds = _open_dataset_with_backends(self.cfg.init_uri)
        lat_name = "lat" if "lat" in ds.coords else "latitude"
        lon_name = "lon" if "lon" in ds.coords else "longitude"

        return _select_nearest_point(ds, lat_name=lat_name, lon_name=lon_name, latitude=latitude, longitude=longitude)


def write_netcdf(ds: xr.Dataset, out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    ds.to_netcdf(out_path)
