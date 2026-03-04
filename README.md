# WOMBAT-1D

Python scaffold for a 1-D water-column implementation of **WOMBAT-lite** and **WOMBAT-mid** with explicit **GOTM** integration for physical mixing/advection.

## What is implemented

- Runtime selection of `wombat-lite` vs `wombat-mid` biogeochemistry.
- User-configurable longitude, latitude, and year.
- Surface forcing loader for JRA55-compatible datasets (light/momentum/heat/freshwater fields can be read from user-provided files).
- Initialization loader for observational/climatological profile datasets at nearest lat/lon.
- GOTM integration:
  - supports running a compiled GOTM executable from this workflow,
  - stages a user-provided GOTM setup directory into a run directory,
  - reads GOTM output (`output.nc` by default) and uses diagnosed vertical diffusivity to mix tracers,
  - falls back to internal diffusion if GOTM is disabled/unavailable.
- Automation script + GitHub Action to sync canonical Fortran WOMBAT source files from:
  - `generic_WOMBATlite.F90`
  - `generic_WOMBATmid.F90`

## GOTM integration workflow

1. Install build dependencies (Ubuntu example):

```bash
sudo apt-get update
sudo apt-get install -y cmake gfortran build-essential libnetcdf-dev libnetcdff-dev
```

2. Build GOTM from source (official repo: https://github.com/gotm-model/code):

```bash
./scripts/install_gotm.sh
```

3. Prepare a GOTM setup directory (e.g., `./gotm_setup`) containing your GOTM namelists/YAML and forcing files.
4. Point model config keys to GOTM paths:
   - `model.gotm_executable`
   - `model.gotm_setup_dir`
   - `model.gotm_run_dir`
   - `model.gotm_output_path`

When `model.use_gotm: true`, WOMBAT-1D runs GOTM and applies vertical mixing using GOTM-reported diffusivity (`nuh/num/Kz/...`).

5. Validate GOTM was built correctly:

```bash
./scripts/check_gotm.sh ./third_party/gotm-code/build/gotm
```


## CI validation

A GitHub Actions workflow (`.github/workflows/gotm-smoke-test.yml`) now builds GOTM from the official repository and runs a binary smoke check on every push/PR.

## Quick start

```bash
python -m venv .venv
source .venv/bin/activate
python -m ensurepip --upgrade
python -m pip install --upgrade pip setuptools wheel
python -m pip install -e .
python scripts/sync_wombat_sources.py
python scripts/create_example_data.py
womcol-run config.example.yml --output output/womcol.nc
```


### Troubleshooting (HPC / older pip)

If `pip install -e .` fails with errors like:
- `No module named pip` during `setup.py develop`, or
- very old pip (e.g. `21.2.3`) warnings,

run:

```bash
source .venv/bin/activate
python -m ensurepip --upgrade
python -m pip install --upgrade pip setuptools wheel
python -m pip install -e .
```

Using `python -m pip` ensures the `pip` attached to the active virtual environment is used.

If you see `ModuleNotFoundError: No module named 'scipy'` during initialization interpolation, update to the latest code and reinstall:

```bash
git pull
python -m pip install -e .
```

If `womcol-run` fails with an xarray backend error like "found the following matches ... but their dependencies may not be installed", install a NetCDF backend in the same venv:

```bash
python -m pip install netCDF4
```

(Alternatives: `h5netcdf` or `scipy`, depending on the file format.)

You can replace `example_data/jra55_point.nc` with your own JRA55 forcing NetCDF file for a real location/year; that is often the easiest path for production runs.

If you generated example files before this fix and still see `KeyError: "no index found for coordinate 'lat'"`, regenerate them with:

```bash
python scripts/create_example_data.py
```


### Upstream link/symlink workflow

You can link this repo directly to a local clone of `ACCESS-NRI/GFDL-generic-tracers` (instead of downloading raw files), which is useful for seamless co-development:

```bash
# copy synced files from local upstream clone
python scripts/sync_wombat_sources.py --upstream-dir /path/to/GFDL-generic-tracers

# OR symlink to upstream files so this repo always points to latest local upstream checkout
python scripts/sync_wombat_sources.py --upstream-dir /path/to/GFDL-generic-tracers --symlink
```

In CI, the sync workflow now checks out the upstream repository directly and syncs from that checkout.

### WOMBAT-lite process alignment check

To compare the Python WOMBAT-lite implementation against a local copy of the canonical Fortran source, run:

```bash
python scripts/check_wombat_lite_alignment.py womcol/upstream/generic_WOMBATlite.F90 --strict
```

This provides a quick tracer/process token sanity check while iterating toward one-to-one parity.

## Notes

- This repository now includes concrete runtime hooks for GOTM execution and uptake of GOTM turbulence diagnostics.
- Full parity with WOMBAT Fortran source terms still requires completing one-to-one tendency translation/bindings in `womcol/wombat.py`.
