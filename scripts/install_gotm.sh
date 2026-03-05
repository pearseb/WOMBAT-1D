#!/usr/bin/env bash
set -euo pipefail

# Build GOTM from source: https://github.com/gotm-model/code
ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
THIRD_PARTY_DIR="${ROOT_DIR}/third_party"
GOTM_SRC_DIR="${THIRD_PARTY_DIR}/gotm-code"
BUILD_DIR="${GOTM_SRC_DIR}/build"

mkdir -p "${THIRD_PARTY_DIR}"

if [[ ! -d "${GOTM_SRC_DIR}" ]]; then
  git clone --recurse-submodules https://github.com/gotm-model/code.git "${GOTM_SRC_DIR}"
else
  git -C "${GOTM_SRC_DIR}" pull --ff-only
fi

# GOTM requires submodules (e.g., extern/gsw). Ensure they are initialized/updated.
git -C "${GOTM_SRC_DIR}" submodule update --init --recursive

# NetCDF C discovery (important on HPC systems where headers/libs are in module prefixes).
NETCDF_HINT="${NetCDF_ROOT:-${NETCDF_ROOT:-}}"
NETCDF_INCLUDE_DIR=""
NETCDF_LIB_DIR=""

if command -v nc-config >/dev/null 2>&1; then
  if [[ -z "${NETCDF_HINT}" ]]; then
    NETCDF_HINT="$(nc-config --prefix)"
  fi
  NETCDF_INCLUDE_DIR="$(nc-config --includedir || true)"
  NETCDF_LIB_DIR="$(nc-config --libs 2>/dev/null | tr ' ' '\n' | sed -n 's/^-L//p' | head -n1)"
fi

# NetCDF Fortran discovery (needed for netcdf.mod).
NETCDFF_HINT="${NetCDF_Fortran_ROOT:-${NETCDF_FORTRAN_ROOT:-}}"
NETCDFF_INCLUDE_DIR=""
NETCDFF_LIB_DIR="${NETCDF_FORTRAN_LIB_DIR:-}"
NETCDFF_MOD_DIR="${NETCDF_FORTRAN_MOD_DIR:-}"
NETCDFF_LIBRARY=""
NF_FFLAGS=""

if command -v nf-config >/dev/null 2>&1; then
  if [[ -z "${NETCDFF_HINT}" ]]; then
    NETCDFF_HINT="$(nf-config --prefix)"
  fi
  NETCDFF_INCLUDE_DIR="$(nf-config --includedir || true)"
  NETCDFF_LIB_DIR="$(nf-config --flibs 2>/dev/null | tr ' ' '\n' | sed -n 's/^-L//p' | head -n1)"
  NF_FFLAGS="$(nf-config --fflags 2>/dev/null || true)"
fi


# Try to locate actual netcdff library file for explicit linking hints.
if command -v nf-config >/dev/null 2>&1; then
  while IFS= read -r ldir; do
    if [[ -z "${NETCDFF_LIB_DIR}" ]]; then
      NETCDFF_LIB_DIR="${ldir}"
    fi
    for cand in "${ldir}/libnetcdff.so" "${ldir}/libnetcdff.a" "${ldir}/libnetcdff.dylib"; do
      if [[ -f "${cand}" ]]; then
        NETCDFF_LIBRARY="${cand}"
        break 2
      fi
    done
  done < <(nf-config --flibs 2>/dev/null | tr ' ' '\n' | sed -n 's/^-L//p')
fi

# Prefer explicit module-dir override if provided.
if [[ -n "${NETCDFF_MOD_DIR}" && -f "${NETCDFF_MOD_DIR}/netcdf.mod" ]]; then
  NETCDFF_INCLUDE_DIR="${NETCDFF_MOD_DIR}"
fi

# If nf-config --includedir does not contain netcdf.mod, parse nf-config --fflags for -I dirs.
if [[ -n "${NF_FFLAGS}" && ( -z "${NETCDFF_INCLUDE_DIR}" || ! -f "${NETCDFF_INCLUDE_DIR}/netcdf.mod" ) ]]; then
  while IFS= read -r idir; do
    if [[ -n "${idir}" && -f "${idir}/netcdf.mod" ]]; then
      NETCDFF_INCLUDE_DIR="${idir}"
      break
    fi
  done < <(printf "%s\n" "${NF_FFLAGS}" | tr " " "\n" | sed -n "s/^-I//p")
fi

# Fallback: netcdf.mod may be colocated with NetCDF C headers on some installs.
if [[ -z "${NETCDFF_INCLUDE_DIR}" && -n "${NETCDF_INCLUDE_DIR}" && -f "${NETCDF_INCLUDE_DIR}/netcdf.mod" ]]; then
  NETCDFF_INCLUDE_DIR="${NETCDF_INCLUDE_DIR}"
fi

# Fail early with actionable guidance instead of failing mid-build with
# "Can't open module file 'netcdf.mod'".
if [[ -z "${NETCDFF_INCLUDE_DIR}" || ! -f "${NETCDFF_INCLUDE_DIR}/netcdf.mod" ]]; then
  cat >&2 <<EOF
ERROR: NetCDF Fortran module file (netcdf.mod) was not found.

Detected:
  nc-config: $(command -v nc-config || echo 'not found')
  nf-config: $(command -v nf-config || echo 'not found')
  NetCDF include dir: ${NETCDF_INCLUDE_DIR:-<unset>}
  NetCDF Fortran include dir: ${NETCDFF_INCLUDE_DIR:-<unset>}
  nf-config --fflags: ${NF_FFLAGS:-<unset>}
  NetCDF Fortran lib dir: ${NETCDFF_LIB_DIR:-<unset>}
  NetCDF Fortran library: ${NETCDFF_LIBRARY:-<unset>}

GOTM/FABM needs NetCDF Fortran headers/modules. Please load/install NetCDF Fortran
(e.g. module load netcdf-fortran) or set one of:
  NetCDF_Fortran_ROOT / NETCDF_FORTRAN_ROOT
  NETCDF_FORTRAN_MOD_DIR
  NETCDF_FORTRAN_LIB_DIR
so that <prefix>/include/netcdf.mod exists (and optionally set NETCDF_FORTRAN_LIB_DIR).
EOF
  exit 2
fi

CMAKE_ARGS=(
  -DCMAKE_BUILD_TYPE=Release
)

if [[ -n "${NETCDF_HINT}" ]]; then
  CMAKE_ARGS+=(
    "-DNetCDF_ROOT=${NETCDF_HINT}"
    "-DNetCDF_DIR=${NETCDF_HINT}"
    "-DCMAKE_PREFIX_PATH=${NETCDF_HINT}${CMAKE_PREFIX_PATH:+;${CMAKE_PREFIX_PATH}}"
  )
  echo "Using NetCDF C hint: ${NETCDF_HINT}"
fi

if [[ -n "${NETCDF_INCLUDE_DIR}" ]]; then
  CMAKE_ARGS+=("-DNetCDF_INCLUDE_DIRS=${NETCDF_INCLUDE_DIR}")
  echo "Using NetCDF C include dir: ${NETCDF_INCLUDE_DIR}"
fi
if [[ -n "${NETCDF_LIB_DIR}" ]]; then
  CMAKE_ARGS+=("-DNetCDF_LIBRARY_DIR=${NETCDF_LIB_DIR}")
  echo "Using NetCDF C lib dir: ${NETCDF_LIB_DIR}"
fi

if [[ -n "${NETCDFF_HINT}" ]]; then
  CMAKE_ARGS+=(
    "-DNetCDF_Fortran_ROOT=${NETCDFF_HINT}"
    "-DNETCDF_FORTRAN_ROOT=${NETCDFF_HINT}"
  )
  echo "Using NetCDF Fortran hint: ${NETCDFF_HINT}"
fi

# Ensure Fortran compiler can find netcdf.mod even when CMake FindNetCDF misses it.
if [[ -n "${NETCDFF_INCLUDE_DIR}" ]]; then
  CMAKE_ARGS+=(
    "-DCMAKE_Fortran_FLAGS=-I${NETCDFF_INCLUDE_DIR} ${CMAKE_Fortran_FLAGS:-}"
    "-DNetCDF_Fortran_INCLUDE_DIR=${NETCDFF_INCLUDE_DIR}"
    "-DNETCDF_FORTRAN_INCLUDE_DIR=${NETCDFF_INCLUDE_DIR}"
  )
  echo "Using NetCDF Fortran include dir: ${NETCDFF_INCLUDE_DIR}"
fi
if [[ -n "${NETCDFF_LIB_DIR}" ]]; then
  CMAKE_ARGS+=("-DNetCDF_Fortran_LIBRARY_DIR=${NETCDFF_LIB_DIR}")
  echo "Using NetCDF Fortran lib dir: ${NETCDFF_LIB_DIR}"
fi
if [[ -n "${NETCDFF_LIBRARY}" ]]; then
  CMAKE_ARGS+=(
    "-DNetCDF_Fortran_LIBRARY=${NETCDFF_LIBRARY}"
    "-DNETCDF_FORTRAN_LIBRARY=${NETCDFF_LIBRARY}"
  )
  echo "Using NetCDF Fortran library: ${NETCDFF_LIBRARY}"
fi

# Help final link step find -lnetcdff on HPC stacks where transitive -L is dropped.
if [[ -n "${NETCDFF_LIB_DIR}" ]]; then
  CMAKE_ARGS+=(
    "-DCMAKE_EXE_LINKER_FLAGS=-L${NETCDFF_LIB_DIR} ${CMAKE_EXE_LINKER_FLAGS:-}"
    "-DCMAKE_SHARED_LINKER_FLAGS=-L${NETCDFF_LIB_DIR} ${CMAKE_SHARED_LINKER_FLAGS:-}"
  )
fi

if [[ -z "${NETCDF_HINT}" && -z "${NETCDF_INCLUDE_DIR}" ]]; then
  echo "No NetCDF C hint detected (nc-config not found; NetCDF_ROOT/NETCDF_ROOT unset)."
fi

cmake -S "${GOTM_SRC_DIR}" -B "${BUILD_DIR}" "${CMAKE_ARGS[@]}"
cmake --build "${BUILD_DIR}" -j

echo "GOTM binary should be at: ${BUILD_DIR}/gotm"
