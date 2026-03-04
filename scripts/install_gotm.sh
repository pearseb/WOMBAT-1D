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

# NetCDF discovery (important on HPC systems where headers/libs are in module prefixes).
# Honor user-provided NetCDF_ROOT/NETCDF_ROOT if present and also infer from nc-config.
NETCDF_HINT="${NetCDF_ROOT:-${NETCDF_ROOT:-}}"
if [[ -z "${NETCDF_HINT}" ]] && command -v nc-config >/dev/null 2>&1; then
  NETCDF_HINT="$(nc-config --prefix)"
fi

CMAKE_ARGS=(
  -DCMAKE_BUILD_TYPE=Release
)

if [[ -n "${NETCDF_HINT}" ]]; then
  CMAKE_ARGS+=(
    "-DNetCDF_ROOT=${NETCDF_HINT}"
    "-DNetCDF_DIR=${NETCDF_HINT}"
  )
  echo "Using NetCDF hint: ${NETCDF_HINT}"
else
  echo "No NetCDF hint detected (nc-config not found; NetCDF_ROOT/NETCDF_ROOT unset)."
  echo "If CMake fails to find NetCDF, set NetCDF_ROOT or load your NetCDF module first."
fi

cmake -S "${GOTM_SRC_DIR}" -B "${BUILD_DIR}" "${CMAKE_ARGS[@]}"
cmake --build "${BUILD_DIR}" -j

echo "GOTM binary should be at: ${BUILD_DIR}/gotm"
