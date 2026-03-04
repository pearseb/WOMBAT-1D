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

cmake -S "${GOTM_SRC_DIR}" -B "${BUILD_DIR}" -DCMAKE_BUILD_TYPE=Release
cmake --build "${BUILD_DIR}" -j

echo "GOTM binary should be at: ${BUILD_DIR}/gotm"
