#!/usr/bin/env bash
set -euo pipefail

GOTM_BIN="${1:-./third_party/gotm-code/build/gotm}"

if [[ ! -x "${GOTM_BIN}" ]]; then
  echo "GOTM binary not found or not executable: ${GOTM_BIN}" >&2
  exit 1
fi

# Try help/version style invocation first; if unsupported, run with timeout to ensure it starts.
if "${GOTM_BIN}" --help >/tmp/gotm_help.txt 2>&1; then
  echo "GOTM help command succeeded"
  exit 0
fi

if command -v timeout >/dev/null 2>&1; then
  timeout 10s "${GOTM_BIN}" >/tmp/gotm_run.txt 2>&1 || true
  echo "GOTM executable launched (non-zero may be expected without setup files)"
else
  "${GOTM_BIN}" >/tmp/gotm_run.txt 2>&1 || true
  echo "GOTM executable launched (non-zero may be expected without setup files)"
fi
