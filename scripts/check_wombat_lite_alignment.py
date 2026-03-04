#!/usr/bin/env python3
"""Heuristic alignment checker between Python WOMBAT-lite and Fortran reference."""

from __future__ import annotations

import argparse
import re
from pathlib import Path

from womcol.wombat import DEFAULT_TRACERS

EXPECTED_PROCESSES = ["growth", "graz", "mort", "remin", "no3", "dfe", "o2"]


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("fortran_file", type=Path)
    parser.add_argument("--strict", action="store_true", help="Exit non-zero if expected process tokens are missing")
    args = parser.parse_args()

    txt = args.fortran_file.read_text().lower()
    py_tr = set(DEFAULT_TRACERS["wombat-lite"])

    print("Python WOMBAT-lite tracers:", sorted(py_tr))

    missing_proc = [p for p in EXPECTED_PROCESSES if re.search(p, txt) is None]
    if missing_proc:
        msg = f"Missing process tokens in Fortran: {missing_proc}"
        if args.strict:
            raise SystemExit(msg)
        print("WARNING:", msg)
    else:
        print("Found expected process tokens in Fortran file.")


if __name__ == "__main__":
    main()
