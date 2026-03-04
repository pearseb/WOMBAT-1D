#!/usr/bin/env python3
"""Sync canonical WOMBAT Fortran files from ACCESS-NRI/GFDL-generic-tracers.

This script downloads generic_WOMBATlite.F90 and generic_WOMBATmid.F90 and stores
snapshots under womcol/upstream/. A generated metadata file is used by CI to detect
upstream changes and keep this repository reproducible.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path
from urllib.request import urlopen

BASE = "https://raw.githubusercontent.com/ACCESS-NRI/GFDL-generic-tracers/main/generic_tracers"
FILES = ["generic_WOMBATlite.F90", "generic_WOMBATmid.F90"]


def main() -> None:
    out_dir = Path("womcol/upstream")
    out_dir.mkdir(parents=True, exist_ok=True)

    metadata: dict[str, dict[str, str]] = {}

    for fname in FILES:
        url = f"{BASE}/{fname}"
        with urlopen(url) as resp:
            text = resp.read().decode("utf-8")

        digest = hashlib.sha256(text.encode("utf-8")).hexdigest()
        local = out_dir / fname
        local.write_text(text)
        metadata[fname] = {"url": url, "sha256": digest}

    (out_dir / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")


if __name__ == "__main__":
    main()
