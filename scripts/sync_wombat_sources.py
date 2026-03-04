#!/usr/bin/env python3
"""Sync canonical WOMBAT Fortran files from ACCESS-NRI/GFDL-generic-tracers.

Supports three modes:
1) Download from GitHub raw URLs (default)
2) Copy from a local upstream clone via --upstream-dir
3) Create symlinks to a local upstream clone via --upstream-dir --symlink
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
from urllib.request import urlopen

BASE = "https://raw.githubusercontent.com/ACCESS-NRI/GFDL-generic-tracers/main/generic_tracers"
FILES = ["generic_WOMBATlite.F90", "generic_WOMBATmid.F90"]


def _sha256_text(text: str) -> str:
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def _download_source(fname: str) -> tuple[str, dict[str, str]]:
    url = f"{BASE}/{fname}"
    with urlopen(url) as resp:
        text = resp.read().decode("utf-8")
    return text, {"url": url, "mode": "download"}


def _read_local_source(upstream_dir: Path, fname: str) -> tuple[str, dict[str, str], Path]:
    src = upstream_dir / "generic_tracers" / fname
    if not src.exists():
        raise FileNotFoundError(f"Missing upstream source file: {src}")
    text = src.read_text()
    return text, {"source": str(src), "mode": "local"}, src


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--upstream-dir", type=Path, default=None, help="Path to local clone of ACCESS-NRI/GFDL-generic-tracers")
    parser.add_argument("--symlink", action="store_true", help="Symlink files from --upstream-dir instead of copying")
    args = parser.parse_args()

    out_dir = Path("womcol/upstream")
    out_dir.mkdir(parents=True, exist_ok=True)

    metadata: dict[str, dict[str, str]] = {}

    for fname in FILES:
        local = out_dir / fname

        if args.upstream_dir is not None:
            text, info, src = _read_local_source(args.upstream_dir, fname)
            if args.symlink:
                if local.exists() or local.is_symlink():
                    local.unlink()
                # relative symlink keeps repo movable
                rel = os.path.relpath(src, start=local.parent)
                local.symlink_to(rel)
                digest = _sha256_text(src.read_text())
                metadata[fname] = {**info, "sha256": digest, "symlink": str(local)}
                continue
            local.write_text(text)
            digest = _sha256_text(text)
            metadata[fname] = {**info, "sha256": digest}
            continue

        text, info = _download_source(fname)
        local.write_text(text)
        digest = _sha256_text(text)
        metadata[fname] = {**info, "sha256": digest}

    (out_dir / "metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")


if __name__ == "__main__":
    main()
