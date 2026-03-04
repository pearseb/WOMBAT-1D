from __future__ import annotations

import argparse
from pathlib import Path

import yaml

from .config import RuntimeConfig
from .simulation import WombatColumnModel


def main() -> None:
    parser = argparse.ArgumentParser(description="Run 1-D WOMBAT-lite/mid water-column model")
    parser.add_argument("config", type=Path, help="YAML config path")
    parser.add_argument("--output", type=Path, default=Path("output/womcol.nc"))
    args = parser.parse_args()

    cfg = RuntimeConfig(**yaml.safe_load(args.config.read_text()))
    WombatColumnModel(cfg).run_to_file(args.output)


if __name__ == "__main__":
    main()
