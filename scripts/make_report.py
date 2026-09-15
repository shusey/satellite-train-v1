#!/usr/bin/env python3
"""Regenerate metrics, figures, and animation from a saved standard suite."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
import sys

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from sat_string.config import load_config
from sat_string.report import build_standard_report


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", required=True, type=Path)
    parser.add_argument("--skip-animation", action="store_true")
    args = parser.parse_args(argv)
    logger = logging.getLogger("sat_string.make_report")
    logger.setLevel(logging.INFO)
    handler = logging.StreamHandler()
    handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
    logger.addHandler(handler)
    file_handler = logging.FileHandler(args.run / "run.log", mode="a", encoding="utf-8")
    file_handler.setFormatter(logging.Formatter("%(asctime)s %(levelname)s %(message)s"))
    logger.addHandler(file_handler)
    try:
        config = load_config(args.run / "normalized_config.yaml")
        logger.info("Regenerating report from saved arrays; no simulation will run")
        build_standard_report(
            args.run,
            config,
            include_animation=(
                config.visualization.animation_enabled and not args.skip_animation
            ),
        )
        logger.info("Report regeneration completed")
        return 0
    except Exception:
        logger.exception("Report regeneration failed")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
