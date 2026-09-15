#!/usr/bin/env python3
"""Run and save one simulation case."""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
import sys
from time import perf_counter

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from sat_string.config import dump_normalized_config, load_config
from sat_string.io import prepare_output_directory, save_case
from sat_string.metrics import compute_case_metrics, save_metrics
from sat_string.simulator import CaseDefinition, simulate_case


def _logger(output: Path) -> logging.Logger:
    logger = logging.getLogger("sat_string.run_case")
    logger.handlers.clear()
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
    file_handler = logging.FileHandler(output / "run.log", encoding="utf-8")
    file_handler.setFormatter(formatter)
    stream_handler = logging.StreamHandler()
    stream_handler.setFormatter(formatter)
    logger.addHandler(file_handler)
    logger.addHandler(stream_handler)
    return logger


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--name", default="controlled_disturbed")
    parser.add_argument("--no-control", action="store_true")
    parser.add_argument("--no-disturbance", action="store_true")
    parser.add_argument("--overwrite", "-overwrite", action="store_true")
    args = parser.parse_args(argv)

    config = load_config(args.config)
    output = prepare_output_directory(args.output, overwrite=args.overwrite)
    dump_normalized_config(config, output / "normalized_config.yaml")
    logger = _logger(output)
    start = perf_counter()
    case = CaseDefinition(args.name, not args.no_control, not args.no_disturbance)

    def progress(time_s: float, final_s: float) -> None:
        elapsed = perf_counter() - start
        fraction = time_s / final_s
        remaining = elapsed * (1.0 - fraction) / fraction if fraction > 0.0 else float("inf")
        logger.info(
            "progress=%.0f%% simulated=%.3f h elapsed=%.1f s estimated_remaining=%.1f s",
            100.0 * fraction,
            time_s / 3600.0,
            elapsed,
            remaining,
        )

    try:
        logger.info(
            "Starting case=%s control=%s disturbance=%s",
            case.name,
            case.control_enabled,
            case.disturbance_enabled,
        )
        result = simulate_case(
            config,
            case,
            progress_callback=progress,
            log_callback=logger.info,
        )
        save_case(result, config, output)
        metrics = compute_case_metrics(result.arrays(), config, result.switch_events)
        save_metrics(metrics, output)
        logger.info("Completed in %.2f s", perf_counter() - start)
        return 0
    except Exception:
        logger.exception("Run failed")
        return 1


if __name__ == "__main__":
    raise SystemExit(main())

