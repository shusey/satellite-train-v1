#!/usr/bin/env python3
"""Run the required three cases and build every v1 report artifact."""

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
from sat_string.report import build_standard_report
from sat_string.simulator import STANDARD_CASES, simulate_case


def _logger(output: Path) -> logging.Logger:
    logger = logging.getLogger("sat_string.standard_suite")
    logger.handlers.clear()
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter("%(asctime)s %(levelname)s %(message)s")
    for handler in (
        logging.FileHandler(output / "run.log", encoding="utf-8"),
        logging.StreamHandler(),
    ):
        handler.setFormatter(formatter)
        logger.addHandler(handler)
    return logger


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--overwrite", "-overwrite", action="store_true")
    parser.add_argument("--skip-animation", action="store_true")
    args = parser.parse_args(argv)

    config = load_config(args.config)
    output = prepare_output_directory(args.output, overwrite=args.overwrite)
    dump_normalized_config(config, output / "normalized_config.yaml")
    logger = _logger(output)
    suite_start = perf_counter()
    try:
        for case_index, case in enumerate(STANDARD_CASES):
            case_start = perf_counter()
            logger.info("Starting case %d/3: %s", case_index + 1, case.name)
            case_log_lines = [
                f"case={case.name}",
                "status=running",
                f"control_enabled={case.control_enabled}",
                f"disturbance_enabled={case.disturbance_enabled}",
            ]

            def progress(time_s: float, final_s: float, name=case.name) -> None:
                elapsed = perf_counter() - case_start
                fraction = time_s / final_s
                remaining = elapsed * (1.0 - fraction) / fraction if fraction > 0.0 else float("inf")
                logger.info(
                    "case=%s progress=%.0f%% simulated=%.3f h estimated_remaining=%.1f s",
                    name,
                    100.0 * fraction,
                    time_s / 3600.0,
                    remaining,
                )
                case_log_lines.append(
                    f"progress={100.0 * fraction:.0f}% simulated_time_s={time_s:.9g} "
                    f"estimated_remaining_s={remaining:.3f}"
                )

            def case_event_log(message: str, name=case.name) -> None:
                logger.info("case=%s %s", name, message)
                case_log_lines.append(message)

            result = simulate_case(
                config,
                case,
                progress_callback=progress,
                log_callback=case_event_log,
            )
            case_output = output / case.name
            save_case(result, config, case_output)
            case_log_lines.extend(
                [
                    "status=completed",
                    f"elapsed_s={perf_counter() - case_start:.6f}",
                ]
            )
            with (case_output / "run.log").open("w", encoding="utf-8", newline="\n") as stream:
                stream.write("\n".join(case_log_lines) + "\n")
            logger.info("Completed case=%s in %.2f s", case.name, perf_counter() - case_start)
            del result

        logger.info("Building metrics, 300 dpi figures, and animation from saved histories")
        build_standard_report(
            output,
            config,
            include_animation=(
                config.visualization.animation_enabled and not args.skip_animation
            ),
            dpi=300,
        )
        logger.info("Standard suite completed in %.2f s", perf_counter() - suite_start)
        return 0
    except Exception:
        logger.exception(
            "Standard suite failed; completed numerical and figure artifacts remain in %s",
            output,
        )
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
