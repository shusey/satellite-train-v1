#!/usr/bin/env python3
"""Sweep controller thresholds and analyze L-to-H switching-wave propagation."""

from __future__ import annotations

import argparse
from dataclasses import replace
import logging
from pathlib import Path
import sys
from time import perf_counter

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT / "src"))

from sat_string.config import dump_normalized_config, load_config
from sat_string.io import prepare_output_directory, write_json
from sat_string.simulator import CaseDefinition, simulate_case
from sat_string.switching_wave import (
    analyze_switching_wave,
    generate_aggregate_plots,
    generate_representative_fit_plot,
    save_sweep_csv,
    sweep_csv_row,
    threshold_scan_values,
)


def _logger(output: Path) -> logging.Logger:
    logger = logging.getLogger("sat_string.switching_wave_sweep")
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


def _nearest_representatives(
    requested: list[float], available: list[dict]
) -> list[dict]:
    selected: list[dict] = []
    used: set[float] = set()
    for target in requested:
        closest = min(available, key=lambda item: abs(float(item["h_on"]) - target))
        value = float(closest["h_on"])
        if value not in used:
            selected.append(closest)
            used.add(value)
    return selected


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument(
        "--center-satellite",
        type=int,
        help="One-based disturbance-center satellite; default: infer nearest to center_s_m",
    )
    parser.add_argument(
        "--representative-h-on",
        action="append",
        type=float,
        dest="representatives",
        help="Representative h_on for fit plot (repeatable; nearest grid value is used)",
    )
    parser.add_argument("--dpi", type=int, default=300)
    parser.add_argument("--overwrite", "-overwrite", action="store_true")
    args = parser.parse_args(argv)

    base_config = load_config(args.config)
    output = prepare_output_directory(args.output, overwrite=args.overwrite)
    dump_normalized_config(base_config, output / "base_normalized_config.yaml")
    logger = _logger(output)
    start = perf_counter()
    thresholds = threshold_scan_values()
    rows: list[dict] = []
    cases: list[dict] = []
    case_definition = CaseDefinition("controlled_disturbed", True, True)

    try:
        logger.info("Starting %d threshold cases", len(thresholds))
        for case_number, h_on in enumerate(thresholds, start=1):
            h_on = float(h_on)
            h_off = h_on - 1.0e-8
            config = replace(
                base_config,
                controller=replace(
                    base_config.controller,
                    threshold_on_m_s2=h_on,
                    threshold_off_m_s2=h_off,
                ),
            )
            case_start = perf_counter()
            logger.info(
                "Case %d/%d h_on=%.3e h_off=%.3e",
                case_number,
                len(thresholds),
                h_on,
                h_off,
            )
            result = simulate_case(config, case_definition)
            analysis = analyze_switching_wave(
                result.switch_events,
                config,
                center_satellite_id=args.center_satellite,
            )
            row = sweep_csv_row(h_on, h_off, analysis)
            rows.append(row)
            cases.append({"h_on": h_on, "h_off": h_off, "analysis": analysis})
            logger.info(
                "Completed h_on=%.3e delay=%.9g s n_fit=(%d,%d) in %.2f s",
                h_on,
                row["T_delay"],
                row["n_fit_left"],
                row["n_fit_right"],
                perf_counter() - case_start,
            )
            del result

        csv_path = save_sweep_csv(rows, output / "switching_wave_sweep.csv")
        plot_paths = generate_aggregate_plots(rows, output, dpi=args.dpi)
        requested = args.representatives or [1.0e-7, 1.4e-7, 2.0e-7]
        representatives = _nearest_representatives(requested, cases)
        fit_plot = generate_representative_fit_plot(
            representatives,
            output / "representative_distance_time_fits.png",
            dpi=args.dpi,
        )
        write_json(
            output / "sweep_metadata.json",
            {
                "threshold_count": len(thresholds),
                "threshold_on_values_m_s2": thresholds.tolist(),
                "threshold_off_rule": "h_off = h_on - 1.0e-8",
                "disturbance_time_s": base_config.disturbance.start_time_s,
                "center_satellite_id": cases[0]["analysis"]["center_satellite_id"],
                "side_definition": {
                    "left": "satellite index i < i_center",
                    "right": "satellite index i > i_center",
                },
                "fit_exclusion": "abs(i - i_center) <= 3",
                "fit_minimum_point_count": 3,
                "adjacent_time_std_ddof": 0,
                "representative_h_on_values_m_s2": [
                    case["h_on"] for case in representatives
                ],
                "artifacts": [
                    csv_path.name,
                    *(path.name for path in plot_paths),
                    fit_plot.name,
                ],
            },
        )
        logger.info("Sweep completed in %.2f s; results=%s", perf_counter() - start, csv_path)
        return 0
    except Exception:
        logger.exception("Switching-wave sweep failed; completed log remains in %s", output)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
