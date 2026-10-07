#!/usr/bin/env python3
"""Evaluate a saved history into a NEW directory; no simulation or overwrite."""
from __future__ import annotations

import argparse
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from sat_string.config import load_config
from sat_string.evaluation import EvaluationSettings, compute_evaluation_metrics
from sat_string.io import load_case_arrays, load_metadata, load_switch_events, prepare_output_directory, write_json
from sat_string.metrics import compute_case_metrics, save_metrics


def evaluate_saved(source, output, *, low=None, epsilon_e=None, epsilon_v=None, window=None,
                   horizon_s=None):
    config = load_config(source / "normalized_config.yaml")
    arrays = load_case_arrays(source)
    events = load_switch_events(source)
    if horizon_s is not None:
        if horizon_s > arrays["time_s"][-1] or horizon_s <= 0:
            raise ValueError("Requested horizon must be within the stored history")
        keep = arrays["time_s"] <= horizon_s + 1e-9
        arrays = {key: value[keep] for key, value in arrays.items()}
    reference = None
    if low is not None:
        low_config = load_config(low / "normalized_config.yaml").to_dict()
        own_config = config.to_dict()
        # Horizon and postprocessing/visualization options do not change physics.
        for data in (low_config, own_config):
            data["simulation"].pop("duration_s")
            data.pop("recovery")
            data.pop("visualization")
        if own_config != low_config or load_metadata(low)["case"]["control_enabled"]:
            raise ValueError("Reference must be a matching all-LOW simulation")
        if load_metadata(source)["case"]["disturbance_enabled"] != load_metadata(low)["case"]["disturbance_enabled"]:
            raise ValueError("Reference and case must have matching disturbance enable flags")
        reference = load_case_arrays(low)
        keep = reference["time_s"] <= arrays["time_s"][-1] + 1e-9
        reference = {key: value[keep] for key, value in reference.items()}
    eps_values = epsilon_e or [config.recovery.position_tolerance_m]
    velocity = config.recovery.velocity_tolerance_m_s if epsilon_v is None else epsilon_v
    hold = config.recovery.hold_time_s if window is None else window
    evaluations = []
    for eps in eps_values:
        settings = EvaluationSettings(eps, velocity, hold)
        evaluations.append(compute_evaluation_metrics(arrays, config, events,
                           settings=settings, low_baseline=reference))
    # Validate inputs before creating the destination; existing inputs stay read-only.
    destination = prepare_output_directory(output)
    metrics = compute_case_metrics(arrays, config, events,
                                  evaluation_settings=EvaluationSettings(eps_values[0], velocity, hold))
    metrics["altitude"].update(evaluations[0]["altitude"])
    save_metrics(metrics, destination)
    write_json(destination / "settling_sensitivity.json", {"evaluations": evaluations})
    write_json(destination / "provenance.json", {
        "source_directory": str(source.resolve()),
        "all_low_reference_directory": str(low.resolve()) if low else None,
        "case": load_metadata(source)["case"],
        "observation_end_time_s": float(arrays["time_s"][-1]),
        "configuration": config.to_dict(),
        "simulation_rerun": False,
    })
    return metrics


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--all-low", type=Path)
    parser.add_argument("--epsilon-e", type=float, nargs="+")
    parser.add_argument("--epsilon-v", type=float)
    parser.add_argument("--hold-window", type=float)
    parser.add_argument("--horizon-hours", type=float)
    args = parser.parse_args()
    m = evaluate_saved(args.run, args.output, low=args.all_low, epsilon_e=args.epsilon_e,
                      epsilon_v=args.epsilon_v, window=args.hold_window,
                      horizon_s=args.horizon_hours * 3600 if args.horizon_hours else None)
    print({key: m[key] for key in ("formation", "settling")})


if __name__ == "__main__":
    main()
