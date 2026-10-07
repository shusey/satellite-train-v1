"""Software regression against an immutable 24-hour reference, not a design claim."""
from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from .config import load_config
from .evaluation import compute_evaluation_metrics
from .io import load_case_arrays, load_metadata, load_switch_events
from .simulator import STANDARD_CASES


def compare_baseline(run_directory, reference_directory) -> dict:
    run, reference = Path(run_directory), Path(reference_directory)
    manifest = json.loads((reference / "reference.json").read_text(encoding="utf-8"))
    config = load_config(reference / "config.yaml")
    low = load_case_arrays(run / "all_low_drag_baseline")
    reports = {}
    for case in STANDARD_CASES:
        folder = run / case.name
        if load_config(folder / "normalized_config.yaml") != config:
            raise ValueError(f"Regression input changed: {case.name}")
        expected_case = {"name": case.name, "control_enabled": case.control_enabled,
                         "disturbance_enabled": case.disturbance_enabled}
        if load_metadata(folder)["case"] != expected_case:
            raise ValueError(f"Regression case flags changed: {case.name}")
        arrays = load_case_arrays(folder)
        original = load_case_arrays(reference / case.name)
        events = load_switch_events(folder)
        expected_events = load_switch_events(reference / case.name)
        timeline = {}
        for key, expected in original.items():
            actual = arrays[key]
            tol = manifest["array_tolerances"][key]
            matches = actual.shape == expected.shape and np.allclose(
                actual, expected, atol=tol["atol"], rtol=tol["rtol"])
            timeline[key] = {
                "passed": bool(matches),
                "max_abs_difference": float(np.max(np.abs(actual.astype(float) - expected.astype(float))))
                                      if actual.shape == expected.shape else None,
                **tol,
            }
        event_match = len(events) == len(expected_events)
        if event_match:
            for actual, expected in zip(events, expected_events, strict=True):
                event_match &= all(actual[key] == expected[key] for key in ("satellite_id", "from_mode", "to_mode"))
                event_match &= all(abs(actual[key] - expected[key]) <= 1e-6
                                   for key in ("command_time", "start_time", "complete_time"))
        m = compute_evaluation_metrics(arrays, config, events,
                                      low_baseline=low if case.disturbance_enabled else None)
        observed = {
            "max_abs_error_m": m["formation"]["max_abs_error_m"],
            "rms_error_m": m["formation"]["rms_error_m"],
            "n_transition_starts": m["switching"]["n_transition_starts"],
            "extra_loss_vs_all_low_m": m["altitude"]["extra_loss_vs_all_low_m"],
            "left_speed_m_s": m["wave"]["left_speed_m_s"],
            "right_speed_m_s": m["wave"]["right_speed_m_s"],
        }
        aggregate = {}
        for key, value in manifest["aggregates"][case.name].items():
            current = observed[key]
            aggregate[key] = {
                "expected": value, "actual": current,
                "passed": current is None if value is None else current is not None and bool(np.isclose(current, value, atol=1e-4, rtol=1e-7)),
            }
        passed = all(v["passed"] for v in timeline.values()) and event_match and all(v["passed"] for v in aggregate.values())
        reports[case.name] = {"passed": passed, "aggregate": aggregate,
                              "events_match": bool(event_match), "histories": timeline}
    return {"passed": all(c["passed"] for c in reports.values()),
            "purpose": "software regression only; gains and thresholds are provisional",
            "reference_source_commit": manifest["source_commit"], "cases": reports}
