"""Post-processing pipeline for a saved standard three-case suite."""

from __future__ import annotations

from pathlib import Path
from typing import Any

import numpy as np

from .config import SimulationConfig
from .io import load_case_arrays, load_switch_events
from .metrics import compute_case_metrics, compute_suite_comparison, save_metrics
from .plotting import generate_formation_animation, generate_standard_plots
from .simulator import STANDARD_CASES


def build_standard_report(
    run_directory: str | Path,
    config: SimulationConfig,
    *,
    include_animation: bool | None = None,
    dpi: int = 300,
) -> dict[str, Any]:
    """Recompute metrics and visuals without rerunning any simulation."""
    run = Path(run_directory)
    arrays: dict[str, dict[str, np.ndarray]] = {}
    case_metrics: dict[str, Any] = {}
    for case in STANDARD_CASES:
        case_directory = run / case.name
        arrays[case.name] = load_case_arrays(case_directory)
        events = load_switch_events(case_directory)
        case_metrics[case.name] = compute_case_metrics(
            arrays[case.name], config, events
        )
        save_metrics(case_metrics[case.name], case_directory)

    comparison_metrics, comparison_arrays = compute_suite_comparison(
        arrays["controlled_disturbed"],
        arrays["controlled_undisturbed"],
        arrays["all_low_drag_baseline"],
    )
    np.savez_compressed(run / "comparison_results.npz", **comparison_arrays)
    aggregate = {"cases": case_metrics, "comparison": comparison_metrics}
    save_metrics(aggregate, run)
    generate_standard_plots(
        arrays["controlled_disturbed"],
        arrays["controlled_undisturbed"],
        arrays["all_low_drag_baseline"],
        config,
        run,
        dpi=dpi,
    )
    animate = config.visualization.animation_enabled if include_animation is None else include_animation
    if animate:
        generate_formation_animation(
            arrays["controlled_disturbed"],
            config,
            run / "formation_animation.mp4",
        )
    return aggregate
