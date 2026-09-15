"""Switching-wave analysis for threshold-sweep simulations."""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path
import tempfile
from typing import Any, Mapping, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .attitude import SwitchEvent
from .config import SimulationConfig


def threshold_scan_values() -> NDArray[np.float64]:
    """Return the requested de-duplicated threshold grid in ascending order."""
    # Construct the grid as integer multiples of 1e-9 to avoid floating-point
    # duplicates at 1.3e-7, 1.4e-7, and 1.5e-7.
    coarse = np.arange(100, 201, 10, dtype=np.int64)
    transition = np.arange(130, 151, 2, dtype=np.int64)
    return np.union1d(coarse, transition).astype(float) * 1.0e-9


def _event_value(event: SwitchEvent | Mapping[str, Any], key: str) -> Any:
    return getattr(event, key) if isinstance(event, SwitchEvent) else event[key]


def first_l_to_h_start_times(
    events: Sequence[SwitchEvent | Mapping[str, Any]],
    satellite_count: int,
    *,
    simulation_end_time_s: float | None = None,
) -> NDArray[np.float64]:
    """Return the first actually-started L-to-H time for every satellite."""
    first = np.full(satellite_count, np.nan, dtype=float)
    for event in events:
        if _event_value(event, "from_mode") != "L" or _event_value(event, "to_mode") != "H":
            continue
        start_time = float(_event_value(event, "start_time"))
        if simulation_end_time_s is not None and start_time > simulation_end_time_s + 1.0e-9:
            continue
        index = int(_event_value(event, "satellite_id")) - 1
        if not 0 <= index < satellite_count:
            raise ValueError(f"Satellite ID is outside 1..{satellite_count}: {index + 1}")
        if not np.isfinite(first[index]) or start_time < first[index]:
            first[index] = start_time
    return first


def infer_center_satellite_index(config: SimulationConfig) -> int:
    """Infer the zero-based satellite nearest the disturbance center at onset."""
    count = config.formation.satellite_count
    positions = (
        np.arange(count, dtype=float) - (count - 1) / 2.0
    ) * config.formation.target_spacing_m
    return int(np.argmin(np.abs(positions - config.disturbance.center_s_m)))


def fit_distance_time(
    time_s: ArrayLike, distance_m: ArrayLike, *, min_points: int = 3
) -> dict[str, float | int]:
    """Fit distance = speed * time + intercept and return fit diagnostics."""
    time = np.asarray(time_s, dtype=float)
    distance = np.asarray(distance_m, dtype=float)
    valid = np.isfinite(time) & np.isfinite(distance)
    time = time[valid]
    distance = distance[valid]
    point_count = int(time.size)
    if point_count < min_points:
        return {
            "speed_m_s": float("nan"),
            "intercept_m": float("nan"),
            "r_squared": float("nan"),
            "rmse_m": float("nan"),
            "point_count": point_count,
        }

    design = np.column_stack((time, np.ones_like(time)))
    speed, intercept = np.linalg.lstsq(design, distance, rcond=None)[0]
    predicted = speed * time + intercept
    residual = distance - predicted
    residual_sum = float(np.sum(residual**2))
    total_sum = float(np.sum((distance - np.mean(distance)) ** 2))
    r_squared = 1.0 - residual_sum / total_sum if total_sum > 0.0 else float("nan")
    return {
        "speed_m_s": float(speed),
        "intercept_m": float(intercept),
        "r_squared": r_squared,
        "rmse_m": float(np.sqrt(np.mean(residual**2))),
        "point_count": point_count,
    }


def _side_analysis(
    first_times_s: NDArray[np.float64],
    center_index: int,
    spacing_m: float,
    *,
    side: str,
) -> dict[str, Any]:
    if side == "left":
        indices = np.arange(center_index - 1, -1, -1, dtype=int)
    elif side == "right":
        indices = np.arange(center_index + 1, first_times_s.size, dtype=int)
    else:
        raise ValueError(f"Unknown side: {side}")

    # Both arrays are ordered from the disturbance center outwards.
    offsets = np.abs(indices - center_index)
    switched = np.isfinite(first_times_s[indices])
    switched_indices = indices[switched]
    switched_offsets = offsets[switched]
    switched_times = first_times_s[switched_indices]
    switched_distances = switched_offsets.astype(float) * spacing_m

    fit_mask = switched_offsets >= 4
    fit = fit_distance_time(
        switched_times[fit_mask], switched_distances[fit_mask], min_points=3
    )

    delta_records: list[dict[str, int | float]] = []
    for inner_position in range(len(switched_indices) - 1):
        # Never bridge across a satellite that did not switch: Delta t is defined
        # only for physically adjacent satellites.
        if switched_offsets[inner_position + 1] - switched_offsets[inner_position] != 1:
            continue
        delta_records.append(
            {
                "inner_satellite_id": int(switched_indices[inner_position] + 1),
                "outer_satellite_id": int(switched_indices[inner_position + 1] + 1),
                "delta_t_s": float(
                    switched_times[inner_position + 1] - switched_times[inner_position]
                ),
            }
        )
    delta_times = np.asarray([record["delta_t_s"] for record in delta_records], dtype=float)
    if delta_times.size:
        mean_delta = float(np.mean(delta_times))
        std_delta = float(np.std(delta_times, ddof=0))
        cv_delta = std_delta / mean_delta if mean_delta != 0.0 else float("nan")
    else:
        mean_delta = std_delta = cv_delta = float("nan")

    return {
        "side": side,
        "switched_satellite_ids": (switched_indices + 1).tolist(),
        "offsets_from_center": switched_offsets.tolist(),
        "distances_m": switched_distances.tolist(),
        "transition_start_times_s": switched_times.tolist(),
        "fit_satellite_ids": (switched_indices[fit_mask] + 1).tolist(),
        "fit": fit,
        "delta_t_records": delta_records,
        "delta_t_s": delta_times.tolist(),
        "mean_delta_t_s": mean_delta,
        "std_delta_t_s": std_delta,
        "cv_delta_t": cv_delta,
    }


def analyze_switching_wave(
    events: Sequence[SwitchEvent | Mapping[str, Any]],
    config: SimulationConfig,
    *,
    center_satellite_id: int | None = None,
) -> dict[str, Any]:
    """Analyze first L-to-H starts on both sides of the disturbance center."""
    count = config.formation.satellite_count
    if center_satellite_id is None:
        center_index = infer_center_satellite_index(config)
    else:
        if not 1 <= center_satellite_id <= count:
            raise ValueError(f"center_satellite_id must be in 1..{count}")
        center_index = center_satellite_id - 1

    first = first_l_to_h_start_times(
        events,
        count,
        simulation_end_time_s=config.simulation.duration_s,
    )
    finite = first[np.isfinite(first)]
    delay = (
        float(np.min(finite) - config.disturbance.start_time_s)
        if finite.size
        else float("nan")
    )
    return {
        "center_satellite_id": center_index + 1,
        "spacing_m": config.formation.target_spacing_m,
        "disturbance_time_s": config.disturbance.start_time_s,
        "first_transition_start_times_s": first.tolist(),
        "delay_s": delay,
        "left": _side_analysis(
            first, center_index, config.formation.target_spacing_m, side="left"
        ),
        "right": _side_analysis(
            first, center_index, config.formation.target_spacing_m, side="right"
        ),
    }


SWEEP_CSV_FIELDS = (
    "h_on",
    "h_off",
    "T_delay",
    "i_center",
    "v_left",
    "R2_left",
    "RMSE_left",
    "n_fit_left",
    "v_right",
    "R2_right",
    "RMSE_right",
    "n_fit_right",
    "delta_t_left",
    "delta_t_pairs_left",
    "mean_delta_t_left",
    "std_delta_t_left",
    "CV_left",
    "delta_t_right",
    "delta_t_pairs_right",
    "mean_delta_t_right",
    "std_delta_t_right",
    "CV_right",
    "n_switched_left",
    "n_switched_right",
)


def sweep_csv_row(h_on: float, h_off: float, analysis: Mapping[str, Any]) -> dict[str, Any]:
    """Flatten one switching-wave result into a one-row CSV representation."""
    left = analysis["left"]
    right = analysis["right"]
    return {
        "h_on": h_on,
        "h_off": h_off,
        "T_delay": analysis["delay_s"],
        "i_center": analysis["center_satellite_id"],
        "v_left": left["fit"]["speed_m_s"],
        "R2_left": left["fit"]["r_squared"],
        "RMSE_left": left["fit"]["rmse_m"],
        "n_fit_left": left["fit"]["point_count"],
        "v_right": right["fit"]["speed_m_s"],
        "R2_right": right["fit"]["r_squared"],
        "RMSE_right": right["fit"]["rmse_m"],
        "n_fit_right": right["fit"]["point_count"],
        "delta_t_left": json.dumps(left["delta_t_s"], separators=(",", ":")),
        "delta_t_pairs_left": json.dumps(left["delta_t_records"], separators=(",", ":")),
        "mean_delta_t_left": left["mean_delta_t_s"],
        "std_delta_t_left": left["std_delta_t_s"],
        "CV_left": left["cv_delta_t"],
        "delta_t_right": json.dumps(right["delta_t_s"], separators=(",", ":")),
        "delta_t_pairs_right": json.dumps(right["delta_t_records"], separators=(",", ":")),
        "mean_delta_t_right": right["mean_delta_t_s"],
        "std_delta_t_right": right["std_delta_t_s"],
        "CV_right": right["cv_delta_t"],
        "n_switched_left": len(left["switched_satellite_ids"]),
        "n_switched_right": len(right["switched_satellite_ids"]),
    }


def save_sweep_csv(rows: Sequence[Mapping[str, Any]], path: str | Path) -> Path:
    """Save one row per h_on, including individual adjacent propagation times."""
    destination = Path(path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    with destination.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=SWEEP_CSV_FIELDS)
        writer.writeheader()
        for row in rows:
            # Spell missing numerical estimates exactly as requested while retaining
            # numeric values in memory for plotting.
            serialized = {
                key: (
                    f"{float(value):.12g}"
                    if key in {"h_on", "h_off"}
                    else "NaN"
                    if isinstance(value, (float, np.floating)) and np.isnan(value)
                    else value
                )
                for key, value in row.items()
            }
            writer.writerow(serialized)
    return destination


def _matplotlib_pyplot():
    os.environ.setdefault(
        "MPLCONFIGDIR", str(Path(tempfile.gettempdir()) / "sat_string_matplotlib")
    )
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    return plt


def _save_line_plot(
    h_on: NDArray[np.float64],
    series: Sequence[tuple[str, NDArray[np.float64]]],
    *,
    ylabel: str,
    path: Path,
    dpi: int,
) -> Path:
    plt = _matplotlib_pyplot()
    figure, axis = plt.subplots(figsize=(8, 5))
    markers = ("o", "s", "^")
    linestyles = ("-", "--", ":")
    for index, (label, values) in enumerate(series):
        axis.plot(
            h_on,
            values,
            marker=markers[index % len(markers)],
            linestyle=linestyles[index % len(linestyles)],
            markerfacecolor="none" if index else None,
            label=label,
        )
    axis.set_xlabel(r"$h_{\rm on}$ (m s$^{-2}$)")
    axis.set_ylabel(ylabel)
    axis.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))
    if ylabel == r"$R^2$":
        axis.ticklabel_format(axis="y", style="plain", useOffset=False)
    if h_on.size > 1:
        margin = 0.02 * float(np.ptp(h_on))
        axis.set_xlim(float(np.min(h_on) - margin), float(np.max(h_on) + margin))
    axis.grid(True, alpha=0.3)
    if len(series) > 1:
        axis.legend()
    figure.savefig(path, dpi=dpi, bbox_inches="tight")
    plt.close(figure)
    return path


def generate_aggregate_plots(
    rows: Sequence[Mapping[str, Any]],
    output_directory: str | Path,
    *,
    dpi: int = 300,
) -> list[Path]:
    """Generate the five requested threshold-sweep aggregate figures."""
    output = Path(output_directory)
    output.mkdir(parents=True, exist_ok=True)
    h_on = np.asarray([row["h_on"] for row in rows], dtype=float)

    def values(key: str) -> NDArray[np.float64]:
        return np.asarray([row[key] for row in rows], dtype=float)

    specifications = (
        (
            "h_on_vs_delay.png",
            [("Delay", values("T_delay"))],
            r"Response delay $T_{\rm delay}$ (s)",
        ),
        (
            "h_on_vs_switching_speed.png",
            [("Left", values("v_left")), ("Right", values("v_right"))],
            r"Switching-wave speed $v_{\rm sw}$ (m s$^{-1}$)",
        ),
        (
            "h_on_vs_r_squared.png",
            [("Left", values("R2_left")), ("Right", values("R2_right"))],
            r"$R^2$",
        ),
        (
            "h_on_vs_rmse.png",
            [("Left", values("RMSE_left")), ("Right", values("RMSE_right"))],
            "Fit RMSE (m)",
        ),
        (
            "h_on_vs_mean_adjacent_time.png",
            [
                ("Left", values("mean_delta_t_left")),
                ("Right", values("mean_delta_t_right")),
            ],
            "Mean adjacent propagation time (s)",
        ),
    )
    return [
        _save_line_plot(h_on, series, ylabel=ylabel, path=output / name, dpi=dpi)
        for name, series, ylabel in specifications
    ]


def generate_representative_fit_plot(
    cases: Sequence[Mapping[str, Any]],
    output_path: str | Path,
    *,
    dpi: int = 300,
) -> Path:
    """Plot d_i versus first L-to-H start time and fitted lines for selected cases."""
    if not cases:
        raise ValueError("At least one representative case is required")
    plt = _matplotlib_pyplot()
    figure, axes = plt.subplots(
        len(cases), 2, figsize=(11, max(4.0, 3.5 * len(cases))), squeeze=False
    )
    for row_index, case in enumerate(cases):
        h_on = float(case["h_on"])
        analysis = case["analysis"]
        for column, side_name in enumerate(("left", "right")):
            axis = axes[row_index, column]
            side = analysis[side_name]
            times = np.asarray(side["transition_start_times_s"], dtype=float)
            distances = np.asarray(side["distances_m"], dtype=float)
            offsets = np.asarray(side["offsets_from_center"], dtype=int)
            excluded = offsets <= 3
            if np.any(excluded):
                axis.scatter(
                    times[excluded], distances[excluded] / 1000.0,
                    marker="x", color="0.5", label=r"Excluded ($|i-i_c|\leq3$)",
                )
            if np.any(~excluded):
                axis.scatter(
                    times[~excluded], distances[~excluded] / 1000.0,
                    marker="o", color="#2563eb", label="Fit candidates",
                )
            fit = side["fit"]
            if int(fit["point_count"]) >= 3 and times[~excluded].size:
                line_time = np.linspace(np.min(times[~excluded]), np.max(times[~excluded]), 100)
                line_distance = (
                    float(fit["speed_m_s"]) * line_time + float(fit["intercept_m"])
                )
                axis.plot(line_time, line_distance / 1000.0, color="#dc2626", label="Least-squares fit")
                annotation = (
                    f"v={fit['speed_m_s']:.3g} m/s\n"
                    f"R²={fit['r_squared']:.4f}, RMSE={fit['rmse_m']:.3g} m\n"
                    f"n={fit['point_count']}"
                )
            else:
                annotation = f"Fit unavailable (n={fit['point_count']} < 3)"
            axis.text(
                0.03, 0.97, annotation, transform=axis.transAxes,
                va="top", bbox={"facecolor": "white", "alpha": 0.8, "edgecolor": "0.8"},
            )
            axis.set_title(f"h_on={h_on:.3e}, {side_name}")
            axis.set_xlabel(r"First L$\to$H start time $t_i^\uparrow$ (s)")
            axis.set_ylabel(r"Distance $d_i$ (km)")
            axis.grid(True, alpha=0.3)
            if times.size:
                axis.legend(fontsize="small", loc="lower right")
    figure.tight_layout()
    destination = Path(output_path)
    destination.parent.mkdir(parents=True, exist_ok=True)
    figure.savefig(destination, dpi=dpi, bbox_inches="tight")
    plt.close(figure)
    return destination
