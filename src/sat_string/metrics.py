"""Metrics computed only from saved time histories and event records."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any, Mapping, Sequence

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .attitude import H, SwitchEvent
from .config import SimulationConfig
from .disturbance import disturbance_center_m, rotating_coordinate_m
from .io import write_json


def detect_order_reversal(
    time_s: ArrayLike, lambda_rad: ArrayLike
) -> dict[str, Any]:
    """Detect the first fixed-index longitude ordering reversal."""
    time = np.asarray(time_s, dtype=float)
    longitude = np.asarray(lambda_rad, dtype=float)
    reversed_mask = np.diff(longitude, axis=1) <= 0.0
    locations = np.argwhere(reversed_mask)
    if locations.size == 0:
        return {
            "order_reversal": False,
            "first_order_reversal_time_s": None,
            "first_order_reversal_gap_number": None,
        }
    first_row = locations[np.argmin(time[locations[:, 0]])]
    return {
        "order_reversal": True,
        "first_order_reversal_time_s": float(time[first_row[0]]),
        "first_order_reversal_gap_number": int(first_row[1] + 1),
    }


def recovery_metric(
    time_s: ArrayLike,
    gap_error_m: ArrayLike,
    gap_rate_m_s: ArrayLike,
    *,
    disturbance_end_time_s: float,
    position_tolerance_m: float,
    velocity_tolerance_m_s: float,
    hold_time_s: float,
) -> dict[str, Any]:
    """Find the first post-disturbance tolerance interval held continuously."""
    time = np.asarray(time_s, dtype=float)
    error = np.asarray(gap_error_m, dtype=float)
    rate = np.asarray(gap_rate_m_s, dtype=float)
    eligible = time >= disturbance_end_time_s - 1e-9
    within = (
        np.all(np.abs(error) <= position_tolerance_m, axis=1)
        & np.all(np.abs(rate) <= velocity_tolerance_m_s, axis=1)
        & eligible
    )
    run_start: float | None = None
    for sample_time, accepted in zip(time, within, strict=True):
        if accepted:
            if run_start is None:
                run_start = float(sample_time)
            if sample_time - run_start >= hold_time_s - 1e-9:
                return {
                    "recovered": True,
                    "recovery_interval_start_s": run_start,
                    "recovery_confirmation_time_s": float(sample_time),
                    "recovery_time_after_disturbance_s": float(
                        sample_time - disturbance_end_time_s
                    ),
                }
        else:
            run_start = None
    return {
        "recovered": False,
        "recovery_interval_start_s": None,
        "recovery_confirmation_time_s": None,
        "recovery_time_after_disturbance_s": None,
    }


def _first_major_peak_time(
    time_s: NDArray[np.float64], values: NDArray[np.float64]
) -> float | None:
    magnitude = np.abs(values)
    maximum = float(np.max(magnitude))
    if maximum <= 0.0:
        return None
    threshold = 0.5 * maximum
    if magnitude.size == 1:
        return float(time_s[0])
    candidates: list[int] = []
    for index in range(magnitude.size):
        left = magnitude[index - 1] if index > 0 else -np.inf
        right = magnitude[index + 1] if index + 1 < magnitude.size else -np.inf
        if magnitude[index] >= threshold and magnitude[index] >= left and magnitude[index] >= right:
            candidates.append(index)
    return float(time_s[candidates[0]]) if candidates else None


def _propagation_regression(
    positions_m: NDArray[np.float64], peak_times_s: Sequence[float | None]
) -> dict[str, Any]:
    valid = np.array([value is not None for value in peak_times_s], dtype=bool)
    if np.count_nonzero(valid) < 3:
        return {
            "status": "not_reliably_estimated",
            "reason": "fewer_than_3_points",
            "point_count": int(np.count_nonzero(valid)),
            "speed_m_s": None,
            "r_squared": None,
        }
    times = np.asarray([value for value in peak_times_s if value is not None], dtype=float)
    positions = positions_m[valid]
    design = np.column_stack((times, np.ones_like(times)))
    slope, intercept = np.linalg.lstsq(design, positions, rcond=None)[0]
    predicted = slope * times + intercept
    residual_sum = float(np.sum((positions - predicted) ** 2))
    total_sum = float(np.sum((positions - np.mean(positions)) ** 2))
    r_squared = 1.0 - residual_sum / total_sum if total_sum > 0.0 else 0.0
    if r_squared < 0.8:
        return {
            "status": "not_reliably_estimated",
            "reason": "r_squared_below_0.8",
            "point_count": int(times.size),
            "speed_m_s": None,
            "r_squared": r_squared,
        }
    return {
        "status": "estimated",
        "reason": None,
        "point_count": int(times.size),
        "speed_m_s": float(slope),
        "r_squared": r_squared,
    }


def _trend(values: NDArray[np.float64], distances: NDArray[np.float64]) -> dict[str, Any]:
    if values.size < 2:
        return {"classification": "insufficient_points", "slope_per_m": None}
    slope = float(np.polyfit(distances, values, 1)[0])
    scale = max(float(np.max(np.abs(values))), 1.0)
    tolerance = 1e-12 * scale / max(float(np.ptp(distances)), 1.0)
    classification = "flat"
    if slope > tolerance:
        classification = "increasing"
    elif slope < -tolerance:
        classification = "decreasing"
    return {"classification": classification, "slope_per_m": slope}


def _event_value(event: SwitchEvent | Mapping[str, Any], key: str) -> Any:
    return getattr(event, key) if isinstance(event, SwitchEvent) else event[key]


def compute_case_metrics(
    arrays: Mapping[str, NDArray],
    config: SimulationConfig,
    switch_events: Sequence[SwitchEvent | Mapping[str, Any]] = (),
) -> dict[str, Any]:
    """Compute all per-case v1 metrics from persisted histories."""
    time = np.asarray(arrays["time_s"], dtype=float)
    error = np.asarray(arrays["gap_error_m"], dtype=float)
    rate = np.asarray(arrays["gap_rate_m_s"], dtype=float)
    distance = np.asarray(arrays["gap_distance_m"], dtype=float)
    longitude = np.asarray(arrays["lambda_rad"], dtype=float)
    radii = np.asarray(arrays["a_m"], dtype=float)
    modes = np.asarray(arrays["mode"], dtype=np.int8)
    areas = np.asarray(arrays["area_m2"], dtype=float)
    maximum_error = np.max(np.abs(error), axis=0)
    minimum_gap_index = np.unravel_index(np.argmin(distance), distance.shape)

    start = config.disturbance.start_time_s
    end = start + config.disturbance.duration_s
    start_index = int(np.searchsorted(time, start, side="left"))
    start_index = min(start_index, len(time) - 1)
    satellite_positions = rotating_coordinate_m(
        longitude[start_index],
        time[start_index],
        config.derived.reference_radius_m,
        config.derived.reference_mean_motion_rad_s,
    )
    center = float(
        disturbance_center_m(
            time[start_index],
            config.disturbance.center_s_m,
            start,
            config.disturbance.propagation_speed_m_s,
        )
    )
    source_satellite = int(np.argmin(np.abs(satellite_positions - center)))
    source_gaps = [
        index
        for index in (source_satellite - 1, source_satellite)
        if 0 <= index < maximum_error.size
    ]
    source_error = max(float(maximum_error[index]) for index in source_gaps)
    amplification = maximum_error / max(source_error, 1e-12)
    gap_positions = 0.5 * (satellite_positions[:-1] + satellite_positions[1:])
    left_indices = np.flatnonzero(gap_positions < center)
    right_indices = np.flatnonzero(gap_positions >= center)
    left_distance = np.abs(gap_positions[left_indices] - center)
    right_distance = np.abs(gap_positions[right_indices] - center)

    after_start = time >= start - 1e-9
    peak_times: list[float | None] = []
    for gap_index in range(error.shape[1]):
        peak_times.append(
            _first_major_peak_time(time[after_start], error[after_start, gap_index])
        )

    transition_count = np.zeros(config.formation.satellite_count, dtype=int)
    for event in switch_events:
        if float(_event_value(event, "start_time")) <= time[-1] + 1e-9:
            transition_count[int(_event_value(event, "satellite_id")) - 1] += 1
    time_delta = np.diff(time)
    high_duration = np.sum((modes[:-1] == H) * time_delta[:, None], axis=0)
    above_low_duration = np.sum(
        (areas[:-1] > config.spacecraft.area_low_m2) * time_delta[:, None], axis=0
    )
    excess_area_integral = np.trapezoid(
        np.maximum(areas - config.spacecraft.area_low_m2, 0.0), time, axis=0
    )
    mean_radius = np.mean(radii, axis=1)
    relative_altitude = radii - mean_radius[:, None]
    altitude_spread = np.ptp(relative_altitude, axis=1)

    metrics: dict[str, Any] = {
        "maximum_gap_error_m_by_gap": maximum_error.tolist(),
        "maximum_gap_error_m": float(np.max(maximum_error)),
        "minimum_adjacent_spacing_m": float(distance[minimum_gap_index]),
        "minimum_adjacent_spacing_time_s": float(time[minimum_gap_index[0]]),
        "minimum_adjacent_spacing_gap_number": int(minimum_gap_index[1] + 1),
        **detect_order_reversal(time, longitude),
        "experimental_propagation_amplification": {
            "source_satellite_number": source_satellite + 1,
            "source_gap_numbers": [index + 1 for index in source_gaps],
            "source_maximum_error_m": source_error,
            "ratio_by_gap": amplification.tolist(),
            "left_outward_trend": _trend(amplification[left_indices], left_distance),
            "right_outward_trend": _trend(amplification[right_indices], right_distance),
        },
        "major_peak_time_s_by_gap": peak_times,
        "propagation_velocity": {
            "left": _propagation_regression(
                gap_positions[left_indices], [peak_times[index] for index in left_indices]
            ),
            "right": _propagation_regression(
                gap_positions[right_indices], [peak_times[index] for index in right_indices]
            ),
        },
        "recovery": {
            **recovery_metric(
                time,
                error,
                rate,
                disturbance_end_time_s=end,
                position_tolerance_m=config.recovery.position_tolerance_m,
                velocity_tolerance_m_s=config.recovery.velocity_tolerance_m_s,
                hold_time_s=config.recovery.hold_time_s,
            ),
            "time_resolution_s": config.simulation.output_step_s,
        },
        "attitude_switching": {
            "transition_start_count_by_satellite": transition_count.tolist(),
            "total_transition_start_count": int(np.sum(transition_count)),
            "high_stable_duration_s_by_satellite": high_duration.tolist(),
            "area_above_low_duration_s_by_satellite": above_low_duration.tolist(),
            "excess_area_time_integral_m2_s_by_satellite": excess_area_integral.tolist(),
        },
        "altitude": {
            "final_mean_altitude_m": float(
                mean_radius[-1] - config.constants.earth_radius_m
            ),
            "maximum_relative_altitude_spread_m": float(np.max(altitude_spread)),
            "final_relative_altitude_spread_m": float(altitude_spread[-1]),
        },
    }
    return metrics


def compute_suite_comparison(
    disturbed: Mapping[str, NDArray],
    undisturbed: Mapping[str, NDArray],
    low_drag_baseline: Mapping[str, NDArray],
) -> tuple[dict[str, Any], dict[str, NDArray[np.float64]]]:
    """Compute disturbance-only and differential-drag altitude comparisons."""
    for other in (undisturbed, low_drag_baseline):
        if not np.array_equal(disturbed["time_s"], other["time_s"]):
            raise ValueError("suite cases must share identical saved times")
    gap_effect = disturbed["gap_error_m"] - undisturbed["gap_error_m"]
    altitude_cost = disturbed["a_m"] - low_drag_baseline["a_m"]
    altitude_loss = -altitude_cost
    mean_radius = np.mean(disturbed["a_m"], axis=1)
    relative_altitude = disturbed["a_m"] - mean_radius[:, None]
    altitude_spread = np.ptp(relative_altitude, axis=1)
    comparison_arrays = {
        "time_s": np.asarray(disturbed["time_s"], dtype=float),
        "gap_error_disturbance_effect_m": gap_effect,
        "altitude_cost_m": altitude_cost,
        "altitude_loss_m": altitude_loss,
        "relative_altitude_m": relative_altitude,
        "altitude_spread_m": altitude_spread,
    }
    metrics = {
        "maximum_disturbance_only_gap_error_m_by_gap": np.max(
            np.abs(gap_effect), axis=0
        ).tolist(),
        "maximum_disturbance_only_gap_error_m": float(np.max(np.abs(gap_effect))),
        "maximum_altitude_loss_vs_all_low_m": float(np.max(altitude_loss)),
        "final_mean_altitude_loss_vs_all_low_m": float(
            np.mean(altitude_loss[-1])
        ),
        "maximum_relative_altitude_spread_m": float(np.max(altitude_spread)),
    }
    return metrics, comparison_arrays


def _flatten(prefix: str, value: Any, rows: list[tuple[str, str]]) -> None:
    if isinstance(value, Mapping):
        for key, item in value.items():
            _flatten(f"{prefix}.{key}" if prefix else str(key), item, rows)
    else:
        rows.append((prefix, json.dumps(value, ensure_ascii=False, allow_nan=False)))


def save_metrics(
    metrics: Mapping[str, Any], output_directory: str | Path
) -> tuple[Path, Path]:
    """Save metrics as structured JSON and a portable two-column CSV."""
    destination = Path(output_directory)
    destination.mkdir(parents=True, exist_ok=True)
    json_path = destination / "metrics.json"
    csv_path = destination / "metrics.csv"
    write_json(json_path, metrics)
    rows: list[tuple[str, str]] = []
    _flatten("", metrics, rows)
    with csv_path.open("w", encoding="utf-8", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(["metric", "value_json"])
        writer.writerows(rows)
    return json_path, csv_path
