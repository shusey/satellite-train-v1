import json

import numpy as np
import pytest

from sat_string.metrics import (
    compute_case_metrics,
    compute_suite_comparison,
    detect_order_reversal,
    recovery_metric,
    save_metrics,
)
from sat_string.simulator import CaseDefinition, simulate_standard_suite


def test_recovery_requires_continuous_hold():
    time = np.arange(0.0, 10.0)
    error = np.zeros((10, 2))
    rate = np.zeros((10, 2))
    error[4, 0] = 100.0
    result = recovery_metric(
        time,
        error,
        rate,
        disturbance_end_time_s=2.0,
        position_tolerance_m=50.0,
        velocity_tolerance_m_s=1e-3,
        hold_time_s=3.0,
    )
    assert result["recovered"] is True
    assert result["recovery_interval_start_s"] == 5.0
    assert result["recovery_confirmation_time_s"] == 8.0
    assert result["recovery_time_after_disturbance_s"] == 6.0


def test_short_tolerance_entry_does_not_count_as_recovery():
    result = recovery_metric(
        np.array([0.0, 1.0, 2.0, 3.0]),
        np.array([[0.0], [0.0], [100.0], [0.0]]),
        np.zeros((4, 1)),
        disturbance_end_time_s=0.0,
        position_tolerance_m=50.0,
        velocity_tolerance_m_s=1e-3,
        hold_time_s=2.0,
    )
    assert result["recovered"] is False
    assert result["recovery_time_after_disturbance_s"] is None


def test_first_order_reversal_is_detected_without_reordering():
    longitude = np.array(
        [
            [0.0, 1.0, 2.0],
            [0.0, 1.1, 1.05],
            [0.2, 0.1, 1.2],
        ]
    )
    result = detect_order_reversal(np.array([0.0, 5.0, 10.0]), longitude)
    assert result == {
        "order_reversal": True,
        "first_order_reversal_time_s": 5.0,
        "first_order_reversal_gap_number": 2,
    }
    assert np.array_equal(longitude[1], [0.0, 1.1, 1.05])


def test_metrics_can_be_recomputed_from_arrays(short_config):
    result = simulate_standard_suite(short_config)["controlled_disturbed"]
    first = compute_case_metrics(result.arrays(), short_config, result.switch_events)
    second = compute_case_metrics(result.arrays(), short_config, result.switch_events)
    assert first == second
    assert len(first["maximum_gap_error_m_by_gap"]) == 4
    assert "recovery" in first
    assert "propagation_velocity" in first


def test_suite_comparison_shapes_and_zero_identity(short_config):
    suite = simulate_standard_suite(short_config)
    metrics, arrays = compute_suite_comparison(
        suite["controlled_disturbed"].arrays(),
        suite["controlled_undisturbed"].arrays(),
        suite["all_low_drag_baseline"].arrays(),
    )
    assert arrays["altitude_cost_m"].shape == suite["controlled_disturbed"].a_m.shape
    assert arrays["gap_error_disturbance_effect_m"].shape == suite["controlled_disturbed"].gap_error_m.shape
    assert metrics["maximum_disturbance_only_gap_error_m"] >= 0.0


def test_metrics_json_and_csv_are_written(tmp_path):
    metrics = {"scalar": 1.5, "nested": {"flag": False, "array": [1, 2]}}
    json_path, csv_path = save_metrics(metrics, tmp_path)
    assert json.loads(json_path.read_text(encoding="utf-8")) == metrics
    text = csv_path.read_text(encoding="utf-8")
    assert "nested.flag" in text
    assert "nested.array" in text
