"""Postprocessing edge cases: delayed/incomplete slews and later re-exits."""
import numpy as np
import pytest
from pathlib import Path

from sat_string.config import load_config
from sat_string.evaluation import EvaluationSettings, compute_evaluation_metrics, settling_metrics, switching_metrics
from sat_string.simulator import simulate_standard_suite


def event(command, start, finish, rising=True, satellite=1):
    return dict(satellite_id=satellite, command_time=command, start_time=start,
                complete_time=finish, from_mode="L" if rising else "H", to_mode="H" if rising else "L")


def test_observed_counts_exclude_future_starts_and_completions():
    events = [event(1, 2, 4), event(6, 7, 10, False), event(9, 11, 13, satellite=2)]
    m = switching_metrics(events, 2, 9)
    assert (m["n_commands"], m["n_transition_starts"], m["n_transition_completions"]) == (3, 2, 1)
    assert m["last_switch_time_s"] == 4
    assert m["min_switch_interval_s"] == 5
    assert m["lhl_cycles_per_satellite"] == [0, 0]
    assert m["pending_or_active_transitions_at_tmax"] == 2


def test_incomplete_secondary_cycle_is_flagged_as_started_only():
    events = [event(0, 0, 2), event(5, 5, 7, False), event(10, 10, 12)]
    m = switching_metrics(events, 1, 11)
    assert m["lhl_cycles_per_satellite"] == [1]
    assert m["secondary_switching_detected"] is False
    assert m["secondary_response_started"] is True


@pytest.mark.parametrize("n", [0, 1, 2, 3])
def test_lhl_counts_completed_cycles_and_classifies(n):
    events = [event(10*k, 10*k, 10*k+2, rising=k % 2 == 0) for k in range(2*n)]
    m = switching_metrics(events, 1, 100)
    assert m["lhl_cycles_per_satellite"] == [n]
    assert m["secondary_switching_detected"] == (n >= 2)
    assert m["multiple_switching_detected"] == (n >= 3)


def test_recovery_precedes_later_exit_and_final_settling():
    t = np.arange(21.)
    e = np.zeros((21, 1))
    e[10] = 10
    switches = switching_metrics([event(13, 13, 15)], 1, 20)
    m = settling_metrics(t, e, np.zeros_like(e), switches, EvaluationSettings(1, .001, 3), evaluation_start_s=2)
    assert m["recovery_entry_time_s"] == 2
    assert m["recovery_time_s"] == 5
    assert m["settling_time_s"] == 15
    assert m["verification_horizon_s"] == 5


def test_pending_transition_prevents_settling():
    t = np.arange(11.)
    e = np.zeros((11, 1))
    switches = switching_metrics([event(9, 12, 14)], 1, 10)
    m = settling_metrics(t, e, e, switches, EvaluationSettings(1, .001, 2), evaluation_start_s=0)
    assert m["recovery_time_s"] == 2
    assert m["settled_by_tmax"] is False


def test_end_suffix_requires_minimum_verification_window():
    t = np.arange(11.)
    e = np.full((11, 1), 2.)
    e[-2:] = 0
    m = settling_metrics(t, e, np.zeros_like(e), switching_metrics([], 1, 10),
                         EvaluationSettings(1, .001, 2), evaluation_start_s=0)
    assert m["settling_time_s"] is None


def test_no_switches_are_null_not_zero():
    m = switching_metrics([], 2, 10)
    assert m["last_switch_time_s"] is None
    assert m["min_switch_interval_s"] is None
    assert m["switches_per_satellite"] == [0, 0]


def test_invalid_non_alternating_events_raise():
    with pytest.raises(ValueError, match="alternate"):
        switching_metrics([event(0, 0, 2), event(5, 5, 7)], 1, 10)


def test_sample_rms_and_reference_sign(short_config):
    suite = simulate_standard_suite(short_config)
    own = suite["controlled_disturbed"]
    low = suite["all_low_drag_baseline"].arrays()
    before = own.a_m.copy()
    m = compute_evaluation_metrics(own.arrays(), short_config, own.switch_events, low_baseline=low)
    assert m["formation"]["rms_error_m"] == pytest.approx(np.sqrt(np.mean(own.gap_error_m**2)))
    assert m["formation"]["rms_error_per_gap_m"] == pytest.approx(np.sqrt(np.mean(own.gap_error_m**2, axis=0)))
    assert m["altitude"]["extra_loss_vs_all_low_m"] == pytest.approx(np.mean(low["a_m"][-1] - own.a_m[-1]))
    assert np.array_equal(own.a_m, before)


def test_wave_has_separate_sides_and_tied_last_arrivals():
    config = load_config(Path(__file__).parent / "fixtures/regression_24h/config.yaml")
    events = [event(100 + abs(i - 11)*10, 100 + abs(i - 11)*10,
                    102 + abs(i - 11)*10, satellite=i) for i in range(1, 22) if i != 11]
    arrays = dict(time_s=np.array([0., 150., 300.]), gap_error_m=np.zeros((3, 20)),
                  gap_rate_m_s=np.zeros((3, 20)), a_m=np.full((3, 21), config.derived.reference_radius_m))
    m = compute_evaluation_metrics(arrays, config, events)
    wave = m["wave"]
    assert wave["left_speed_m_s"] == pytest.approx(10000)
    assert wave["right_speed_m_s"] == pytest.approx(10000)
    assert wave["max_hops"] == 10
    assert wave["last_arrival_satellite_ids"] == [1, 21]
    assert wave["left"]["last_arrival_satellite_ids"] == [1]
    assert wave["right"]["last_arrival_satellite_ids"] == [21]


@pytest.mark.parametrize("settings", [(-1, 1, 1), (1, float("nan"), 1), (1, 1, -1)])
def test_invalid_postprocessing_settings(settings):
    with pytest.raises(ValueError):
        EvaluationSettings(*settings)
