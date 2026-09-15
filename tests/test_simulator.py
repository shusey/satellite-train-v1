import numpy as np
import pytest

from sat_string.atmosphere import background_density
from sat_string.simulator import CaseDefinition, simulate_case, simulate_standard_suite


def test_all_low_undisturbed_preserves_zero_gap_error(short_config):
    result = simulate_case(
        short_config,
        CaseDefinition("all_low_undisturbed", False, False),
    )
    assert result.time_s[0] == 0.0
    assert result.time_s[-1] == short_config.simulation.duration_s
    assert np.max(np.abs(result.gap_error_m)) < 1e-6
    assert np.max(np.abs(result.gap_rate_m_s)) < 1e-12
    assert len(result.switch_events) == 0
    assert np.all(result.mode == 0)


def test_controlled_undisturbed_has_no_unnecessary_switches(short_config):
    result = simulate_case(
        short_config,
        CaseDefinition("controlled_undisturbed", True, False),
    )
    assert len(result.switch_events) == 0
    assert np.all(result.area_m2 == pytest.approx(short_config.spacecraft.area_low_m2))


def test_standard_pulse_peaks_at_center_and_midtime(short_config):
    result = simulate_case(
        short_config,
        CaseDefinition("disturbed", False, True),
    )
    background = background_density(
        result.altitude_m,
        short_config.formation.reference_altitude_m,
        short_config.atmosphere.reference_density_kg_m3,
        short_config.atmosphere.effective_scale_height_m,
    )
    relative_increment = result.density_kg_m3 / background - 1.0
    peak = np.unravel_index(np.argmax(relative_increment), relative_increment.shape)
    assert result.time_s[peak[0]] == pytest.approx(30.0)
    assert peak[1] == 2
    assert relative_increment[peak] == pytest.approx(
        short_config.disturbance.amplitude_fraction, rel=1e-6
    )


def test_noncommensurate_schedules_save_each_requested_time(short_config):
    result = simulate_case(
        short_config,
        CaseDefinition("schedule", True, False),
    )
    assert result.time_s == pytest.approx(np.arange(0.0, 60.0 + 5.0, 5.0))
    assert np.all(np.isfinite(result.a_m))
    assert np.all(result.a_m > short_config.constants.earth_radius_m)


def test_three_standard_cases_are_independent(short_config):
    suite = simulate_standard_suite(short_config)
    assert set(suite) == {
        "controlled_disturbed",
        "controlled_undisturbed",
        "all_low_drag_baseline",
    }
    assert suite["controlled_disturbed"].a_m is not suite["all_low_drag_baseline"].a_m
    assert np.max(
        np.abs(
            suite["controlled_disturbed"].density_kg_m3
            - suite["controlled_undisturbed"].density_kg_m3
        )
    ) > 0.0

