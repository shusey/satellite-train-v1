from pathlib import Path

import numpy as np
import pytest

from sat_string.config import load_config
from sat_string.controller import gap_kinematics, initial_formation
from sat_string.dynamics import mean_motion_rad_s, orbital_rates, state_derivative


CONFIG = load_config(Path(__file__).parents[1] / "configs" / "baseline.yaml")


def test_initial_formation_has_zero_gap_error_and_rate_for_odd_and_even_counts():
    for count in (20, 21):
        radius, longitude = initial_formation(
            count,
            CONFIG.derived.reference_radius_m,
            CONFIG.derived.target_spacing_rad,
        )
        error, rate, distance = gap_kinematics(
            radius,
            longitude,
            reference_radius_m=CONFIG.derived.reference_radius_m,
            target_spacing_rad=CONFIG.derived.target_spacing_rad,
            mu_m3_s2=CONFIG.constants.mu_m3_s2,
        )
        assert error == pytest.approx(np.zeros(count - 1), abs=1e-9)
        assert rate == pytest.approx(np.zeros(count - 1), abs=1e-15)
        assert distance == pytest.approx(
            np.full(count - 1, CONFIG.formation.target_spacing_m), abs=1e-9
        )


def test_common_decay_rates_are_identical():
    radius = np.full(3, CONFIG.derived.reference_radius_m)
    density = np.full(3, CONFIG.atmosphere.reference_density_kg_m3)
    area = np.full(3, CONFIG.spacecraft.area_low_m2)
    radius_rate, longitude_rate = orbital_rates(
        radius,
        density,
        area,
        mu_m3_s2=CONFIG.constants.mu_m3_s2,
        drag_coefficient=CONFIG.spacecraft.drag_coefficient,
        mass_kg=CONFIG.spacecraft.mass_kg,
    )
    assert np.ptp(radius_rate) == 0.0
    assert np.ptp(longitude_rate) == 0.0


def test_high_drag_decreases_radius_faster_and_increases_mean_motion():
    radius = np.array([CONFIG.derived.reference_radius_m] * 2)
    density = np.array([CONFIG.atmosphere.reference_density_kg_m3] * 2)
    area = np.array([CONFIG.spacecraft.area_low_m2, CONFIG.spacecraft.area_high_m2])
    radius_rate, _ = orbital_rates(
        radius,
        density,
        area,
        mu_m3_s2=CONFIG.constants.mu_m3_s2,
        drag_coefficient=CONFIG.spacecraft.drag_coefficient,
        mass_kg=CONFIG.spacecraft.mass_kg,
    )
    assert abs(radius_rate[1]) > abs(radius_rate[0])
    advanced_radius = radius + radius_rate * 100.0
    motion = mean_motion_rad_s(CONFIG.constants.mu_m3_s2, advanced_radius)
    assert advanced_radius[1] < advanced_radius[0]
    assert motion[1] > motion[0]


def test_invalid_orbit_state_stops_immediately():
    state = np.array(
        [[CONFIG.constants.earth_radius_m, CONFIG.derived.reference_radius_m], [0.0, 1.0]]
    )
    with pytest.raises(FloatingPointError, match="crossed Earth"):
        state_derivative(
            0.0,
            state,
            CONFIG,
            lambda _: np.full(2, CONFIG.spacecraft.area_low_m2),
            disturbance_enabled=False,
        )

