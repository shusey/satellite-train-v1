from math import e

import numpy as np
import pytest

from sat_string.atmosphere import background_density
from sat_string.disturbance import (
    disturbance_center_m,
    fwhm_to_sigma,
    relative_density_increment,
    spatial_envelope,
    temporal_envelope,
)


def test_density_at_reference_altitude():
    rho0 = 3.584968e-12
    assert background_density(400e3, 400e3, rho0, 58_310.0) == pytest.approx(rho0)


def test_density_one_scale_height_higher():
    rho0 = 3.584968e-12
    scale = 58_310.0
    assert background_density(400e3 + scale, 400e3, rho0, scale) == pytest.approx(rho0 / e)


def test_temporal_envelope_endpoints_and_midpoint():
    values = temporal_envelope(np.array([10.0, 15.0, 20.0, 20.1]), 10.0, 10.0)
    assert values == pytest.approx([0.0, 1.0, 0.0, 0.0], abs=1e-15)


def test_spatial_fwhm_definition():
    fwhm = 100_000.0
    sigma = fwhm_to_sigma(fwhm)
    values = spatial_envelope(
        np.array([0.0, fwhm / 2.0]), 0.0, 0.0, 0.0, sigma, 0.0
    )
    assert values == pytest.approx([1.0, 0.5])


@pytest.mark.parametrize("speed", [-25.0, 40.0])
def test_moving_disturbance_center_and_direction(speed):
    times = np.array([100.0, 104.0, 111.0])
    centers = disturbance_center_m(times, 500.0, 100.0, speed)
    assert centers == pytest.approx(500.0 + speed * (times - 100.0))
    for time, center in zip(times, centers, strict=True):
        assert spatial_envelope(center, time, 500.0, 100.0, 1000.0, speed) == pytest.approx(1.0)


def test_disabled_disturbance_is_zero():
    result = relative_density_increment(
        np.array([-1.0, 0.0, 1.0]),
        5.0,
        enabled=False,
        amplitude_fraction=0.5,
        center_at_start_m=0.0,
        start_time_s=0.0,
        duration_s=10.0,
        sigma_m=1.0,
        propagation_speed_m_s=0.0,
    )
    assert np.array_equal(result, np.zeros(3))

