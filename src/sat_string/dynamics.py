"""Orbit-averaged near-circular differential-drag dynamics."""

from __future__ import annotations

from collections.abc import Callable

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .atmosphere import background_density
from .config import SimulationConfig
from .disturbance import relative_density_increment, rotating_coordinate_m


def mean_motion_rad_s(mu_m3_s2: float, radius_m: ArrayLike) -> NDArray[np.float64]:
    """Return sqrt(mu/a^3)."""
    radius = np.asarray(radius_m, dtype=float)
    return np.sqrt(mu_m3_s2 / radius**3)


def orbital_rates(
    radius_m: ArrayLike,
    density_kg_m3: ArrayLike,
    area_m2: ArrayLike,
    *,
    mu_m3_s2: float,
    drag_coefficient: float,
    mass_kg: float,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Return da/dt and unwrapped mean-longitude rate for each satellite."""
    radius = np.asarray(radius_m, dtype=float)
    density = np.asarray(density_kg_m3, dtype=float)
    area = np.asarray(area_m2, dtype=float)
    radius_rate = (
        -density
        * drag_coefficient
        * area
        / mass_kg
        * np.sqrt(mu_m3_s2 * radius)
    )
    return radius_rate, mean_motion_rad_s(mu_m3_s2, radius)


def density_at_state(
    time_s: float,
    radius_m: ArrayLike,
    lambda_rad: ArrayLike,
    config: SimulationConfig,
    *,
    disturbance_enabled: bool,
) -> NDArray[np.float64]:
    """Evaluate altitude and pulse-dependent density for an arbitrary RK stage."""
    radius = np.asarray(radius_m, dtype=float)
    altitude = radius - config.constants.earth_radius_m
    background = background_density(
        altitude,
        config.formation.reference_altitude_m,
        config.atmosphere.reference_density_kg_m3,
        config.atmosphere.effective_scale_height_m,
    )
    along_track = rotating_coordinate_m(
        lambda_rad,
        time_s,
        config.derived.reference_radius_m,
        config.derived.reference_mean_motion_rad_s,
    )
    disturbance = config.disturbance
    increment = relative_density_increment(
        along_track,
        time_s,
        enabled=disturbance_enabled and disturbance.enabled,
        amplitude_fraction=disturbance.amplitude_fraction,
        center_at_start_m=disturbance.center_s_m,
        start_time_s=disturbance.start_time_s,
        duration_s=disturbance.duration_s,
        sigma_m=config.derived.disturbance_sigma_m,
        propagation_speed_m_s=disturbance.propagation_speed_m_s,
    )
    return background * (1.0 + increment)


def state_derivative(
    time_s: float,
    state: NDArray[np.float64],
    config: SimulationConfig,
    area_provider: Callable[[float], ArrayLike],
    *,
    disturbance_enabled: bool,
) -> NDArray[np.float64]:
    """Evaluate [da/dt, dlambda/dt] with density and area at an RK stage."""
    state_array = np.asarray(state, dtype=float)
    if state_array.ndim != 2 or state_array.shape[0] != 2:
        raise ValueError("state must have shape (2, N)")
    radius, longitude = state_array
    if (
        np.any(radius <= config.constants.earth_radius_m)
        or not np.all(np.isfinite(state_array))
    ):
        raise FloatingPointError(
            f"Invalid orbital state at t={time_s:.9g} s: radius crossed Earth or contains NaN/Inf"
        )
    density = density_at_state(
        time_s,
        radius,
        longitude,
        config,
        disturbance_enabled=disturbance_enabled,
    )
    area = np.asarray(area_provider(time_s), dtype=float)
    radius_rate, longitude_rate = orbital_rates(
        radius,
        density,
        area,
        mu_m3_s2=config.constants.mu_m3_s2,
        drag_coefficient=config.spacecraft.drag_coefficient,
        mass_kg=config.spacecraft.mass_kg,
    )
    derivative = np.vstack((radius_rate, longitude_rate))
    if not np.all(np.isfinite(derivative)):
        raise FloatingPointError(f"Nonfinite derivative at t={time_s:.9g} s")
    return derivative

