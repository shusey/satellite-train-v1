"""Formation kinematics and simultaneous binary-drag command selection."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray

from .dynamics import mean_motion_rad_s


LOW_MODE = 0
HIGH_MODE = 2
NO_COMMAND = -1


def initial_formation(
    satellite_count: int,
    reference_radius_m: float,
    target_spacing_rad: float,
) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    """Return symmetric initial radii and unwrapped longitudes for odd or even N."""
    indices = np.arange(satellite_count, dtype=float)
    longitude = (indices - (satellite_count - 1) / 2.0) * target_spacing_rad
    radius = np.full(satellite_count, reference_radius_m, dtype=float)
    return radius, longitude


def gap_kinematics(
    radius_m: ArrayLike,
    lambda_rad: ArrayLike,
    *,
    reference_radius_m: float,
    target_spacing_rad: float,
    mu_m3_s2: float,
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    """Return gap error, gap rate, and physical gap using the reference radius."""
    radius = np.asarray(radius_m, dtype=float)
    longitude = np.asarray(lambda_rad, dtype=float)
    longitude_difference = np.diff(longitude)
    error = reference_radius_m * (longitude_difference - target_spacing_rad)
    rate = reference_radius_m * np.diff(mean_motion_rad_s(mu_m3_s2, radius))
    distance = reference_radius_m * longitude_difference
    return error, rate, distance


def switching_function(
    gap_error_m: ArrayLike,
    gap_rate_m_s: ArrayLike,
    kp_s2_inv: float,
    kd_s_inv: float,
) -> NDArray[np.float64]:
    """Compute q for all satellites on a finite open path."""
    error = np.asarray(gap_error_m, dtype=float)
    rate = np.asarray(gap_rate_m_s, dtype=float)
    if error.ndim != 1 or rate.shape != error.shape or error.size < 1:
        raise ValueError("gap error and rate must be same-shaped nonempty vectors")
    edge_signal = kp_s2_inv * error + kd_s_inv * rate
    q = np.empty(error.size + 1, dtype=float)
    q[0] = edge_signal[0]
    q[-1] = -edge_signal[-1]
    if error.size > 1:
        q[1:-1] = edge_signal[1:] - edge_signal[:-1]
    return q


def select_commands(
    q_m_s2: ArrayLike,
    modes: ArrayLike,
    command_allowed: ArrayLike,
    threshold_on_m_s2: float,
    threshold_off_m_s2: float,
) -> NDArray[np.int8]:
    """Select all targets from one immutable pre-transition state snapshot."""
    q = np.asarray(q_m_s2, dtype=float)
    mode = np.asarray(modes, dtype=np.int8)
    allowed = np.asarray(command_allowed, dtype=bool)
    if q.shape != mode.shape or q.shape != allowed.shape:
        raise ValueError("q, modes, and command_allowed must have the same shape")
    commands = np.full(q.shape, NO_COMMAND, dtype=np.int8)
    commands[allowed & (mode == LOW_MODE) & (q >= threshold_on_m_s2)] = HIGH_MODE
    commands[allowed & (mode == HIGH_MODE) & (q <= threshold_off_m_s2)] = LOW_MODE
    return commands
