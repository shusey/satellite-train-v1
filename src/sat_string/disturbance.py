"""Single moving density-pulse envelopes in the reference rotating frame."""

from __future__ import annotations

from math import log, sqrt

import numpy as np
from numpy.typing import ArrayLike, NDArray


def fwhm_to_sigma(fwhm_m: float) -> float:
    """Convert a Gaussian full width at half maximum to standard deviation."""
    return fwhm_m / (2.0 * sqrt(2.0 * log(2.0)))


def rotating_coordinate_m(
    lambda_rad: ArrayLike,
    time_s: float | ArrayLike,
    reference_radius_m: float,
    reference_mean_motion_rad_s: float,
) -> NDArray[np.float64]:
    """Return along-track position in the unwrapped reference rotating frame."""
    return reference_radius_m * (
        np.asarray(lambda_rad, dtype=float)
        - reference_mean_motion_rad_s * np.asarray(time_s, dtype=float)
    )


def disturbance_center_m(
    time_s: float | ArrayLike,
    center_at_start_m: float,
    start_time_s: float,
    propagation_speed_m_s: float,
) -> NDArray[np.float64]:
    """Return s_c(t), allowing either sign of propagation speed."""
    return center_at_start_m + propagation_speed_m_s * (
        np.asarray(time_s, dtype=float) - start_time_s
    )


def temporal_envelope(
    time_s: float | ArrayLike,
    start_time_s: float,
    duration_s: float,
) -> NDArray[np.float64]:
    """Raised-cosine single-pulse time envelope."""
    time = np.asarray(time_s, dtype=float)
    phase = 2.0 * np.pi * (time - start_time_s) / duration_s
    inside = (time >= start_time_s) & (time <= start_time_s + duration_s)
    return np.where(inside, 0.5 * (1.0 - np.cos(phase)), 0.0)


def spatial_envelope(
    along_track_m: ArrayLike,
    time_s: float | ArrayLike,
    center_at_start_m: float,
    start_time_s: float,
    sigma_m: float,
    propagation_speed_m_s: float,
) -> NDArray[np.float64]:
    """Gaussian spatial envelope evaluated about the moving center."""
    along_track = np.asarray(along_track_m, dtype=float)
    center = disturbance_center_m(
        time_s, center_at_start_m, start_time_s, propagation_speed_m_s
    )
    return np.exp(-0.5 * ((along_track - center) / sigma_m) ** 2)


def relative_density_increment(
    along_track_m: ArrayLike,
    time_s: float,
    *,
    enabled: bool,
    amplitude_fraction: float,
    center_at_start_m: float,
    start_time_s: float,
    duration_s: float,
    sigma_m: float,
    propagation_speed_m_s: float,
) -> NDArray[np.float64]:
    """Return A_rho G_s G_t for each sampling position."""
    along_track = np.asarray(along_track_m, dtype=float)
    if not enabled:
        return np.zeros_like(along_track)
    return amplitude_fraction * spatial_envelope(
        along_track,
        time_s,
        center_at_start_m,
        start_time_s,
        sigma_m,
        propagation_speed_m_s,
    ) * temporal_envelope(time_s, start_time_s, duration_s)

