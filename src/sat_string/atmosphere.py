"""Altitude-dependent exponential background atmosphere."""

from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray


def background_density(
    altitude_m: ArrayLike,
    reference_altitude_m: float,
    reference_density_kg_m3: float,
    effective_scale_height_m: float,
) -> NDArray[np.float64]:
    """Return the local exponential density approximation in kg/m^3."""
    altitude = np.asarray(altitude_m, dtype=float)
    return reference_density_kg_m3 * np.exp(
        -(altitude - reference_altitude_m) / effective_scale_height_m
    )

