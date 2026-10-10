"""Offline assessment against hypothetical mission tolerances; no dynamics."""
from copy import deepcopy

import numpy as np


def exceedance_metrics(time_s, max_error_m, tolerance_m):
    """Strict > tolerance on saved samples.

    Durations use the left sample on [t[k], t[k+1]); the terminal sample has
    zero duration. First exceedance/ever also include the terminal observation.
    This is a sampled estimate, not a guarantee between samples. Nonuniform
    grids are weighted by actual interval lengths, never by sample counts.
    """
    t, y = np.asarray(time_s, float), np.asarray(max_error_m, float)
    if (t.ndim != 1 or len(t) < 2 or y.shape != t.shape
            or not np.all(np.isfinite(t)) or not np.all(np.isfinite(y))
            or np.any(np.diff(t) <= 0) or np.any(y < 0)
            or not np.isfinite(tolerance_m) or tolerance_m < 0):
        raise ValueError("Finite increasing times, nonnegative errors/tolerance required")
    over = y > tolerance_m
    ix = np.flatnonzero(over)
    edges = np.diff(np.r_[False, over[:-1], False].astype(int))
    starts, stops = np.flatnonzero(edges == 1), np.flatnonzero(edges == -1)
    durations = t[stops] - t[starts]
    total = float(np.sum(np.diff(t)[over[:-1]]))
    return dict(tolerance_m=float(tolerance_m), ever_exceeded=bool(ix.size),
                first_exceedance_time_s=float(t[ix[0]]) if ix.size else None,
                exceedance_duration_s=total, exceedance_fraction=total / (t[-1] - t[0]),
                longest_continuous_exceedance_s=float(durations.max()) if len(durations) else 0.,
                exceeds_at_observation_end=bool(over[-1]),
                maximum_sample_interval_s=float(np.diff(t).max()),
                duration_convention="left_sample_hold; strict >; final sample zero duration")


def validate_comparison_config(actual, reference, kp_factor):
    """Require precisely Kp alone to differ from the disturbed 60-day reference."""
    expected = deepcopy(reference)
    if kp_factor != 1:
        expected['controller']['kp_s2_inv'] = kp_factor * reference['derived']['controller_kp_s2_inv']
    expected['derived']['controller_kp_s2_inv'] *= kp_factor
    if actual != expected:
        raise ValueError("Saved input differs beyond the allowed Kp factor")


def improvement(controlled, reference):
    """Absolute and percent reduction, keeping both original magnitudes."""
    return dict(all_low_m=reference, controlled_m=controlled,
                reduction_m=reference-controlled,
                reduction_percent=100*(reference-controlled)/reference if reference else None)
