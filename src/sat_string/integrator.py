"""Fixed-step classical RK4 helpers with explicit event-boundary support."""

from __future__ import annotations

from collections.abc import Callable, Iterable

import numpy as np
from numpy.typing import NDArray


Rhs = Callable[[float, NDArray[np.float64]], NDArray[np.float64]]


def rk4_step(
    rhs: Rhs,
    time_s: float,
    state: NDArray[np.float64],
    step_s: float,
) -> NDArray[np.float64]:
    """Advance one fixed classical fourth-order Runge-Kutta step."""
    if step_s <= 0.0:
        raise ValueError("RK4 step must be positive")
    y = np.asarray(state, dtype=float)
    k1 = rhs(time_s, y)
    k2 = rhs(time_s + 0.5 * step_s, y + 0.5 * step_s * k1)
    k3 = rhs(time_s + 0.5 * step_s, y + 0.5 * step_s * k2)
    k4 = rhs(time_s + step_s, y + step_s * k3)
    result = y + (step_s / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)
    if not np.all(np.isfinite(result)):
        raise FloatingPointError(f"RK4 produced NaN/Inf at t={time_s + step_s:.9g} s")
    return result


def next_step_end(
    time_s: float,
    nominal_step_s: float,
    final_time_s: float,
    event_times_s: Iterable[float] = (),
    atol: float = 1e-9,
) -> float:
    """Clip a normal step at the earliest future discrete event."""
    if nominal_step_s <= 0.0:
        raise ValueError("nominal step must be positive")
    end = min(time_s + nominal_step_s, final_time_s)
    for event_time in event_times_s:
        if time_s + atol < event_time < end - atol:
            end = float(event_time)
        elif abs(event_time - end) <= atol:
            end = float(event_time)
    return end


def integrate_interval(
    rhs: Rhs,
    time_s: float,
    state: NDArray[np.float64],
    target_time_s: float,
    max_step_s: float,
) -> NDArray[np.float64]:
    """Integrate exactly to a caller-provided event boundary."""
    if target_time_s < time_s:
        raise ValueError("target time precedes current time")
    current = float(time_s)
    result = np.asarray(state, dtype=float)
    while current < target_time_s - 1e-12:
        step = min(max_step_s, target_time_s - current)
        result = rk4_step(rhs, current, result, step)
        current += step
    return result
