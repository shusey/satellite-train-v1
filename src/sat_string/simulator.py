"""Hybrid event-driven simulation orchestration."""

from __future__ import annotations

from collections.abc import Callable
from dataclasses import dataclass

import numpy as np
from numpy.typing import NDArray

from .attitude import AttitudeStateMachine, SwitchEvent
from .config import SimulationConfig
from .controller import (
    gap_kinematics,
    initial_formation,
    select_commands,
    switching_function,
)
from .dynamics import density_at_state, state_derivative
from .integrator import next_step_end, rk4_step


@dataclass(frozen=True)
class CaseDefinition:
    name: str
    control_enabled: bool
    disturbance_enabled: bool


@dataclass
class SimulationResult:
    case: CaseDefinition
    time_s: NDArray[np.float64]
    a_m: NDArray[np.float64]
    lambda_rad: NDArray[np.float64]
    altitude_m: NDArray[np.float64]
    area_m2: NDArray[np.float64]
    mode: NDArray[np.int8]
    density_kg_m3: NDArray[np.float64]
    q_m_s2: NDArray[np.float64]
    gap_error_m: NDArray[np.float64]
    gap_rate_m_s: NDArray[np.float64]
    gap_distance_m: NDArray[np.float64]
    switch_events: list[SwitchEvent]

    def arrays(self) -> dict[str, NDArray]:
        return {
            "time_s": self.time_s,
            "a_m": self.a_m,
            "lambda_rad": self.lambda_rad,
            "altitude_m": self.altitude_m,
            "area_m2": self.area_m2,
            "mode": self.mode,
            "density_kg_m3": self.density_kg_m3,
            "q_m_s2": self.q_m_s2,
            "gap_error_m": self.gap_error_m,
            "gap_rate_m_s": self.gap_rate_m_s,
            "gap_distance_m": self.gap_distance_m,
        }


STANDARD_CASES = (
    CaseDefinition("controlled_disturbed", True, True),
    CaseDefinition("controlled_undisturbed", True, False),
    CaseDefinition("all_low_drag_baseline", False, True),
)


def _formation_quantities(
    state: NDArray[np.float64], config: SimulationConfig
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    error, rate, distance = gap_kinematics(
        state[0],
        state[1],
        reference_radius_m=config.derived.reference_radius_m,
        target_spacing_rad=config.derived.target_spacing_rad,
        mu_m3_s2=config.constants.mu_m3_s2,
    )
    q = switching_function(
        error,
        rate,
        config.derived.controller_kp_s2_inv,
        config.derived.controller_kd_s_inv,
    )
    return error, rate, distance, q


def simulate_case(
    config: SimulationConfig,
    case: CaseDefinition,
    *,
    progress_callback: Callable[[float, float], None] | None = None,
    log_callback: Callable[[str], None] | None = None,
) -> SimulationResult:
    """Run one independent case using the required hybrid event ordering."""
    count = config.formation.satellite_count
    radius, longitude = initial_formation(
        count,
        config.derived.reference_radius_m,
        config.derived.target_spacing_rad,
    )
    state = np.vstack((radius, longitude))
    attitude = AttitudeStateMachine(
        count,
        config.spacecraft.area_low_m2,
        config.spacecraft.area_high_m2,
        config.controller.slew_time_s,
        config.controller.dwell_time_s,
        config.controller.command_delay_s,
    )

    stored: dict[str, list[NDArray | float]] = {
        key: []
        for key in (
            "time_s",
            "a_m",
            "lambda_rad",
            "altitude_m",
            "area_m2",
            "mode",
            "density_kg_m3",
            "q_m_s2",
            "gap_error_m",
            "gap_rate_m_s",
            "gap_distance_m",
        )
    }
    time_s = 0.0
    final_time = config.simulation.duration_s
    next_control = 0.0 if case.control_enabled else float("inf")
    next_output = 0.0
    tolerance = 1e-9
    next_progress_fraction = 0.05

    def log(message: str) -> None:
        if log_callback is not None:
            log_callback(message)

    def advance_schedule(current: float, period: float) -> float:
        candidate = current + period
        if candidate > final_time + tolerance:
            return float("inf")
        return min(candidate, final_time)

    while True:
        # 1. Start delayed commands, then 2. finalize completed transitions.
        for index in attitude.process_pending_starts(time_s, tolerance):
            log(f"t={time_s:.9g} s satellite={index + 1} transition started")
        for index in attitude.process_completions(time_s, tolerance):
            log(
                f"t={time_s:.9g} s satellite={index + 1} transition completed; "
                f"next command at t>={attitude.next_allowed_time[index]:.9g} s"
            )

        error, rate, distance, q = _formation_quantities(state, config)

        # 3-4. Compute from one snapshot and then apply all commands together.
        if abs(time_s - next_control) <= tolerance:
            targets = select_commands(
                q,
                attitude.mode.copy(),
                attitude.command_allowed(time_s, tolerance),
                config.controller.threshold_on_m_s2,
                config.controller.threshold_off_m_s2,
            )
            old_event_count = len(attitude.events)
            attitude.apply_commands(time_s, targets)
            for event in attitude.events[old_event_count:]:
                log(
                    f"t={event.command_time:.9g} s satellite={event.satellite_id} "
                    f"command {event.from_mode}->{event.to_mode}; "
                    f"start={event.start_time:.9g} complete={event.complete_time:.9g} s"
                )
            next_control = advance_schedule(
                next_control, config.controller.control_period_s
            )

        # 5. Save the post-event state at scheduled output times.
        if abs(time_s - next_output) <= tolerance:
            density = density_at_state(
                time_s,
                state[0],
                state[1],
                config,
                disturbance_enabled=case.disturbance_enabled,
            )
            stored["time_s"].append(float(time_s))
            stored["a_m"].append(state[0].copy())
            stored["lambda_rad"].append(state[1].copy())
            stored["altitude_m"].append(
                state[0].copy() - config.constants.earth_radius_m
            )
            stored["area_m2"].append(attitude.areas(time_s))
            stored["mode"].append(attitude.mode.copy())
            stored["density_kg_m3"].append(density)
            stored["q_m_s2"].append(q.copy())
            stored["gap_error_m"].append(error.copy())
            stored["gap_rate_m_s"].append(rate.copy())
            stored["gap_distance_m"].append(distance.copy())
            next_output = advance_schedule(
                next_output, config.simulation.output_step_s
            )

        if time_s >= final_time - tolerance:
            break

        # 6. No RK4 step may cross a control, output, delayed-start, or completion event.
        attitude_event = attitude.next_event_time(time_s, tolerance)
        step_end = next_step_end(
            time_s,
            config.simulation.integration_step_s,
            final_time,
            (next_control, next_output, attitude_event),
            tolerance,
        )
        if step_end <= time_s + tolerance:
            raise RuntimeError(
                f"Event scheduler failed to advance at t={time_s:.12g} s "
                f"(control={next_control}, output={next_output}, attitude={attitude_event})"
            )

        def rhs(stage_time: float, stage_state: NDArray[np.float64]) -> NDArray[np.float64]:
            return state_derivative(
                stage_time,
                stage_state,
                config,
                attitude.areas,
                disturbance_enabled=case.disturbance_enabled,
            )

        state = rk4_step(rhs, time_s, state, step_end - time_s)
        time_s = step_end
        if (
            np.any(state[0] <= config.constants.earth_radius_m)
            or not np.all(np.isfinite(state))
        ):
            raise FloatingPointError(
                f"Invalid orbital state after integration at t={time_s:.9g} s"
            )

        fraction = time_s / final_time
        if progress_callback is not None and fraction + tolerance >= next_progress_fraction:
            progress_callback(time_s, final_time)
            while next_progress_fraction <= fraction + tolerance:
                next_progress_fraction += 0.05

    arrays = {key: np.asarray(value) for key, value in stored.items()}
    return SimulationResult(
        case=case,
        switch_events=list(attitude.events),
        **arrays,
    )


def simulate_standard_suite(
    config: SimulationConfig,
    *,
    progress_callback: Callable[[str, float, float], None] | None = None,
    log_callback: Callable[[str, str], None] | None = None,
) -> dict[str, SimulationResult]:
    """Run all three standard cases independently."""
    results: dict[str, SimulationResult] = {}
    for case in STANDARD_CASES:
        progress = None
        logger = None
        if progress_callback is not None:
            progress = lambda time, final, name=case.name: progress_callback(name, time, final)
        if log_callback is not None:
            logger = lambda message, name=case.name: log_callback(name, message)
        results[case.name] = simulate_case(
            config,
            case,
            progress_callback=progress,
            log_callback=logger,
        )
    return results

