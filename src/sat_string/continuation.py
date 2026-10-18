"""Resume a saved post-event snapshot using unchanged core primitives.

No force, control law, actuator interpolation, or RK4 implementation is copied.
The small scheduling loop follows simulator.simulate_case event ordering.
Only regular output/control boundaries can be resumed, avoiding a step split
that a longer uninterrupted run would not have made.
"""
from __future__ import annotations

import numpy as np

from .attitude import AttitudeStateMachine, H, H_TO_L, L, L_TO_H, SwitchEvent
from .dynamics import density_at_state, state_derivative
from .controller import select_commands
from .integrator import next_step_end, rk4_step
from .simulator import SimulationResult, _formation_quantities


def restore_attitude(config, events, end_s):
    attitude = AttitudeStateMachine(config.formation.satellite_count,
        config.spacecraft.area_low_m2, config.spacecraft.area_high_m2,
        config.controller.slew_time_s, config.controller.dwell_time_s,
        config.controller.command_delay_s,
        dwell_high_time_s=config.controller.dwell_high_time_s,
        dwell_low_time_s=config.controller.dwell_low_time_s)
    for value in sorted(events, key=lambda e: e.command_time if isinstance(e, SwitchEvent) else e["command_time"]):
        event = value if isinstance(value, SwitchEvent) else SwitchEvent(**value)
        if event.command_time > end_s + 1e-9:
            continue
        i = event.satellite_id - 1
        target = H if event.to_mode == "H" else L
        attitude.events.append(event)
        if event.complete_time <= end_s + 1e-9:
            attitude.mode[i] = target
            attitude.next_allowed_time[i] = event.complete_time + (config.controller.effective_dwell_high_s
                if target == H else config.controller.effective_dwell_low_s)
            attitude.transition_start_time[i] = np.nan
            attitude.transition_complete_time[i] = np.nan
        elif event.start_time <= end_s + 1e-9:
            attitude.mode[i] = L_TO_H if target == H else H_TO_L
            attitude.transition_start_time[i] = event.start_time
            attitude.transition_complete_time[i] = event.complete_time
        else:
            attitude.pending_target[i] = target
            attitude.pending_start_time[i] = event.start_time
    return attitude


def continue_case(arrays, source_config, config, case, events, *, progress_callback=None, density_provider=None):
    """Return joined zero-origin history; source arrays are never modified."""
    source, destination = source_config.to_dict(), config.to_dict()
    for data in (source, destination):
        data["simulation"].pop("duration_s")
    if source != destination:
        raise ValueError("Continuation may change only simulation.duration_s")
    time = float(arrays["time_s"][-1])
    if not np.isclose(time, source_config.simulation.duration_s, atol=1e-9, rtol=0):
        raise ValueError("Source history must reach its configured endpoint")
    if config.simulation.duration_s <= time:
        raise ValueError("Continuation endpoint must be later than source")
    periods = [config.simulation.output_step_s]
    if case.control_enabled:
        periods.append(config.controller.control_period_s)
    if any(not np.isclose(time / period, round(time / period), atol=1e-8, rtol=0) for period in periods):
        raise ValueError("Resume endpoint must be a regular output/control boundary")
    state = np.vstack((arrays["a_m"][-1], arrays["lambda_rad"][-1])).copy()
    attitude = restore_attitude(config, events, time)
    if not np.array_equal(attitude.mode, arrays["mode"][-1]) or not np.allclose(
            attitude.areas(time), arrays["area_m2"][-1], atol=1e-14, rtol=0):
        raise ValueError("Saved actuator state and event history disagree")
    final = config.simulation.duration_s
    next_control = time + config.controller.control_period_s if case.control_enabled else float("inf")
    next_output = time + config.simulation.output_step_s
    initial_time = time
    stored = {key: [] for key in arrays}
    next_fraction = .05

    def advance(current, period):
        candidate = current + period
        return float("inf") if candidate > final + 1e-9 else min(candidate, final)

    while time < final - 1e-9:
        end = next_step_end(time, config.simulation.integration_step_s, final,
                           (next_control, next_output, attitude.next_event_time(time, 1e-9)), 1e-9)
        if end <= time + 1e-9:
            raise RuntimeError("Continuation scheduler did not advance")
        state = rk4_step(lambda t, x: state_derivative(t, x, config, attitude.areas,
                         disturbance_enabled=case.disturbance_enabled, density_provider=density_provider), time, state, end-time)
        time = end
        if np.any(state[0] <= config.constants.earth_radius_m) or not np.all(np.isfinite(state)):
            raise FloatingPointError("Invalid continued orbital state")
        attitude.process_pending_starts(time, 1e-9)
        attitude.process_completions(time, 1e-9)
        error, rate, distance, q = _formation_quantities(state, config)
        if abs(time-next_control) <= 1e-9:
            targets = select_commands(q, attitude.mode.copy(), attitude.command_allowed(time, 1e-9),
                                      config.controller.threshold_on_m_s2, config.controller.threshold_off_m_s2)
            attitude.apply_commands(time, targets)
            next_control = advance(next_control, config.controller.control_period_s)
        if abs(time-next_output) <= 1e-9:
            snapshot = dict(time_s=float(time), a_m=state[0].copy(), lambda_rad=state[1].copy(),
                            altitude_m=state[0]-config.constants.earth_radius_m, area_m2=attitude.areas(time),
                            mode=attitude.mode.copy(), density_kg_m3=density_at_state(time,state[0],state[1],config,
                            disturbance_enabled=case.disturbance_enabled, density_provider=density_provider), q_m_s2=q.copy(),
                            gap_error_m=error.copy(), gap_rate_m_s=rate.copy(), gap_distance_m=distance.copy())
            for key, value in snapshot.items():
                stored[key].append(value)
            next_output = advance(next_output, config.simulation.output_step_s)
        fraction = (time-initial_time)/(final-initial_time)
        if progress_callback is not None and fraction >= next_fraction - 1e-9:
            progress_callback(time, final)
            next_fraction += .05
    joined = {key: np.concatenate((value, np.asarray(stored[key], dtype=value.dtype)), axis=0)
              for key, value in arrays.items()}
    return SimulationResult(case=case, switch_events=list(attitude.events), **joined)
