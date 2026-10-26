"""Profile-aware relay experiments using unchanged dynamics/actuator primitives.

Each batch member has its own open-path controller and attitude machine. Only
density evaluation is vectorized across members. B multipliers persist through
LOW, HIGH and slews. Large histories are local memory maps, not duplicated NPZ.
"""
from pathlib import Path
from time import perf_counter
import numpy as np
from .attitude import AttitudeStateMachine
from .continuation import restore_attitude
from .controller import gap_kinematics, select_commands, switching_function
from .dynamics import density_at_state, orbital_rates, node_acceleration_radius_rate
from .integrator import next_step_end, rk4_step
from .io import prepare_output_directory, write_json
from .control_cost import clean_json
from .passive_dispersion import initial_profile_state


def profile_rhs(t, state, config, provider, factors, areas, disturbance_enabled=False,
                node_acceleration=None):
    """The original orbital_rates with fixed satellite-specific Cd*A/m factors."""
    a, u = state
    if not np.all(np.isfinite(state)) or np.any(a <= config.constants.earth_radius_m):
        raise FloatingPointError('Invalid profile orbital state')
    rho = density_at_state(t, a, u, config, disturbance_enabled=disturbance_enabled,
                           density_provider=provider)
    da, du = orbital_rates(a, rho, areas(t), mu_m3_s2=config.constants.mu_m3_s2,
        drag_coefficient=config.spacecraft.drag_coefficient*factors,
        mass_kg=config.spacecraft.mass_kg)
    if node_acceleration is not None:
        da = da + node_acceleration_radius_rate(a, node_acceleration(t),
            reference_radius_m=config.derived.reference_radius_m,
            mu_m3_s2=config.constants.mu_m3_s2)
    return np.vstack((da, du))


def integrate_control_profiles(config, profiles, provider, output, *, checkpoint=None,
                               disturbance_enabled=False, progress=None,
                               node_acceleration=None, node_input_metadata=None):
    """Same pre-event/control/post-event-output ordering as simulate_case.

    Resume a post-event regular output/control boundary. Prefix is never changed;
    the suffix includes a duplicate boundary sample which loaders drop. Config,
    profiles and disturbance flags are checked before accepting a checkpoint.
    """
    n = config.formation.satellite_count
    if (node_acceleration is None) != (node_input_metadata is None):
        raise ValueError('Node input and reproducible metadata must be supplied together')
    if checkpoint is not None and checkpoint.get('node_input') != node_input_metadata:
        raise ValueError('Checkpoint node input differs')
    state, factors = initial_profile_state(config, [p['delta_a_m'] for p in profiles],
                                          [p['B_multiplier'] for p in profiles])
    end = config.simulation.duration_s
    save = config.simulation.output_step_s
    period = config.controller.control_period_s
    start = 0. if checkpoint is None else float(checkpoint['time_s'])
    if end <= start or end/save != round(end/save) or start/save != round(start/save):
        raise ValueError('Continuation must use regular output boundaries')
    if checkpoint is not None:
        old, new = checkpoint['configuration'].copy(), config.to_dict()
        old['simulation'] = dict(old['simulation']); old['simulation'].pop('duration_s')
        new['simulation'].pop('duration_s')
        if old != new or checkpoint['profiles'] != profiles or checkpoint['disturbance_enabled'] != disturbance_enabled:
            raise ValueError('Checkpoint inputs differ')
        if start/period != round(start/period):
            raise ValueError('Resume must be a regular control boundary')
        state = np.asarray(checkpoint['state'], float)
        if state.shape != (2, len(profiles)*n) or not np.all(np.isfinite(state)):
            raise ValueError('Invalid checkpoint state')
    attitudes = []
    for i in range(len(profiles)):
        if checkpoint is None:
            c = config.controller; s = config.spacecraft
            attitude = AttitudeStateMachine(n, s.area_low_m2, s.area_high_m2,
                c.slew_time_s, c.dwell_time_s, c.command_delay_s,
                dwell_high_time_s=c.dwell_high_time_s, dwell_low_time_s=c.dwell_low_time_s)
        else:
            attitude = restore_attitude(config, checkpoint['events'][i], start)
            if attitude.mode.tolist() != checkpoint['modes'][i] or not np.allclose(
                    attitude.areas(start), checkpoint['areas'][i], atol=1e-14, rtol=0):
                raise ValueError('Checkpoint actuator state disagrees with events')
        attitudes.append(attitude)
    out = prepare_output_directory(output)
    times = np.arange(start, end+save/2, save)
    np.save(out/'time_s.npy', times)
    shape = (len(times), len(profiles), n)
    arrays = {key: np.lib.format.open_memmap(out/f'{key}.npy', mode='w+',
              dtype='int8' if key == 'mode' else 'float64', shape=shape)
              for key in ('a_m', 'lambda_rad', 'density_kg_m3', 'area_m2', 'mode')}
    write_json(out/'metadata.json', dict(configuration=config.to_dict(), profiles=profiles,
               start_s=start, end_s=end, disturbance_enabled=disturbance_enabled,
               **({'node_input': node_input_metadata} if node_input_metadata is not None else {})))
    stable_areas = np.concatenate([a.areas(start) for a in attitudes])
    active_slews = any(np.any((a.mode==1)|(a.mode==3)) for a in attitudes)
    def areas(t):
        # Stable attitudes have time-independent area; avoid repeated no-op array
        # construction at four RK stages. During slews use the ORIGINAL actuator.
        return np.concatenate([a.areas(t) for a in attitudes]) if active_slews else stable_areas
    rhs = lambda t, x: profile_rhs(t, x, config, provider, factors, areas, disturbance_enabled,
                                  node_acceleration)
    t = start; k = 0; began = perf_counter(); next_progress = start+86400.
    next_control = start+(period if checkpoint is not None else 0.)
    if not any(p['control_enabled'] for p in profiles): next_control = float('inf')
    next_output = start
    minimum_gap = float('inf')
    next_actuator = min(a.next_event_time(t) for a in attitudes)
    while True:
        actuator_due = abs(t-next_actuator)<=1e-9
        control_due = abs(t-next_control)<=1e-9
        update = actuator_due or control_due
        if update and (checkpoint is None or t > start+1e-9):
            for i, attitude in enumerate(attitudes):
                if actuator_due:
                    attitude.process_pending_starts(t); attitude.process_completions(t)
                if profiles[i]['control_enabled'] and abs(t-next_control) <= 1e-9:
                    member = state[:, i*n:(i+1)*n]
                    e, v, _ = gap_kinematics(*member, reference_radius_m=config.derived.reference_radius_m,
                        target_spacing_rad=config.derived.target_spacing_rad, mu_m3_s2=config.constants.mu_m3_s2)
                    q = switching_function(e, v, config.derived.controller_kp_s2_inv, config.derived.controller_kd_s_inv)
                    attitude.apply_commands(t, select_commands(q, attitude.mode.copy(), attitude.command_allowed(t),
                        config.controller.threshold_on_m_s2, config.controller.threshold_off_m_s2))
            stable_areas = np.concatenate([a.areas(t) for a in attitudes])
            active_slews = any(np.any((a.mode==1)|(a.mode==3)) for a in attitudes)
            next_actuator = min(a.next_event_time(t) for a in attitudes)
        if abs(t-next_control) <= 1e-9:
            next_control += period
            if next_control > end+1e-9: next_control = float('inf')
        if abs(t-next_output) <= 1e-9:
            arrays['a_m'][k] = state[0].reshape(shape[1:])
            arrays['lambda_rad'][k] = state[1].reshape(shape[1:])
            arrays['density_kg_m3'][k] = density_at_state(t, *state, config,
                disturbance_enabled=disturbance_enabled, density_provider=provider).reshape(shape[1:])
            arrays['area_m2'][k] = areas(t).reshape(shape[1:])
            arrays['mode'][k] = np.asarray([a.mode for a in attitudes])
            minimum_gap = min(minimum_gap, float(np.min(np.diff(arrays['lambda_rad'][k], axis=-1)*config.derived.reference_radius_m)))
            k += 1; next_output += save
        if t >= end-1e-9: break
        stop = next_step_end(t, config.simulation.integration_step_s, end,
            [next_control, next_output, next_actuator] +
            ([] if node_input_metadata is None else node_input_metadata.get('boundaries_s', [])))
        if stop <= t+1e-9: raise RuntimeError('Scheduler did not advance')
        state = rk4_step(rhs, t, state, stop-t); t = stop
        if t >= next_progress:
            for a in arrays.values(): a.flush()
            if progress: progress(t, end, perf_counter()-began,
                [len(a.events) for a in attitudes], minimum_gap)
            next_progress += 86400.
    for a in arrays.values(): a.flush()
    events = clean_json([[e.to_dict() for e in a.events] for a in attitudes])
    write_json(out/'checkpoint.json', dict(configuration=config.to_dict(), profiles=profiles,
        disturbance_enabled=disturbance_enabled, time_s=t, state=state.tolist(), events=events,
        modes=[a.mode.tolist() for a in attitudes], areas=[a.areas(t).tolist() for a in attitudes],
        **({'node_input': node_input_metadata} if node_input_metadata is not None else {})))
    write_json(out/'completion.json', dict(complete=True, start_s=start, end_s=t,
        samples=k, elapsed_s=perf_counter()-began, event_counts=[len(a.events) for a in attitudes],
        minimum_signed_gap_m=minimum_gap))
    return state, events


def derived_histories(t, a, u, area, mode, rho, config):
    """Adapt memory-mapped histories to existing evaluation/residual interfaces."""
    from .long_horizon import decompose_q
    e, v, gap = gap_kinematics(a, u, reference_radius_m=config.derived.reference_radius_m,
        target_spacing_rad=config.derived.target_spacing_rad, mu_m3_s2=config.constants.mu_m3_s2)
    p, d = decompose_q(e, v, config.derived.controller_kp_s2_inv, config.derived.controller_kd_s_inv)
    return dict(time_s=t, a_m=a, lambda_rad=u, area_m2=area, mode=mode,
        density_kg_m3=rho, altitude_m=a-config.constants.earth_radius_m,
        gap_error_m=e, gap_rate_m_s=v, gap_distance_m=gap, q_m_s2=p+d)
