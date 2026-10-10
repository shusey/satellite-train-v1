"""All-LOW comparison with passive RK4 exposure quadratures.

Physical derivatives and RK4 are reused unchanged. Extra rows integrate
diagnostics and never feed back into a or lambda.
"""
from __future__ import annotations

import numpy as np

from .atmosphere import background_density
from .controller import initial_formation
from .dynamics import density_at_state, state_derivative
from .integrator import next_step_end, rk4_step
from .simulator import CaseDefinition, SimulationResult, _formation_quantities


def local_background(radius, config):
    return background_density(
        radius - config.constants.earth_radius_m,
        config.formation.reference_altitude_m,
        config.atmosphere.reference_density_kg_m3,
        config.atmosphere.effective_scale_height_m,
    )


def exposure_rhs(time_s, state, config):
    """Rows: a, lambda, J_rho, J_delta, J_drag, J_drag_delta.

    J_drag integrates signed tangential acceleration -rho*Cd*A*mu/(2*m*a),
    in m/s. Background is evaluated at the SAME actual satellite state,
    not a separate undisturbed trajectory. Density quadratures are kg*s/m^3.
    """
    area = np.full(config.formation.satellite_count, config.spacecraft.area_low_m2)
    physical = state_derivative(time_s, state[:2], config, lambda _: area,
                                disturbance_enabled=True)
    coefficient = config.spacecraft.drag_coefficient * area / config.spacecraft.mass_kg
    # Recover the same stage density from da/dt, avoiding a second model call.
    rho = -physical[0] / (coefficient * np.sqrt(config.constants.mu_m3_s2 * state[0]))
    delta = rho - local_background(state[0], config)
    if not config.disturbance.enabled:
        # Disabled forcing is exactly zero; do not integrate subtraction roundoff.
        delta = np.zeros_like(rho)
    drag_factor = -0.5 * coefficient * config.constants.mu_m3_s2 / state[0]
    return np.vstack((physical, rho, delta, drag_factor * rho, drag_factor * delta))


def simulate_all_low_exposure(config, progress=None):
    """Use the original initial state, max step and output-event boundaries."""
    if config.atmosphere.density_model != "exponential":
        raise ValueError("Passive exposure runner requires the exponential atmosphere")
    n = config.formation.satellite_count
    a, longitude = initial_formation(n, config.derived.reference_radius_m,
                                    config.derived.target_spacing_rad)
    state = np.vstack((a, longitude, np.zeros((4, n))))
    times, states = [], []
    t, next_output, last_progress = 0.0, 0.0, -1
    end = config.simulation.duration_s
    while True:
        if abs(t - next_output) < 1e-9:
            times.append(t)
            states.append(state.copy())
            next_output = min(next_output + config.simulation.output_step_s, end)
        if t >= end - 1e-9:
            break
        stop = next_step_end(t, config.simulation.integration_step_s, end, (next_output,))
        state = rk4_step(lambda ts, y: exposure_rhs(ts, y, config), t, state, stop - t)
        t = stop
        bucket = int(t / end * 10)
        if progress and bucket > last_progress:
            progress(t, end)
            last_progress = bucket
    times, states = np.asarray(times), np.asarray(states)
    error, rate, distance, q = zip(*[_formation_quantities(s[:2], config) for s in states])
    result = SimulationResult(
        CaseDefinition("all_low_" + config.disturbance.type, False, True),
        times, states[:, 0], states[:, 1],
        states[:, 0] - config.constants.earth_radius_m,
        np.full((len(times), n), config.spacecraft.area_low_m2),
        np.zeros((len(times), n), dtype=np.int8),
        np.asarray([density_at_state(ts, s[0], s[1], config, disturbance_enabled=True)
                    for ts, s in zip(times, states)]),
        np.asarray(q), np.asarray(error), np.asarray(rate), np.asarray(distance), [],
    )
    exposure = dict(time_s=times, J_rho=states[:, 2], J_delta=states[:, 3],
                    J_drag=states[:, 4], J_drag_delta=states[:, 5])
    return result, exposure


def snapshot_metrics(error, rate):
    return dict(max_abs_error_m=float(np.max(np.abs(error))),
                rms_error_m=float(np.sqrt(np.mean(error**2))),
                max_abs_rate_m_s=float(np.max(np.abs(rate))),
                rms_rate_m_s=float(np.sqrt(np.mean(rate**2))),
                gap_error_m=error.tolist(), gap_rate_m_s=rate.tolist())


def evaluate_motion(arrays, exposure, config):
    """Gaussian passage end: all satellites downstream by at least 6 sigma.

    Tail/A <= exp(-18), not literal finite support. Exact comoving has
    no passage end within this observation; ongoing forcing is kept separate.
    Forecasts freeze the sampled gap rate and are NOT new physical runs.
    """
    t = arrays["time_s"]
    e, v = arrays["gap_error_m"], arrays["gap_rate_m_s"]
    x = config.derived.reference_radius_m * arrays["lambda_rad"]
    d = config.disturbance
    offset = x - (d.center_s_m + d.propagation_speed_m_s * t[:, None])
    passed = np.all(offset >= 6 * config.derived.disturbance_sigma_m, axis=1)
    indices = np.flatnonzero(passed & np.logical_and.accumulate(passed[::-1])[::-1])
    finish = int(indices[0]) if indices.size else None
    relative_speed = (config.derived.reference_radius_m
                      * config.derived.reference_mean_motion_rad_s - d.propagation_speed_m_s)
    snapshots = {}
    if finish is not None:
        for name, delay in (("passage_end", 0), ("three_hours_after", 10800),
                            ("one_day_after", 86400)):
            target = t[finish] + delay
            idx = int(np.searchsorted(t, target))
            snapshots[name] = (dict(time_s=float(t[idx]), **snapshot_metrics(e[idx], v[idx]))
                               if idx < len(t) else None)
    growth = None
    if finish is not None:
        fit_mask = t >= t[finish] + 10800
        if np.count_nonzero(fit_mask) >= 2:
            slopes = np.polyfit(t[fit_mask]-t[fit_mask][0], e[fit_mask], 1)[0]
            growth = dict(window_start_s=float(t[fit_mask][0]), window_end_s=float(t[-1]),
                          fitted_gap_growth_rate_m_s=slopes.tolist(),
                          max_abs_fitted_growth_rate_m_s=float(np.max(np.abs(slopes))),
                          caution="Very small slopes can be floating-point accumulation, not resolved physical drift")
    end = snapshot_metrics(e[-1], v[-1])
    forecast = e[-1] + v[-1] * (30 * 86400 - t[-1])
    jr, jd = exposure["J_rho"][-1], exposure["J_delta"][-1]
    mean_delta = float(np.mean(jd))
    drag, delta_drag = exposure["J_drag"][-1], exposure["J_drag_delta"][-1]
    predicted_radius = (np.sqrt(config.derived.reference_radius_m)
                        - 0.5 * config.spacecraft.drag_coefficient
                        * config.spacecraft.area_low_m2 / config.spacecraft.mass_kg
                        * np.sqrt(config.constants.mu_m3_s2) * jr)**2
    predicted_rate = config.derived.reference_radius_m * np.diff(
        np.sqrt(config.constants.mu_m3_s2 / predicted_radius**3))
    increment = arrays["density_kg_m3"] - local_background(arrays["a_m"], config)
    nominal_impulse = (config.atmosphere.reference_density_kg_m3
                       * d.amplitude_fraction * np.sqrt(2*np.pi)
                       * config.derived.disturbance_sigma_m / abs(relative_speed)
                       if relative_speed else None)
    return dict(
        formation=dict(max_abs_error_m=float(np.max(np.abs(e))),
                       max_error_percent_of_spacing=float(np.max(np.abs(e)) / config.formation.target_spacing_m * 100),
                       rms_error_m=float(np.sqrt(np.mean(e**2))), final=end,
                       max_error_per_gap_m=np.max(np.abs(e), axis=0).tolist()),
        motion=dict(disturbance_speed_m_s=d.propagation_speed_m_s,
                    initial_satellite_speed_m_s=relative_speed+d.propagation_speed_m_s,
                    relative_speed_m_s=relative_speed,
                    fwhm_encounter_time_s=d.fwhm_m/abs(relative_speed) if relative_speed else None,
                    train_transit_delay_s=(config.formation.satellite_count-1)*config.formation.target_spacing_m/abs(relative_speed) if relative_speed else None,
                    passage_end_s=float(t[finish]) if finish is not None else None,
                    forcing_ongoing_at_end=finish is None),
        exposure=dict(J_rho_kg_s_m3=jr.tolist(), J_delta_kg_s_m3=jd.tolist(),
                      J_delta_spread_kg_s_m3=float(np.ptp(jd)),
                      mean_J_delta_kg_s_m3=mean_delta,
                      relative_J_delta_spread=float(np.ptp(jd)/mean_delta) if mean_delta else None,
                      J_rho_spread_kg_s_m3=float(np.ptp(jr)),
                      adjacent_J_rho_kg_s_m3=np.diff(jr).tolist(),
                      frozen_background_full_passage_J_delta_kg_s_m3=nominal_impulse,
                      max_instant_increment_difference_kg_m3=float(np.max(np.ptp(increment,axis=1))),
                      max_instant_adjacent_increment_difference_kg_m3=float(np.max(np.abs(np.diff(increment,axis=1)))),
                      max_instant_gap_rate_m_s=float(np.max(np.abs(v))),
                      adjacent_J_delta_kg_s_m3=np.diff(jd).tolist(),
                      J_drag_m_s=drag.tolist(), J_drag_delta_m_s=delta_drag.tolist(),
                      J_drag_spread_m_s=float(np.ptp(drag)),
                      J_drag_delta_spread_m_s=float(np.ptp(delta_drag)),
                      adjacent_J_drag_delta_m_s=np.diff(delta_drag).tolist(),
                      total_exposure_radius_identity_max_error_m=float(np.max(np.abs(predicted_radius-arrays["a_m"][-1]))),
                      total_exposure_predicted_gap_rate_m_s=predicted_rate.tolist(),
                      total_exposure_rate_identity_max_error_m_s=float(np.max(np.abs(predicted_rate-v[-1]))),
                      minus_three_adjacent_incremental_drag_m_s=(-3*np.diff(delta_drag)).tolist()),
        post_passage=snapshots,
        post_passage_growth=growth,
        frozen_rate_forecast_30d=dict(
            applicable_only_after_passage=finish is not None,
            max_abs_error_m=float(np.max(np.abs(forecast))) if finish is not None else None,
            gap_error_m=forecast.tolist() if finish is not None else None),
    )
