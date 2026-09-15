"""Offline formation/control/orbital-cost analysis; no changes to simulation core."""
from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any
import numpy as np

from .attitude import SwitchEvent
from .config import SimulationConfig
from .switching_wave import analyze_switching_wave, first_l_to_h_start_times


def clean_json(value):
    """Represent unavailable estimates as JSON null rather than nonstandard NaN."""
    if isinstance(value, Mapping):
        return {key: clean_json(item) for key, item in value.items()}
    if isinstance(value, (list, tuple, np.ndarray)):
        return [clean_json(item) for item in value]
    if isinstance(value, (float, np.floating)):
        return float(value) if np.isfinite(value) else None
    if isinstance(value, np.integer):
        return int(value)
    return value


def event_activity(events: Sequence, count: int, end_s: float) -> dict:
    """Exact integrals on [0,end] for the core's initially-LOW smoothstep actuator.

    Switch count = actual transition starts, including a start exactly at end.
    Stable HIGH excludes both slews. Equivalent HIGH = integral((A-AL)/(AH-AL)).
    Pending starts beyond end do not count. Partial final slews are integrated.
    """
    records = [e.to_dict() if isinstance(e, SwitchEvent) else dict(e) for e in events]
    high = np.zeros(count)
    equivalent = np.zeros(count)
    switches = np.zeros(count, dtype=int)
    for index in range(count):
        own = sorted((e for e in records if int(e["satellite_id"]) == index + 1),
                     key=lambda e: float(e["start_time"]))
        last = 0.0
        level = 0.0
        for event in own:
            start = float(event["start_time"])
            complete = float(event["complete_time"])
            if start > end_s:
                break
            if start < last - 1e-9 or complete <= start:
                raise ValueError("Overlapping or nonpositive-duration transitions")
            expected = "H" if level else "L"
            if event["from_mode"] != expected or event["to_mode"] == expected:
                raise ValueError("Events must form an initially-LOW alternating sequence")
            switches[index] += 1
            high[index] += level * (start - last)
            equivalent[index] += level * (start - last)
            duration = complete - start
            z = min((end_s - start) / duration, 1.0)
            # Integral_0^z (3u^2 - 2u^3) du = z^3 - z^4/2.
            integral = duration * (z**3 - 0.5 * z**4)
            equivalent[index] += integral if level == 0 else duration * z - integral
            last = min(complete, end_s)
            if complete > end_s:
                break
            level = 1.0 - level
        if last < end_s:
            high[index] += level * (end_s - last)
            equivalent[index] += level * (end_s - last)
    return {
        "switch_count_by_satellite": switches.tolist(),
        "switch_count": int(switches.sum()),
        "high_stable_s_by_satellite": high.tolist(),
        "high_stable_total_s": float(high.sum()),
        "high_equivalent_s_by_satellite": equivalent.tolist(),
        "high_equivalent_total_s": float(equivalent.sum()),
    }


def propagation_metrics(events, config, center_satellite_id=None):
    """Report observed reach, finite-horizon censoring, and vetted first-front fits.

    No arbitrary density cutoff is used to claim a causally unaffected satellite.
    'Distal' below means |i-center|>=4, matching the existing fit exclusion.
    """
    wave = analyze_switching_wave(events, config, center_satellite_id=center_satellite_id)
    count = config.formation.satellite_count
    center = wave["center_satellite_id"] - 1
    first = first_l_to_h_start_times(events, count,
                                    simulation_end_time_s=config.simulation.duration_s)
    reached = np.flatnonzero(np.isfinite(first))
    offsets = np.abs(reached - center)
    result = {
        "center_satellite_id": center + 1,
        "first_switch_start_s_by_satellite": first.tolist(),
        "reached_satellite_count": int(reached.size),
        "range_hops": int(offsets.max()) if offsets.size else 0,
        "range_m": float(offsets.max() * config.formation.target_spacing_m) if offsets.size else 0.0,
        "distal_satellite_ids": (np.flatnonzero(np.abs(np.arange(count) - center) >= 4) + 1).tolist(),
        "time_censoring": "No arrival by simulation end is not proof of permanent propagation failure.",
    }
    for side in ("left", "right"):
        indices = (np.arange(center - 1, -1, -1) if side == "left"
                   else np.arange(center + 1, count))
        contiguous = []
        first_missing = None
        for i in indices:
            if not np.isfinite(first[i]):
                first_missing = int(i + 1)
                break
            contiguous.append(int(i + 1))
        fit = wave[side]["fit"]
        times = np.asarray(wave[side]["transition_start_times_s"])
        fit_ids = wave[side]["fit_satellite_ids"]
        contiguous_fit = all(i in contiguous for i in fit_ids)
        reasons = []
        if fit["point_count"] < 3:
            reasons.append("fewer_than_3_distal_points")
        if times.size > 1 and not np.all(np.diff(times) > 0):
            reasons.append("nonmonotone_outward_first_arrivals")
        if not contiguous_fit:
            reasons.append("gap_in_front")
        if not np.isfinite(fit["r_squared"]) or fit["r_squared"] < 0.8:
            reasons.append("R2_below_0.8_or_undefined")
        if not np.isfinite(fit["speed_m_s"]) or fit["speed_m_s"] <= 0:
            reasons.append("nonpositive_or_undefined_speed")
        result[side] = {
            "switched_satellite_ids": wave[side]["switched_satellite_ids"],
            "contiguous_last_satellite_id": contiguous[-1] if contiguous else None,
            "first_unreached_satellite_id": first_missing,
            "boundary_reached": bool(indices.size and indices[-1] + 1 in contiguous),
            "observed_front_status": ("boundary_reached" if first_missing is None and indices.size
                                      else "not_reached_by_end"),
            "speed_m_s": fit["speed_m_s"] if not reasons else None,
            "fit_status": "estimated" if not reasons else ";".join(reasons),
            "raw_fit": fit,
            "fit_satellite_ids": fit_ids,
            "adjacent_delays": wave[side]["delta_t_records"],
        }
    return clean_json(result)


def evaluate_cost(arrays, config: SimulationConfig, events=(), low_baseline=None):
    """Compute signed orbital changes and endpoint losses (positive means loss)."""
    t = np.asarray(arrays["time_s"])
    if t[0] != 0 or not np.isclose(t[-1], config.simulation.duration_s):
        raise ValueError("Cost analysis requires a complete zero-origin history")
    a = np.asarray(arrays["a_m"])
    e = np.asarray(arrays["gap_error_m"])
    duration = float(t[-1])
    activity = event_activity(events, a.shape[1], duration)
    wave = propagation_metrics(events, config)
    delta_a = a[-1] - a[0]
    # Algebraic difference avoids cancellation between ~-3e7 J/kg values.
    delta_eps = config.constants.mu_m3_s2 / 2 * (a[-1] - a[0]) / (a[0] * a[-1])
    equivalent_sampled = np.trapezoid(
        (arrays["area_m2"] - config.spacecraft.area_low_m2)
        / (config.spacecraft.area_high_m2 - config.spacecraft.area_low_m2), t, axis=0)
    rms_by_gap = np.sqrt(np.trapezoid(e**2, t, axis=0) / duration)
    out = {
        "maximum_gap_error_m": float(np.max(np.abs(e))),
        "maximum_gap_error_m_by_gap": np.max(np.abs(e), axis=0).tolist(),
        "rms_gap_error_m": float(np.sqrt(np.mean(rms_by_gap**2))),
        "rms_gap_error_m_by_gap": rms_by_gap.tolist(),
        "final_maximum_gap_error_m": float(np.max(np.abs(e[-1]))),
        "final_maximum_gap_rate_m_s": float(np.max(np.abs(arrays["gap_rate_m_s"][-1]))),
        "minimum_gap_m": float(np.min(arrays["gap_distance_m"])),
        "activity": activity,
        "propagation": wave,
        "delta_a_m_by_satellite": delta_a.tolist(),
        "delta_altitude_m_by_satellite": delta_a.tolist(),
        "mean_a_loss_m": float(-delta_a.mean()),
        "delta_specific_energy_J_kg_by_satellite": delta_eps.tolist(),
        "delta_energy_J_by_satellite": (config.spacecraft.mass_kg * delta_eps).tolist(),
        "total_energy_loss_J": float(-config.spacecraft.mass_kg * delta_eps.sum()),
        "sampled_vs_exact_equivalent_max_abs_s": float(np.max(np.abs(
            equivalent_sampled - np.asarray(activity["high_equivalent_s_by_satellite"])))),
    }
    distal = np.asarray(wave["distal_satellite_ids"], dtype=int) - 1
    out["distal_high_equivalent_s"] = float(np.asarray(activity["high_equivalent_s_by_satellite"])[distal].sum())
    out["distal_switch_count"] = int(np.asarray(activity["switch_count_by_satellite"])[distal].sum())
    if low_baseline is not None:
        if not np.array_equal(t, low_baseline["time_s"]):
            raise ValueError("Comparison time grids differ")
        ref_a = low_baseline["a_m"]
        extra = ref_a[-1] - a[-1]
        extra_eps = config.constants.mu_m3_s2 / 2 * extra / (a[-1] * ref_a[-1])
        # da/dt=-k(state)*A. Split loss into extra area at controlled state,
        # and residual density/radius changes (also forcing differences for Case C).
        k = (config.spacecraft.drag_coefficient / config.spacecraft.mass_kg
             * arrays["density_kg_m3"] * np.sqrt(config.constants.mu_m3_s2 * a))
        ref_k = (config.spacecraft.drag_coefficient / config.spacecraft.mass_kg
                 * low_baseline["density_kg_m3"] * np.sqrt(config.constants.mu_m3_s2 * ref_a))
        direct = np.trapezoid(k * (arrays["area_m2"] - config.spacecraft.area_low_m2), t, axis=0)
        feedback = np.trapezoid((k - ref_k) * config.spacecraft.area_low_m2, t, axis=0)
        out.update({
            "extra_a_loss_vs_low_m_by_satellite": extra.tolist(),
            "extra_mean_a_loss_vs_low_m": float(extra.mean()),
            "extra_max_a_loss_vs_low_m": float(extra.max()),
            "extra_energy_loss_vs_low_J": float(config.spacecraft.mass_kg * extra_eps.sum()),
            "extra_energy_loss_vs_low_J_by_satellite": (config.spacecraft.mass_kg * extra_eps).tolist(),
            "direct_extra_area_loss_m_by_satellite": direct.tolist(),
            "density_radius_residual_loss_m_by_satellite": feedback.tolist(),
            "cost_balance_max_residual_m": float(np.max(np.abs(extra - direct - feedback))),
            "distal_extra_a_loss_sum_m": float(extra[distal].sum()),
        })
    return clean_json(out)


def linear_modes(config: SimulationConfig):
    """Frozen tangent model about unforced common all-LOW decay at a=a0.

    Ideal signed incremental area delta A=q/g (not implemented or feasible at
    the LOW boundary). No hybrid stability claim follows from these eigenvalues.
    """
    a = config.derived.reference_radius_m
    rho = config.atmosphere.reference_density_kg_m3
    mu = config.constants.mu_m3_s2
    area = config.spacecraft.area_low_m2
    f = -rho * config.spacecraft.drag_coefficient * area / config.spacecraft.mass_kg * np.sqrt(mu*a)
    b = f * (-1/config.atmosphere.effective_scale_height_m - 2/a)
    g = 1.5 * config.derived.reference_radius_m * mu / a**2 * rho * config.spacecraft.drag_coefficient / config.spacecraft.mass_kg
    kp = config.derived.controller_kp_s2_inv
    kd = config.derived.controller_kd_s_inv
    count = config.formation.satellite_count
    rows = []
    for k in range(count):
        eigenvalue = 2 * (1 - np.cos(k*np.pi/count))
        roots = np.roots([1, kd*eigenvalue-b, kp*eigenvalue])
        rows.append({
            "mode": k, "laplacian_eigenvalue": float(eigenvalue),
            "real_1_s_inv": float(roots[0].real), "imag_1_s_inv": float(roots[0].imag),
            "real_2_s_inv": float(roots[1].real), "imag_2_s_inv": float(roots[1].imag),
        })
    return {
        "reference": "frozen tangent about unforced common LOW decay, t=0",
        "background_decay_m_s": float(f),
        "b_s_inv": float(b), "g_m_inv_s2": float(g),
        "high_minus_low_virtual_acceleration_m_s2": float(g * (config.spacecraft.area_high_m2-area)),
        "stability_condition": "Kp>0 and Kd*lambda_2>b, disagreement modes only",
        "disagreement_stable": bool(kp>0 and kd*rows[1]["laplacian_eigenvalue"]>b),
        "modes": rows,
    }
