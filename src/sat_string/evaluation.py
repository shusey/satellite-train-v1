"""Finite-observation metrics, computed from histories without changing dynamics."""
from __future__ import annotations

from dataclasses import dataclass, replace
from math import isfinite

import numpy as np

from .attitude import SwitchEvent
from .config import SimulationConfig
from .control_cost import event_activity, propagation_metrics


@dataclass(frozen=True)
class EvaluationSettings:
    """Post-processing tolerances; W is also the minimum settling verification."""

    epsilon_e_m: float
    epsilon_v_m_s: float
    hold_window_s: float

    def __post_init__(self):
        values = (self.epsilon_e_m, self.epsilon_v_m_s, self.hold_window_s)
        if not all(isfinite(v) and v >= 0 for v in values):
            raise ValueError("Evaluation tolerances and hold window must be finite and nonnegative")

    @classmethod
    def from_config(cls, config):
        r = config.recovery
        return cls(r.position_tolerance_m, r.velocity_tolerance_m_s, r.hold_time_s)


def switching_metrics(events, satellite_count: int, end_s: float) -> dict:
    """Count observed events, never scheduled completions beyond the horizon.

    Standard switches = transition starts (legacy convention). A cycle requires
    two COMPLETED alternating transitions from the core's initial LOW state.
    Secondary means >=2 completed L-H-L cycles; multiple means >=3.
    """
    records = [e.to_dict() if isinstance(e, SwitchEvent) else dict(e) for e in events]
    by_sat = [[] for _ in range(satellite_count)]
    for e in records:
        index = int(e["satellite_id"]) - 1
        times = [float(e[k]) for k in ("command_time", "start_time", "complete_time")]
        if not 0 <= index < satellite_count or not all(isfinite(t) for t in times):
            raise ValueError("Invalid switch event satellite ID or time")
        if not 0 <= times[0] <= times[1] < times[2]:
            raise ValueError("Expected 0 <= command <= start < completion")
        if times[0] <= end_s + 1e-9:
            by_sat[index].append(e)
    counts = {k: [] for k in ("command_time", "start_time", "complete_time")}
    intervals = {k: [] for k in counts}
    cycles, last_completed, last_started, high_starts = [], [], [], []
    pending = 0
    for own in by_sat:
        own.sort(key=lambda e: e["command_time"])
        level, last_end, n_cycles = "L", -1.0, 0
        for e in own:
            if e["from_mode"] != level or e["to_mode"] != ("H" if level == "L" else "L"):
                raise ValueError("Events must alternate from the initial LOW state")
            if e["command_time"] < last_end - 1e-9:
                raise ValueError("Commands cannot overlap an unfinished transition")
            level, last_end = e["to_mode"], e["complete_time"]
            if last_end <= end_s + 1e-9 and level == "L":
                n_cycles += 1
            if last_end > end_s + 1e-9:
                pending += 1
        cycles.append(n_cycles)
        high_starts.append(sum(e["to_mode"] == "H" and e["start_time"] <= end_s + 1e-9 for e in own))
        for key in counts:
            observed = [float(e[key]) for e in own if e[key] <= end_s + 1e-9]
            counts[key].append(len(observed))
            intervals[key].extend(np.diff(observed).tolist())
            if key == "complete_time":
                last_completed.extend(observed)
            if key == "start_time":
                last_started.extend(observed)
    maximum = max(cycles, default=0)
    return {
        "standard_switch_definition": "observed_transition_start",
        "n_commands": sum(counts["command_time"]),
        "n_transition_starts": sum(counts["start_time"]),
        "n_transition_completions": sum(counts["complete_time"]),
        "commands_per_satellite": counts["command_time"],
        "switches_per_satellite": counts["start_time"],
        "completions_per_satellite": counts["complete_time"],
        "last_switch_time_s": max(last_completed, default=None),
        "last_transition_start_time_s": max(last_started, default=None),
        "pending_or_active_transitions_at_tmax": pending,
        "min_switch_interval_s": min(intervals["start_time"], default=None),
        "min_command_interval_s": min(intervals["command_time"], default=None),
        "min_completion_interval_s": min(intervals["complete_time"], default=None),
        "lhl_cycles_per_satellite": cycles,
        "max_lhl_cycles": maximum,
        "secondary_switching_detected": maximum >= 2,
        "multiple_switching_detected": maximum >= 3,
        "secondary_satellite_ids": [i + 1 for i, n in enumerate(cycles) if n >= 2],
        "l_to_h_starts_per_satellite": high_starts,
        "secondary_response_started": any(n >= 2 for n in high_starts),
    }


def settling_metrics(time_s, error, rate, switching: dict, settings: EvaluationSettings,
                     *, evaluation_start_s: float) -> dict:
    """Recovery is first confirmed hold; settling is an accepted final suffix.

    Suffix starts after the last completed transition, contains no active/pending
    transition, and must last at least W. Conditions are checked at saved samples,
    using <= for compatibility with legacy recovery. Neither proves permanence.
    """
    t = np.asarray(time_s, dtype=float)
    e, v = np.asarray(error), np.asarray(rate)
    if (t.ndim != 1 or t.size < 2 or not np.all(np.isfinite(t))
            or np.any(np.diff(t) <= 0) or e.shape != v.shape
            or e.ndim != 2 or e.shape[0] != t.size
            or not np.all(np.isfinite(e)) or not np.all(np.isfinite(v))):
        raise ValueError("Expected finite, increasing time and matching (time,gap) arrays")
    within = ((np.max(np.abs(e), axis=1) <= settings.epsilon_e_m)
              & (np.max(np.abs(v), axis=1) <= settings.epsilon_v_m_s)
              & (t >= evaluation_start_s - 1e-9))
    entry = confirmation = None
    run_start = None
    for time, accepted in zip(t, within, strict=True):
        if not accepted:
            run_start = None
        else:
            if run_start is None:
                run_start = float(time)
            if time - run_start >= settings.hold_window_s - 1e-9:
                entry, confirmation = run_start, float(time)
                break
    suffix = np.logical_and.accumulate(within[::-1])[::-1]
    last = switching["last_switch_time_s"]
    eligible = suffix & (t >= (last if last is not None else evaluation_start_s) - 1e-9)
    eligible &= t[-1] - t >= settings.hold_window_s - 1e-9
    if switching["pending_or_active_transitions_at_tmax"]:
        eligible[:] = False
    indices = np.flatnonzero(eligible)
    settled = float(t[indices[0]]) if indices.size else None
    return {
        "recovery_time_s": confirmation,
        "recovery_entry_time_s": entry,
        "recovery_time_after_disturbance_s": (confirmation - evaluation_start_s
                                               if confirmation is not None else None),
        "settling_time_s": settled,
        "settled_by_tmax": settled is not None,
        "verification_horizon_s": float(t[-1] - settled) if settled is not None else None,
        "epsilon_e_m": settings.epsilon_e_m,
        "epsilon_v_m_s": settings.epsilon_v_m_s,
        "hold_window_s": settings.hold_window_s,
        "evaluation_start_time_s": evaluation_start_s,
        "observation_end_time_s": float(t[-1]),
        "maximum_sample_interval_s": float(np.max(np.diff(t))),
        "interpretation": "finite observed suffix; no guarantee beyond tmax or between samples",
    }


def compute_evaluation_metrics(arrays, config: SimulationConfig, events=(), *,
                               settings: EvaluationSettings | None = None,
                               low_baseline=None, evaluation_start_s=None) -> dict:
    """Additive schema; supports saved prefixes without rerunning a simulation."""
    t, e = np.asarray(arrays["time_s"]), np.asarray(arrays["gap_error_m"])
    a = np.asarray(arrays["a_m"])
    if not np.isclose(t[0], 0.0, atol=1e-9, rtol=0):
        raise ValueError("Evaluation requires a zero-origin history")
    end = float(t[-1])
    settings = settings or EvaluationSettings.from_config(config)
    switching = switching_metrics(events, config.formation.satellite_count, end)
    activity = event_activity(events, config.formation.satellite_count, end)
    switching.update({k: v for k, v in activity.items() if k.startswith("high_")})
    # Sample RMS: every saved time and every adjacent gap has equal weight.
    formation = {
        "max_abs_error_m": float(np.max(np.abs(e))),
        "rms_error_m": float(np.sqrt(np.mean(e**2))),
        "rms_error_per_gap_m": np.sqrt(np.mean(e**2, axis=0)).tolist(),
        "final_max_abs_error_m": float(np.max(np.abs(e[-1]))),
        "final_max_abs_rate_m_s": float(np.max(np.abs(arrays["gap_rate_m_s"][-1]))),
        "rms_definition": "sqrt(mean over all saved samples and adjacent gaps of error squared)",
        "time_weighted_rms_error_m": float(np.sqrt(np.mean(np.trapezoid(e**2, t, axis=0) / end))),
    }
    horizon_config = replace(config, simulation=replace(config.simulation, duration_s=end))
    wave = propagation_metrics(events, horizon_config)
    first = wave["first_switch_start_s_by_satellite"]
    arrivals = [(i + 1, v) for i, v in enumerate(first) if v is not None]
    last_arrival = max((v for _, v in arrivals), default=None)
    center_id = wave["center_satellite_id"]
    for side in ("left", "right"):
        own = [(i, v) for i, v in arrivals if (i < center_id if side == "left" else i > center_id)]
        final = max((v for _, v in own), default=None)
        wave[side].update(
            max_hops=max((abs(i - center_id) for i, _ in own), default=0),
            last_arrival_time_s=final,
            last_arrival_satellite_ids=[i for i, v in own if v == final],
        )
    wave.update({
        "left_speed_m_s": wave["left"]["speed_m_s"],
        "right_speed_m_s": wave["right"]["speed_m_s"],
        "max_hops": wave["range_hops"],
        "last_arrival_time_s": last_arrival,
        "last_arrival_satellite_ids": [i for i, v in arrivals if v == last_arrival],
    })
    losses = a[0] - a[-1]
    altitude = {
        "mean_altitude_loss_m": float(np.mean(losses)),
        "altitude_loss_per_satellite_m": losses.tolist(),
        "evaluation_time_s": end,
        "loss_definition": "initial minus final altitude; positive is loss; arithmetic satellite mean",
        "extra_loss_vs_all_low_m": None,
        "extra_loss_per_satellite_vs_all_low_m": None,
        "extra_loss_definition": "mean(a_all_low(T)-a_case(T)); matched initial conditions, pulse and endpoint",
    }
    if low_baseline is not None:
        if not np.array_equal(t, low_baseline["time_s"]) or a.shape != low_baseline["a_m"].shape:
            raise ValueError("All-LOW reference must have matching time grid and satellite count")
        if not np.allclose(a[0], low_baseline["a_m"][0], atol=1e-9, rtol=0):
            raise ValueError("All-LOW reference must have matching initial radii")
        extra = low_baseline["a_m"][-1] - a[-1]
        altitude.update(extra_loss_vs_all_low_m=float(np.mean(extra)),
                        extra_loss_per_satellite_vs_all_low_m=extra.tolist())
    start = (config.disturbance.start_time_s + config.disturbance.duration_s
             if evaluation_start_s is None else float(evaluation_start_s))
    return {
        "evaluation_schema_version": 1,
        "formation": formation, "switching": switching,
        "settling": settling_metrics(t, e, arrays["gap_rate_m_s"], switching,
                                      settings, evaluation_start_s=start),
        "wave": wave, "altitude": altitude,
    }
