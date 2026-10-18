"""Offline episode and q analysis. Episodes are temporal groups, not causal proofs."""
from __future__ import annotations

from dataclasses import replace
import numpy as np

from .attitude import SwitchEvent
from .control_cost import clean_json
from .switching_wave import analyze_switching_wave, infer_center_satellite_index


def decompose_q(error, rate, kp, kd):
    """Same open-path incidence operation as controller.switching_function."""
    def incidence(edges):
        return np.column_stack((edges[:, 0], np.diff(edges, axis=1), -edges[:, -1]))
    e, v = np.asarray(error, dtype=float), np.asarray(rate, dtype=float)
    if e.ndim != 2 or e.shape != v.shape or not e.shape[1]:
        raise ValueError("Expected matching (time,gap) arrays")
    return incidence(kp*e), incidence(kd*v)


def _records(events):
    return sorted((e.to_dict() if isinstance(e, SwitchEvent) else dict(e) for e in events),
                  key=lambda e: (e["start_time"], e["satellite_id"]))


def infer_episode_gap(events):
    """Learn the natural gap separation from the saved 14-day reference.

    Geometric midpoint at largest ratio of sorted positive L->H arrival gaps.
    A ratio >=10 is an analysis quality gate, not a new physical parameter.
    Freeze the inferred gap for all longer horizons; do not retune each run.
    """
    starts = sorted(e["start_time"] for e in _records(events) if e["to_mode"] == "H")
    gaps = np.sort(np.diff(starts))
    gaps = gaps[gaps > 0]
    if gaps.size < 2:
        raise ValueError("Insufficient arrivals to infer episode separation")
    ratios = gaps[1:]/gaps[:-1]
    index = int(np.argmax(ratios))
    if ratios[index] < 10:
        raise ValueError("No clear gap separation; specify an analysis gap explicitly")
    return dict(gap_s=float(np.sqrt(gaps[index]*gaps[index+1])),
                largest_intra_gap_s=float(gaps[index]), smallest_inter_gap_s=float(gaps[index+1]),
                separation_ratio=float(ratios[index]))


def group_episodes(events, end_s, gap_s):
    """Group rising starts; attach each return to its own satellite's rise.

    Unlike per-satellite cycle number this correctly includes a late first
    response of the central satellite in the second GLOBAL episode.
    """
    if not np.isfinite(gap_s) or gap_s <= 0:
        raise ValueError("Episode gap must be positive")
    groups, active = [], {}
    last_rise = None
    for e in _records(events):
        if e["command_time"] > end_s + 1e-9:
            continue
        i = int(e["satellite_id"])
        if e["to_mode"] == "H":
            if e["start_time"] > end_s + 1e-9:
                continue
            if i in active:
                raise ValueError("Missing return before the next rise")
            if last_rise is None or e["start_time"]-last_rise > gap_s:
                groups.append([])
            active[i] = len(groups)-1
            groups[-1].append(e)
            last_rise = e["start_time"]
        elif e["to_mode"] == "L":
            if i not in active:
                raise ValueError("Return without an observed rise")
            groups[active.pop(i)].append(e)
        else:
            raise ValueError("Unknown event target")
    return groups


def _fit_outer_branches(records, config, end_s):
    c = replace(config, simulation=replace(config.simulation, duration_s=end_s))
    wave = analyze_switching_wave(records, c)
    result = {}
    first = np.asarray(wave["first_transition_start_times_s"])
    for side in ("left", "right"):
        s = wave[side]
        fit = s["fit"]
        ids = s["fit_satellite_ids"]
        arrivals = first[np.asarray(ids, dtype=int)-1]
        valid = (len(ids) >= 3 and np.all(np.diff(arrivals)>0)
                 and np.all(np.abs(np.diff(ids)) == 1)
                 and np.isfinite(fit["r_squared"]) and fit["r_squared"] >= .8
                 and np.isfinite(fit["speed_m_s"]) and fit["speed_m_s"]>0)
        result[side] = dict(speed_m_s=fit["speed_m_s"] if valid else None,
                           fit=fit, fit_satellite_ids=ids,
                           interpretation="distal outward branch, excludes three center hops; not one front covering inner delayed arrivals")
    return clean_json(result)


def _sample_index(t, time):
    index = int(np.searchsorted(t, time))
    return index if index < len(t) and abs(t[index]-time)<=1e-8 else None


def analyze_long_horizon(arrays, config, events, low_arrays, gap_s):
    t = np.asarray(arrays["time_s"])
    e, v = arrays["gap_error_m"], arrays["gap_rate_m_s"]
    if not np.array_equal(t, low_arrays["time_s"]):
        raise ValueError("Matched all-LOW time grid is required")
    qp, qd = decompose_q(e, v, config.derived.controller_kp_s2_inv, config.derived.controller_kd_s_inv)
    residual = float(np.max(np.abs(qp+qd-arrays["q_m_s2"])))
    if residual > 1e-20:
        raise ValueError("P/D decomposition disagrees with stored q")
    extra = np.mean(low_arrays["a_m"]-arrays["a_m"], axis=1)
    groups = group_episodes(events, float(t[-1]), gap_s)
    summaries, triggers, endpoints = [], [], []
    count, center = config.formation.satellite_count, infer_center_satellite_index(config)+1
    all_events = _records(events)
    previous_end = None
    for number, records in enumerate(groups, 1):
        rising = [r for r in records if r["to_mode"]=="H" and r["start_time"]<=t[-1]+1e-9]
        returning = [r for r in records if r["to_mode"]=="L" and r["start_time"]<=t[-1]+1e-9]
        complete = len(returning)==len(rising) and all(r["complete_time"]<=t[-1]+1e-9 for r in records)
        start = min(r["start_time"] for r in rising)
        scheduled_end = max(r["complete_time"] for r in records)
        end = scheduled_end if complete else float(t[-1])
        ks = _sample_index(t, start); ke = _sample_index(t,end)
        if ks is None or ke is None:
            raise ValueError("Episode boundary is not on saved output grid")
        active = slice(ks,ke+1)
        before = max(ks-1,0)
        first = np.full(count,np.nan)
        for r in rising:
            i=r["satellite_id"]-1
            first[i]=min(first[i],r["start_time"]) if np.isfinite(first[i]) else r["start_time"]
        origins=[]
        for i, time in enumerate(first):
            neighbors=[first[j] if np.isfinite(first[j]) else np.inf for j in (i-1,i+1) if 0<=j<count]
            if np.isfinite(time) and all(time<=n for n in neighbors) and any(time<n for n in neighbors):
                origins.append(i+1)
        seed_ids=[r["satellite_id"] for r in rising if r["start_time"]==start]
        per_sat=[sum(r["satellite_id"]==i+1 and r["start_time"]<=t[-1]+1e-9 for r in records) for i in range(count)]
        summary=dict(wave_number=number,start_time_s=float(start),
                     first_satellite_ids=seed_ids, origin_candidates=origins,
                     complete=complete, end_time_s=float(end) if complete else None,
                     observed_duration_s=float(end-start), transition_starts=sum(per_sat),
                     commands=len(records), transition_completions=sum(r["complete_time"]<=t[-1]+1e-9 for r in records),
                     l_to_h_satellite_ids=sorted({r["satellite_id"] for r in rising}),
                     h_to_l_satellite_ids=sorted({r["satellite_id"] for r in returning}),
                     l_to_h_satellite_count=len({r["satellite_id"] for r in rising}),
                     h_to_l_satellite_count=len({r["satellite_id"] for r in returning}),
                     switches_per_satellite=per_sat,first_l_to_h_times_s=first.tolist(),
                     max_hops_from_disturbance_center=max(abs(r["satellite_id"]-center) for r in rising),
                     max_graph_distance_from_first_trigger=max(abs(r["satellite_id"]-seed_ids[0]) for r in rising),
                     error_before_start_m=float(np.max(np.abs(e[before]))),
                     E_n_active_max_error_m=float(np.max(np.abs(e[active]))),
                     active_rms_error_m=float(np.sqrt(np.mean(e[active]**2))),
                     error_at_observed_end_m=float(np.max(np.abs(e[ke]))),
                     extra_loss_before_start_m=float(extra[before]),
                     extra_loss_at_end_m=float(extra[ke]),
                     extra_loss_active_increment_m=float(extra[ke]-extra[before]),
                     extra_loss_since_previous_wave_end_m=float(extra[ke]-extra[_sample_index(t,previous_end)]) if previous_end is not None else float(extra[ke]),
                     outer_branch_fits=_fit_outer_branches(records,config,float(t[-1])))
        rise_counts=[sum(r["satellite_id"]==i+1 for r in rising) for i in range(count)]
        summary["l_to_h_starts_per_satellite"]=rise_counts
        summary["contains_repeated_local_cycles"]=max(rise_counts,default=0)>1
        if number<len(groups):
            next_start=min(r["start_time"] for r in groups[number] if r["to_mode"]=="H")
            next_k=_sample_index(t,next_start)
            summary["next_start_interval_s"]=float(next_start-start)
            summary["quiet_interval_max_error_m"]=float(np.max(np.abs(e[ke:next_k]))) if next_k>ke else None
            summary["associated_start_to_next_max_error_m"]=float(np.max(np.abs(e[ks:next_k])))
            summary["overlaps_next_episode"]=next_start<end
        else:
            summary["next_start_interval_s"]=None
            summary["quiet_interval_max_error_m"]=float(np.max(np.abs(e[ke:]))) if complete else None
            summary["associated_start_to_next_max_error_m"]=None
            summary["overlaps_next_episode"]=False
        summaries.append(summary)
        if complete:
            endpoints.append(dict(wave_number=number,time_s=float(end),error_m=e[ke].tolist(),
                                   rate_m_s=v[ke].tolist(),relative_a_m=(arrays["a_m"][ke]-np.mean(arrays["a_m"][ke])).tolist()))
        for r in rising:
            k=_sample_index(t,r["command_time"]);j=r["satellite_id"]-1
            if k is None:
                continue
            previous=[z["complete_time"] for z in all_events if z["satellite_id"]==j+1 and z["complete_time"]<r["command_time"]-1e-9]
            own_return=next((z for z in records if z["satellite_id"]==j+1 and z["to_mode"]=="L" and z["command_time"]>r["command_time"]),None)
            record=dict(wave_number=number,satellite_id=j+1,command_time_s=r["command_time"],
                        start_time_s=r["start_time"],complete_time_s=r["complete_time"],
                        q_P_m_s2=float(qp[k,j]),q_D_m_s2=float(qd[k,j]),q_m_s2=float(arrays["q_m_s2"][k,j]),
                        P_fraction_of_q=float(qp[k,j]/arrays["q_m_s2"][k,j]),
                        D_fraction_of_q=float(qd[k,j]/arrays["q_m_s2"][k,j]),
                        previous_sample_time_s=float(t[k-1]) if k else None,
                        previous_sample_q_m_s2=float(arrays["q_m_s2"][k-1,j]) if k else None,
                        previous_same_satellite_completion_s=max(previous,default=None),
                        return_command_s=own_return["command_time"] if own_return else None,
                        earliest_return_allowed_s=r["complete_time"]+config.controller.effective_dwell_high_s)
            record["local_cycle_number"]=sum(z["satellite_id"]==j+1 and z["to_mode"]=="H"
                and z["command_time"]<=r["command_time"] for z in all_events)
            record["P_magnitude_exceeds_D_at_trigger"]=bool(abs(qp[k,j])>abs(qd[k,j]))
            record["return_completion_s"]=own_return["complete_time"] if own_return else None
            record["cycle_complete_by_tmax"]=bool(own_return is not None and own_return["complete_time"]<=t[-1]+1e-9)
            if previous:
                prev_own=_sample_index(t,max(previous))
                if prev_own is not None:
                    record.update(previous_cycle_end_P_m_s2=float(qp[prev_own,j]),
                                  previous_cycle_end_D_m_s2=float(qd[prev_own,j]))
            stop=_sample_index(t,own_return["command_time"]) if own_return else len(t)-1
            crossed=np.flatnonzero((arrays["q_m_s2"][k:stop+1,j]<=config.controller.threshold_off_m_s2))
            high_cross=np.flatnonzero((arrays["q_m_s2"][k:stop+1,j]<=config.controller.threshold_off_m_s2)&(arrays["mode"][k:stop+1,j]==2))
            record["first_below_off_after_start_s"]=float(t[k+crossed[0]]) if crossed.size else None
            record["first_high_state_below_off_s"]=float(t[k+high_cross[0]]) if high_cross.size else None
            if previous_end is not None and not any(z["satellite_id"]==j+1 and z["to_mode"]=="H"
                    and previous_end<z["command_time"]<r["command_time"] for z in all_events):
                prev=_sample_index(t,previous_end)
                record.update(quiet_start_s=previous_end,quiet_start_P_m_s2=float(qp[prev,j]),
                              quiet_start_D_m_s2=float(qd[prev,j]),
                              quiet_delta_P_m_s2=float(qp[k,j]-qp[prev,j]),
                              quiet_delta_D_m_s2=float(qd[k,j]-qd[prev,j]))
                # Exact kinematic identity: d(q_P)/dt=(Kp/Kd)*q_D.
                # The prediction alone additionally freezes q_D at quiet start.
                ratio=config.derived.controller_kp_s2_inv/config.derived.controller_kd_s_inv
                slope=ratio*qd[prev,j]
                wait=(config.controller.threshold_on_m_s2-qp[prev,j]-qd[prev,j])/slope if slope>0 else -1
                record["quiet_start_P_rate_m_s3"]=float(slope)
                record["constant_D_ramp_prediction_start_s"]=float(previous_end+wait) if wait>0 else None
            if own_return:
                kr=_sample_index(t,own_return["command_time"])
                if kr is not None:
                    record.update(return_q_m_s2=float(arrays["q_m_s2"][kr,j]),
                                  return_P_m_s2=float(qp[kr,j]),return_D_m_s2=float(qd[kr,j]))
            triggers.append(record)
        if complete:
            previous_end=float(end)
    differences=[]
    for left,right in zip(endpoints,endpoints[1:]):
        de=np.asarray(right["error_m"])-left["error_m"]
        dv=np.asarray(right["rate_m_s"])-left["rate_m_s"]
        da=np.asarray(right["relative_a_m"])-left["relative_a_m"]
        # Mixed-unit raw norm is meaningless: retain physical norms separately
        # and a documented dimensionless norm scaled by existing tolerances.
        se=config.recovery.position_tolerance_m;sv=config.recovery.velocity_tolerance_m_s
        scaled=float(np.linalg.norm(np.r_[de/se,dv/sv])) if se>0 and sv>0 else None
        differences.append(dict(from_wave=left["wave_number"],to_wave=right["wave_number"],
                                error_l2_difference_m=float(np.linalg.norm(de)),
                                rate_l2_difference_m_s=float(np.linalg.norm(dv)),
                                relative_a_l2_difference_m=float(np.linalg.norm(da)),
                                tolerance_scaled_state_difference=scaled))
    intervals=[]
    for satellite_id in range(1,count+1):
        own=[r for r in triggers if r["satellite_id"]==satellite_id]
        for left,right in zip(own,own[1:]):
            intervals.append(dict(satellite_id=satellite_id,from_local_cycle=left["local_cycle_number"],
                to_local_cycle=right["local_cycle_number"],start_time_s=right["start_time_s"],
                interval_s=right["start_time_s"]-left["start_time_s"]))
    blocks=[]
    # Fixed seven-day reporting bins are observation summaries, not input changes.
    for start in np.arange(0,t[-1],7*86400.):
        end=min(start+7*86400.,t[-1]);mask=(t>=start)&(t<end)
        if end==t[-1]:mask|=t==end
        k0=_sample_index(t,start);k1=_sample_index(t,end)
        observed=[r for r in all_events if start<=r["start_time"]<end or
                  (end==t[-1] and r["start_time"]==end)]
        blocks.append(dict(start_s=float(start),end_s=float(end),
            max_error_m=float(np.max(np.abs(e[mask]))),rms_error_m=float(np.sqrt(np.mean(e[mask]**2))),
            transition_starts=len(observed),l_to_h_starts=sum(r["to_mode"]=="H" for r in observed),
            extra_mean_altitude_loss_increment_m=float(extra[k1]-extra[k0])))
    grouping_sensitivity=[]
    for factor in [.5,1.,2.]:
        alternate=group_episodes(events,float(t[-1]),gap_s*factor)
        grouping_sensitivity.append(dict(gap_factor=factor,episode_count=len(alternate),
            starts_per_episode=[sum(r["start_time"]<=t[-1]+1e-9 for r in g) for g in alternate]))
    report=dict(observation_end_s=float(t[-1]),wave_gap_s=gap_s,q_decomposition_max_residual_m_s2=residual,
                E_n_definition="maximum |e| during active episode start through final return completion; quiet intervals reported separately",
                state_norm_scales=dict(error_m=config.recovery.position_tolerance_m,rate_m_s=config.recovery.velocity_tolerance_m_s),
                waves=summaries,triggers=triggers,cycle_end_states=endpoints,cycle_end_differences=differences,
                local_cycle_intervals=intervals,seven_day_blocks=blocks,grouping_sensitivity=grouping_sensitivity)
    return clean_json(report),dict(q_P_m_s2=qp,q_D_m_s2=qd,extra_mean_loss_m=extra)
