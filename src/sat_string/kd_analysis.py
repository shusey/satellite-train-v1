"""Offline Kd comparison: physical velocity distribution and threshold contributions."""
from __future__ import annotations
import numpy as np
from .long_horizon import decompose_q,group_episodes,_fit_outer_branches
from .residual import cycle_pairs,quiet_packets,predict_crossing
from .control_cost import clean_json
from .control_cost import event_activity


def velocity_distribution(error,rate):
    """All-gap statistics; e*v>0 increases squared error, zero is neutral."""
    e=np.asarray(error,dtype=float);v=np.asarray(rate,dtype=float)
    if e.ndim!=1 or e.shape!=v.shape or not e.size or not np.all(np.isfinite(e+v)):
        raise ValueError("Finite matching nonempty gap vectors required")
    ev=e*v
    magnitude=float(np.sum(np.abs(ev)))
    return dict(V_res_all_m_s=float(np.max(np.abs(v))),
        V_res_rms_all_m_s=float(np.sqrt(np.mean(v*v))),
        V_res_median_abs_all_m_s=float(np.median(np.abs(v))),
        error_rms_at_end_m=float(np.sqrt(np.mean(e*e))),
        worsening_gap_fraction=float(np.mean(ev>0)),recovering_gap_fraction=float(np.mean(ev<0)),
        neutral_gap_fraction=float(np.mean(ev==0)),
        error_squared_mean_rate_m2_s=float(2*np.mean(ev)),
        worsening_weight_fraction=float(np.sum(np.maximum(ev,0))/magnitude) if magnitude else None,
        quadrant_counts=dict(e_positive_v_positive=int(np.sum((e>0)&(v>0))),
            e_positive_v_negative=int(np.sum((e>0)&(v<0))),
            e_negative_v_negative=int(np.sum((e<0)&(v<0))),
            e_negative_v_positive=int(np.sum((e<0)&(v>0))),zero_product=int(np.sum(ev==0))))


def contribution(qp,qd,dominance_ratio=2.):
    """Magnitude dominance (2:1 default), with signs/cancellation kept separately."""
    if dominance_ratio<=1 or not np.isfinite(dominance_ratio):raise ValueError("Dominance ratio must exceed one")
    magnitude=abs(qp)+abs(qd)
    if magnitude==0:kind="zero"
    elif abs(qp)>=dominance_ratio*abs(qd):kind="P_dominant"
    elif abs(qd)>=dominance_ratio*abs(qp):kind="D_dominant"
    else:kind="mixed"
    return dict(q_P_m_s2=float(qp),q_D_m_s2=float(qd),q_m_s2=float(qp+qd),
        D_magnitude_fraction=float(abs(qd)/magnitude) if magnitude else None,
        contribution_class=kind,opposed_signs=bool(qp*qd<0),
        cancellation_ratio=float(abs(qp+qd)/magnitude) if magnitude else None)


def match_all_low(config,low_config):
    """Only Kd may differ; disabled-controller status is checked by the runner."""
    left,right=config.to_dict(),low_config.to_dict()
    for item in [left,right]:
        item["controller"].pop("kd_s_inv",None)
        item["derived"].pop("controller_kd_s_inv")
    return left==right


def threshold_events(arrays,config,events,dominance_ratio=2.):
    """Final uninterrupted above-on segment in the current LOW interval vs command.

    Already above at LOW entry is not labeled a new upward crossing. Saved
    samples are post-event; q itself is computed before the command changes area.
    """
    t=arrays["time_s"];on=config.controller.threshold_on_m_s2
    kp,kd=config.derived.controller_kp_s2_inv,config.derived.controller_kd_s_inv
    packets=quiet_packets(cycle_pairs(events,float(t[-1])),float(t[-1]))
    packet_starts={x["start_time_s"] for x in packets}
    previous_low={};rows=[]
    for event in sorted(events,key=lambda x:(x["command_time"],x["satellite_id"])):
        sat=event["satellite_id"];j=sat-1
        if event["to_mode"]=="L":
            previous_low[sat]=event["complete_time"];continue
        time=event["command_time"]
        if time>t[-1]+1e-9:continue
        start=previous_low.get(sat,0.)
        k0=int(np.searchsorted(t,start-1e-9));k=int(np.searchsorted(t,time-1e-9))
        if k>=len(t) or abs(t[k]-time)>1e-8:raise ValueError("Command must be on saved grid")
        q=arrays["q_m_s2"][k0:k+1,j]
        if q[-1]<on:raise ValueError("HIGH command below threshold")
        below=np.flatnonzero(q<on);cross=k0+int(below[-1])+1 if below.size else k0
        row=dict(satellite_id=sat,command_time_s=float(time),LOW_entry_time_s=float(start),
            quiet_packet_onset=bool(time in packet_starts),
            sampled_crossing_time_s=float(t[cross]),command_minus_crossing_s=float(time-t[cross]),
            crossing_kind="upward_from_below" if below.size else "already_above_at_LOW_entry")
        for label,index in [("command",k),("crossing",cross)]:
            p,d=decompose_q(arrays["gap_error_m"][index:index+1],arrays["gap_rate_m_s"][index:index+1],kp,kd)
            row.update({label+"_"+key:value for key,value in contribution(p[0,j],d[0,j],dominance_ratio).items()})
        rows.append(row)
    return rows


def fixed_state_forecasts(error,rate,time_s,allowed,config,factors=(.3,1.,3.)):
    """Counterfactual signal evaluation only: SAME physical state, no simulation.

    P slope is Kp*B*v independent of Kd; the initial qD offset still changes
    the threshold margin even when the physical velocity is held fixed.
    """
    kp=config.derived.controller_kp_s2_inv;kd0=2/21600.;result=[]
    for factor in factors:
        kd=kd0*factor;p,d=decompose_q(np.asarray(error)[None,:],np.asarray(rate)[None,:],kp,kd)
        candidates=[]
        for j in range(len(p[0])):
            pred=predict_crossing(p[0,j],d[0,j],time_s,allowed[j],kp,kd,
                config.controller.threshold_on_m_s2,config.controller.control_period_s)
            if pred["predicted_command_s"] is not None:candidates.append((pred["predicted_command_s"],j,pred))
        if candidates:
            _,j,pred=min(candidates,key=lambda x:(x[0],x[1]))
            result.append(dict(kd_factor=factor,satellite_id=j+1,q_P_m_s2=float(p[0,j]),q_D_m_s2=float(d[0,j]),**pred))
        else:result.append(dict(kd_factor=factor,satellite_id=None,predicted_command_s=None))
    return result


def return_gates(arrays,config,events):
    """Include unfinished HIGH pulses, so delayed/censored returns are visible."""
    t=arrays["time_s"];off=config.controller.threshold_off_m_s2;period=config.controller.control_period_s
    rows=[]
    for pair in cycle_pairs(events,float(t[-1])):
        rise=pair["rise"];back=pair["return_event"];j=pair["satellite_id"]-1
        gate=np.ceil((rise["complete_time"]+config.controller.dwell_time_s-1e-9)/period)*period
        k=int(np.searchsorted(t,gate-1e-9));observed=k<len(t)
        rows.append(dict(satellite_id=j+1,local_cycle_number=pair["local_cycle_number"],
            high_complete_s=rise["complete_time"],earliest_return_gate_s=float(gate),
            gate_observed=observed,q_at_return_gate_m_s2=float(arrays["q_m_s2"][k,j]) if observed else None,
            off_not_met_at_gate=bool(arrays["q_m_s2"][k,j]>off) if observed else None,
            low_command_s=back["command_time"] if back else None,
            low_complete_s=back["complete_time"] if back else None,
            observed_hold_s=float(back["command_time"]-rise["complete_time"]) if back else None,
            return_right_censored=back is None,
            completed_by_tmax=bool(back is not None and back["complete_time"]<=t[-1]+1e-9)))
    return rows


def enhance_residual(report,states,config):
    """Add statistics to existing endpoint rows; full signed vectors stay in NPZ."""
    for section in ["local_cycles","group_cycles","quiet_cycles"]:
        for row in report[section]:
            k=row["end_state_row"]
            row.update(velocity_distribution(states["gap_error_m"][k],states["gap_rate_m_s"][k]))
            row["Q_D_res_div_Kd_m_s"]=row["Q_D_res_all_m_s2"]/config.derived.controller_kd_s_inv
    report["summary"].update(Kp_s2_inv=config.derived.controller_kp_s2_inv,
        Kd_s_inv=config.derived.controller_kd_s_inv,
        velocity_distribution_definition="all gaps at endpoint; strict e*v signs; no engineering zero tolerance")
    return report


def first_day_metrics(arrays,config,events):
    """Same first-day window for all gains, independent of temporal grouping."""
    t=arrays["time_s"];end=min(86400.,float(t[-1]));e=arrays["gap_error_m"][t<=end+1e-9]
    return clean_json(dict(end_time_s=end,max_error_m=float(np.max(np.abs(e))),
        rms_error_m=float(np.sqrt(np.mean(e**2))),final_max_error_m=float(np.max(np.abs(e[-1]))),
        transition_starts=sum(x["start_time"]<=end+1e-9 for x in events),
        HIGH_activity=event_activity(events,config.formation.satellite_count,end)))


def first_response(arrays,config,events,gap_s,report,states):
    """Temporal first group and each satellite's first local response, kept distinct."""
    t=arrays["time_s"];groups=group_episodes(events,float(t[-1]),gap_s)
    first=[pair for pair in cycle_pairs(events,float(t[-1])) if pair["local_cycle_number"]==1]
    starts=[pair["rise"] for pair in first]
    if not starts:return dict(first_command_time_s=None,first_satellite_ids=[],definition="no observed response")
    onset=min(x["command_time"] for x in starts);records=groups[0]
    rises=[x for x in records if x["to_mode"]=="H"];backs=[x for x in records if x["to_mode"]=="L"]
    complete=len(rises)==len(backs) and all(x["complete_time"]<=t[-1] for x in records)
    finish=max(x["complete_time"] for x in records) if complete else float(t[-1])
    k0=int(np.searchsorted(t,onset-1e-9));k1=int(np.searchsorted(t,finish+1e-9))
    local=[dict(x) for x in report["local_cycles"] if x["local_cycle_number"]==1]
    # A fixed-gap group may contain multiple cycles; fit FIRST per-satellite arrivals separately.
    fits=_fit_outer_branches(starts,config,float(t[-1]))
    return clean_json(dict(first_command_time_s=onset,
        first_satellite_ids=sorted(x["satellite_id"] for x in starts if x["command_time"]==onset),
        first_arrivals_s_by_satellite={str(x["satellite_id"]):x["start_time"] for x in starts},
        first_arrival_outer_fits=fits,first_local_cycles=local,first_24h=first_day_metrics(arrays,config,events),
        first_group=dict(complete=complete,end_time_s=finish if complete else None,
            observed_end_s=finish,transition_starts=len(records),
            participant_satellite_ids=sorted({x["satellite_id"] for x in rises}),
            contains_repeated_local_cycles=len(rises)>len({x["satellite_id"] for x in rises}),
            max_error_m=float(np.max(np.abs(arrays["gap_error_m"][k0:k1]))),
            rms_error_m=float(np.sqrt(np.mean(arrays["gap_error_m"][k0:k1]**2)))),
        first_quiet_endpoint=report["quiet_cycles"][0] if report["quiet_cycles"] else None,
        definition="frozen-gap first group may merge responses; first local cycles and first quiet packet are separate"))
