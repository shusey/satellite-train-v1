"""Post-cycle residual states and falsifiable constant-D forecasts; no dynamics changes."""
from __future__ import annotations
from math import ceil
import numpy as np
from .control_cost import clean_json
from .long_horizon import decompose_q,group_episodes
from .continuation import restore_attitude


def adjacent_gaps(satellites,count):
    """Zero-based edges incident on one-based satellite IDs; v lives on edges."""
    return sorted({g for j in satellites for g in (j-2,j-1) if 0<=g<count-1})


def cycle_pairs(events,end_s):
    """Local cycle = one satellite's L->H then H->L; preserve incomplete rises."""
    own={}
    for event in sorted(events,key=lambda z:(z["command_time"],z["satellite_id"])):
        if event["command_time"]<=end_s+1e-9:
            own.setdefault(event["satellite_id"],[]).append(event)
    result=[]
    for sat,records in own.items():
        active=None;number=0
        for event in records:
            if event["to_mode"]=="H":
                if active is not None:raise ValueError("Missing return before rise")
                number+=1;active=event
            elif event["to_mode"]=="L":
                if active is None:raise ValueError("Return without rise")
                result.append(dict(satellite_id=sat,local_cycle_number=number,rise=active,return_event=event))
                active=None
            else:raise ValueError("Unknown target")
        if active is not None:
            result.append(dict(satellite_id=sat,local_cycle_number=number,rise=active,return_event=None))
    return sorted(result,key=lambda z:(z["rise"]["command_time"],z["satellite_id"]))


def predict_crossing(qp,qd,time_s,allowed_s,kp,kd,h_on,period_s):
    """q=qP0+qD0+(Kp/Kd)qD0*(t-t0); command includes LOW dwell/grid.

    This extrapolation assumes constant D and no intervening actuation. No
    upward crossing is predicted for nonpositive D below the HIGH threshold.
    """
    if not all(np.isfinite(x) for x in [qp,qd,time_s,allowed_s,kp,kd,h_on,period_s]) or min(kp,kd,period_s)<=0:
        raise ValueError("Finite state and positive gains/period required")
    slope=kp/kd*qd
    if qp+qd>=h_on:
        crossing=time_s;status="already_at_or_above_on"
    elif slope<=0:
        return dict(status="no_upward_crossing_constant_D",predicted_crossing_s=None,
                    predicted_command_s=None,predicted_wait_s=None,q_P_rate_m_s3=float(slope))
    else:
        crossing=time_s+(h_on-qp-qd)/slope;status="future_upward_crossing"
    command=ceil((max(crossing,allowed_s)-1e-9)/period_s)*period_s
    if qp+qd+slope*(command-time_s)<h_on-1e-20:
        return dict(status="falls_below_on_before_allowed",predicted_crossing_s=float(crossing),
                    predicted_command_s=None,predicted_wait_s=None,q_P_rate_m_s3=float(slope))
    return dict(status=status,predicted_crossing_s=float(crossing),predicted_command_s=float(command),
                predicted_wait_s=float(command-time_s),q_P_rate_m_s3=float(slope))


def node_component(edges,j,gain):
    if j==0:return gain*edges[:,0]
    if j==edges.shape[1]:return -gain*edges[:,-1]
    return gain*(edges[:,j]-edges[:,j-1])


def quiet_packets(pairs,end_s):
    """Merge overlapping command-to-LOW-completion intervals, with no gap cutoff.

    A completed packet ends at an all-LOW, no-pending boundary. This retains
    short quiet gaps that the older temporal response grouping may merge.
    """
    packets=[];active=[];finish=None;complete=True
    for pair in sorted(pairs,key=lambda p:(p["rise"]["command_time"],p["satellite_id"])):
        rise=pair["rise"];back=pair["return_event"];start=rise["command_time"]
        done=back is not None and back["complete_time"]<=end_s+1e-9
        stop=back["complete_time"] if done else end_s
        if active and start>finish+1e-9:
            packets.append((active,finish,complete));active=[];complete=True
        active.append(pair);finish=max(finish if len(active)>1 else start,stop);complete=complete and done
    if active:packets.append((active,finish,complete))
    result=[]
    for number,(records,finish,complete) in enumerate(packets,1):
        start=min(p["rise"]["command_time"] for p in records)
        result.append(dict(wave_number=number,complete=complete,start_time_s=start,
            end_time_s=finish if complete else None,
            first_satellite_ids=sorted({p["satellite_id"] for p in records if p["rise"]["command_time"]==start}),
            l_to_h_satellite_ids=sorted({p["satellite_id"] for p in records}),transition_starts=2*len(records)))
    return result


def end_state(arrays,time_s,participants,seeds,config):
    """First saved sample at/after completion, with explicit sampling offset.

    V is a gap velocity, not the orbital velocity of an individual satellite.
    Seed/participant V uses the union of incident adjacent gaps.
    """
    t=arrays["time_s"];k=int(np.searchsorted(t,time_s-1e-9))
    if k>=len(t):return None,None
    e=arrays["gap_error_m"][k:k+1];v=arrays["gap_rate_m_s"][k:k+1]
    p,d=decompose_q(e,v,config.derived.controller_kp_s2_inv,config.derived.controller_kd_s_inv)
    q=arrays["q_m_s2"][k]
    if np.max(np.abs(p[0]+d[0]-q))>1e-20:raise ValueError("Stored q disagrees with decomposition")
    count=len(q);gp=adjacent_gaps(participants,count);gs=adjacent_gaps(seeds,count)
    summary=dict(completion_time_s=float(time_s),state_sample_time_s=float(t[k]),
        sample_offset_s=float(t[k]-time_s),source_sample_index=k,
        participant_satellite_ids=list(participants),seed_satellite_ids=list(seeds),
        participant_gap_ids=[g+1 for g in gp],seed_gap_ids=[g+1 for g in gs],
        V_res_all_m_s=float(np.max(np.abs(v))),
        V_res_participants_m_s=float(np.max(np.abs(v[0,gp]))),
        V_res_seeds_m_s=float(np.max(np.abs(v[0,gs]))),
        Q_D_res_all_m_s2=float(np.max(np.abs(d))),
        Q_D_res_participants_m_s2=float(np.max(np.abs(d[0,np.asarray(participants)-1]))),
        Q_D_res_seeds_m_s2=float(np.max(np.abs(d[0,np.asarray(seeds)-1]))),
        all_low_at_sample=bool(np.all(arrays["mode"][k]==0)),
        nonlow_satellite_ids=(np.flatnonzero(arrays["mode"][k]!=0)+1).tolist())
    state=dict(time_s=float(t[k]),gap_error_m=e[0],gap_rate_m_s=v[0],
               q_P_m_s2=p[0],q_D_m_s2=d[0],q_m_s2=q,mode=arrays["mode"][k])
    return summary,state


def _compare(arrays,events,summary,state,j,prediction,actual_s,config,relative_D_tolerance):
    t=arrays["time_s"];t0=summary["state_sample_time_s"];end=float(actual_s if actual_s is not None else t[-1])
    k0=summary["source_sample_index"];k1=int(np.searchsorted(t,end+1e-9))
    td=t[k0:k1]
    kp,kd=config.derived.controller_kp_s2_inv,config.derived.controller_kd_s_inv
    ds=node_component(arrays["gap_rate_m_s"][k0:k1],j,kd)
    ps=node_component(arrays["gap_error_m"][k0:k1],j,kp)
    d0=state["q_D_m_s2"][j]
    variation=float(np.max(np.abs(ds-d0))/abs(d0)) if d0!=0 else None
    crossed=np.flatnonzero(arrays["q_m_s2"][k0:k1,j]>=config.controller.threshold_on_m_s2)
    intervening=[z for z in events if t0+1e-9<z["command_time"]<end-1e-9]
    pending=[z for z in events if z["command_time"]<=t0+1e-9<z["complete_time"]]
    reasons=[]
    if not summary["all_low_at_sample"]:reasons.append("other_satellites_not_LOW_at_end")
    if pending:reasons.append("pending_or_active_transition_at_end")
    if intervening:reasons.append("intervening_commands")
    if variation is None or variation>relative_D_tolerance:reasons.append("D_not_constant_within_analysis_tolerance")
    if prediction["predicted_command_s"] is None:reasons.append("no_finite_upward_forecast")
    valid=not reasons
    predicted=prediction["predicted_command_s"]
    error=float(actual_s-predicted) if actual_s is not None and predicted is not None else None
    wait=float(actual_s-t0) if actual_s is not None else None
    return dict(actual_next_command_s=actual_s,actual_next_wait_s=wait,
        actual_threshold_first_sample_s=float(td[crossed[0]]) if crossed.size else None,
        observation_window_end_s=end,right_censored=actual_s is None,
        max_relative_D_change=variation,intervening_command_count=len(intervening),
        prediction_valid_on_observed_window=valid,invalid_reasons=reasons,
        prediction_error_actual_minus_predicted_s=error,
        relative_prediction_error=error/wait if error is not None and wait and wait>0 else None,
        delta_q_P_m_s2=float(ps[-1]-ps[0]),delta_q_D_m_s2=float(ds[-1]-ds[0]),
        q_P_change_from_integrated_D_m_s2=float(kp/kd*np.trapezoid(ds,td)),
        validated_prediction_error_s=error if valid else None)


def _off_activity(arrays,rise,back,config):
    t=arrays["time_s"];j=rise["satellite_id"]-1
    mask=(t>=rise["start_time"]-1e-9)&(t<=back["command_time"]+1e-9)
    tt=t[mask];q=arrays["q_m_s2"][mask,j];m=arrays["mode"][mask,j]
    below=np.flatnonzero(q<=config.controller.threshold_off_m_s2)
    high=np.flatnonzero((q<=config.controller.threshold_off_m_s2)&(m==2))
    first=float(tt[below[0]]) if below.size else None
    allowed=rise["complete_time"]+config.controller.dwell_time_s
    earliest=ceil((allowed-1e-9)/config.controller.control_period_s)*config.controller.control_period_s
    hold=back["command_time"]-rise["complete_time"]
    high_stable=back["start_time"]-rise["complete_time"]
    equivalent=high_stable+.5*(rise["complete_time"]-rise["start_time"]+back["complete_time"]-back["start_time"])
    high_k=int(np.searchsorted(t,rise["complete_time"]-1e-9))
    low_k=int(np.searchsorted(t,back["command_time"]-1e-9))
    return dict(high_command_s=rise["command_time"],high_slew_start_s=rise["start_time"],
        high_complete_s=rise["complete_time"],low_command_s=back["command_time"],
        low_slew_start_s=back["start_time"],low_complete_s=back["complete_time"],
        high_complete_to_low_command_s=float(hold),high_stable_time_s=float(high_stable),
        high_equivalent_time_s=float(equivalent),first_q_below_off_s=first,
        first_HIGH_q_below_off_s=float(tt[high[0]]) if high.size else None,
        off_cross_to_low_command_s=float(back["command_time"]-first) if first is not None else None,
        earliest_low_command_s=float(earliest),
        q_at_HIGH_complete_m_s2=float(arrays["q_m_s2"][high_k,j]),
        q_at_LOW_command_m_s2=float(arrays["q_m_s2"][low_k,j]),
        return_driver="dwell_limited" if abs(back["command_time"]-earliest)<=1e-9 else "threshold_limited_after_dwell")


def residual_analysis(arrays,config,events,*,groups=None,gap_s=None,relative_D_tolerance=.1,low_arrays=None):
    """Local endpoint states plus complete-response endpoints with separate forecasts.

    relative_D_tolerance is an explicit postprocessing quality setting, not a
    controller or physical parameter. Raw errors survive even when invalid.
    """
    if not np.isfinite(relative_D_tolerance) or relative_D_tolerance<0:raise ValueError("Invalid analysis tolerance")
    end=float(arrays["time_s"][-1]);kp=config.derived.controller_kp_s2_inv;kd=config.derived.controller_kd_s_inv
    frames=[];local=[];global_rows=[];quiet_rows=[];incomplete=[]
    def snapshot(time,participants,seeds):
        summary,state=end_state(arrays,time,participants,seeds,config)
        if state is not None:
            summary["end_state_row"]=len(frames);frames.append(state)
        return summary,state
    pairs=cycle_pairs(events,end)
    for pair in pairs:
        rise,back=pair["rise"],pair["return_event"];sat=pair["satellite_id"];j=sat-1
        if back is None or back["complete_time"]>end+1e-9:
            incomplete.append(dict(satellite_id=sat,local_cycle_number=pair["local_cycle_number"],start_time_s=rise["start_time"]))
            continue
        row,state=snapshot(back["complete_time"],[sat],[sat])
        row.update(satellite_id=sat,local_cycle_number=pair["local_cycle_number"],**_off_activity(arrays,rise,back,config))
        for key in ["q_P_m_s2","q_D_m_s2","q_m_s2"]:row["own_"+key]=float(state[key][j])
        first_k=int(np.searchsorted(arrays["time_s"],rise["command_time"]-1e-9))
        last_k=row["source_sample_index"];incident=adjacent_gaps([sat],config.formation.satellite_count)
        active_e=arrays["gap_error_m"][first_k:last_k+1]
        row["E_n_active_all_m"]=float(np.max(np.abs(active_e)))
        row["E_n_active_incident_m"]=float(np.max(np.abs(active_e[:,incident])))
        row["active_rms_all_m"]=float(np.sqrt(np.mean(active_e**2)))
        if low_arrays is not None:
            ix=[max(0,first_k-1),last_k]
            cost=low_arrays["a_m"][ix,j]-arrays["a_m"][ix,j]
            row["own_extra_loss_active_increment_m"]=float(cost[1]-cost[0])
        nxt=next((p["rise"] for p in pairs if p["satellite_id"]==sat and p["local_cycle_number"]==pair["local_cycle_number"]+1),None)
        actual=nxt["command_time"] if nxt else None
        pred=predict_crossing(state["q_P_m_s2"][j],state["q_D_m_s2"][j],state["time_s"],
            back["complete_time"]+config.controller.dwell_time_s,kp,kd,config.controller.threshold_on_m_s2,config.controller.control_period_s)
        row.update(prediction=pred,**_compare(arrays,events,row,state,j,pred,actual,config,relative_D_tolerance))
        row["next_same_satellite_start_interval_s"]=actual-rise["command_time"] if actual is not None else None
        local.append(row)
    if groups is None:
        if gap_s is None:raise ValueError("Response grouping or explicit reference gap required")
        groups=[]
        for number,records in enumerate(group_episodes(events,end,gap_s),1):
            rises=[z for z in records if z["to_mode"]=="H"]
            returns=[z for z in records if z["to_mode"]=="L"]
            complete=len(rises)==len(returns) and all(z["complete_time"]<=end+1e-9 for z in records)
            start=min(z["start_time"] for z in rises)
            groups.append(dict(wave_number=number,complete=complete,start_time_s=start,
                end_time_s=max(z["complete_time"] for z in records) if complete else None,
                first_satellite_ids=[z["satellite_id"] for z in rises if z["start_time"]==start],
                l_to_h_satellite_ids=sorted({z["satellite_id"] for z in rises}),transition_starts=len(records)))
    tagged=[("response_group",g) for g in groups]+[("quiet_packet",g) for g in quiet_packets(pairs,end)]
    for scope,group in tagged:
        if not group["complete"] or group["end_time_s"]>end:continue
        participants=group["l_to_h_satellite_ids"];seeds=group["first_satellite_ids"]
        row,state=snapshot(group["end_time_s"],participants,seeds)
        row.update(endpoint_scope=scope,response_group_number=group["wave_number"],group_start_s=group["start_time_s"],
                   group_transition_starts=group["transition_starts"])
        attitude=restore_attitude(config,events,state["time_s"])
        candidates=[]
        for j in range(config.formation.satellite_count):
            if attitude.mode[j]!=0 or attitude.pending_target[j]!=-1:continue
            pred=predict_crossing(state["q_P_m_s2"][j],state["q_D_m_s2"][j],state["time_s"],
                attitude.next_allowed_time[j],kp,kd,config.controller.threshold_on_m_s2,config.controller.control_period_s)
            if pred["predicted_command_s"] is not None:candidates.append((pred["predicted_command_s"],j,pred))
        future=[z for z in events if z["to_mode"]=="H" and state["time_s"]+1e-9<z["command_time"]<=end+1e-9]
        nxt=min(future,key=lambda z:(z["command_time"],z["satellite_id"])) if future else None
        actual=nxt["command_time"] if nxt else None;actual_j=nxt["satellite_id"]-1 if nxt else None
        if candidates:
            _,j,pred=min(candidates,key=lambda z:(z[0],z[1]));row["predicted_satellite_id"]=j+1
        else:
            j=actual_j if actual_j is not None else int(np.argmax(state["q_D_m_s2"]))
            pred=predict_crossing(state["q_P_m_s2"][j],state["q_D_m_s2"][j],state["time_s"],
                attitude.next_allowed_time[j],kp,kd,config.controller.threshold_on_m_s2,config.controller.control_period_s)
            row["predicted_satellite_id"]=None
        row.update(prediction=pred,actual_next_satellite_id=actual_j+1 if actual_j is not None else None,
                   prediction_target_matches_actual=bool(actual_j==j) if actual_j is not None else None)
        row.update(_compare(arrays,events,row,state,j,pred,actual,config,relative_D_tolerance))
        for key in ["q_P_m_s2","q_D_m_s2","q_m_s2"]:
            row["predicted_satellite_"+key]=float(state[key][j])
            row["actual_next_satellite_"+key]=float(state[key][actual_j]) if actual_j is not None else None
        start_k=int(np.searchsorted(arrays["time_s"],group["start_time_s"]))
        end_k=row["source_sample_index"]
        row["E_n_active_max_error_m"]=float(np.max(np.abs(arrays["gap_error_m"][start_k:end_k+1])))
        if low_arrays is not None:
            extra=np.mean(low_arrays["a_m"][[max(0,start_k-1),end_k]]-arrays["a_m"][[max(0,start_k-1),end_k]],axis=1)
            row["extra_loss_increment_m"]=float(extra[1]-extra[0])
        (global_rows if scope=="response_group" else quiet_rows).append(row)
    packed={key:np.asarray([s[key] for s in frames]) for key in frames[0]} if frames else {}
    summary=dict(schema_version=1,observation_end_s=end,dwell_time_s=config.controller.dwell_time_s,
        relative_D_tolerance=relative_D_tolerance,local_completed_cycles=len(local),incomplete_cycles=incomplete,
        complete_response_groups=len(global_rows),
        local_valid_finite_forecasts=sum(r["prediction_valid_on_observed_window"] and not r["right_censored"] for r in local),
        group_valid_finite_forecasts=sum(r["prediction_valid_on_observed_window"] and not r["right_censored"] for r in global_rows),
        complete_quiet_packets=len(quiet_rows),
        quiet_valid_finite_forecasts=sum(r["prediction_valid_on_observed_window"] and not r["right_censored"] for r in quiet_rows),
        dwell_limited_returns=sum(r["return_driver"]=="dwell_limited" for r in local),
        threshold_limited_returns=sum(r["return_driver"]!="dwell_limited" for r in local),
        high_hold_times_s=sorted({r["high_complete_to_low_command_s"] for r in local}),
        local_definition="one satellite L-H-L, sampled at first saved time >= LOW completion; neighbor activity flagged",
        group_definition="existing temporal response group, all returns completed; future events do not define residual values",
        forecast_definition="constant D, then LOW dwell and control-grid gate; raw invalid forecasts retained with reasons")
    return clean_json(dict(summary=summary,local_cycles=local,group_cycles=global_rows,quiet_cycles=quiet_rows)),packed
