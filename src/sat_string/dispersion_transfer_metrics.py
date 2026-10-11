"""Finite-observation accuracy, operating burden and safety for transfer tests."""
import numpy as np
from .control_cost import event_activity
from .evaluation import switching_metrics, EvaluationSettings, settling_metrics
from .necessity import exceedance_metrics, improvement


def quiet_candidate(t, e, v, mode, settings):
    """Final all-LOW accepted suffix of length W; no claim beyond observation."""
    within = ((np.max(np.abs(e), axis=1)<=settings.epsilon_e_m)
        & (np.max(np.abs(v), axis=1)<=settings.epsilon_v_m_s) & np.all(mode==0,axis=1))
    suffix = np.logical_and.accumulate(within[::-1])[::-1]
    indices = np.flatnonzero(suffix & (t[-1]-t>=settings.hold_window_s))
    time = float(t[indices[0]]) if indices.size else None
    return dict(candidate_start_s=time, verification_horizon_s=float(t[-1]-time) if time is not None else None,
        epsilon_e_m=settings.epsilon_e_m,epsilon_v_m_s=settings.epsilon_v_m_s,
        hold_window_s=settings.hold_window_s,definition='all LOW and accepted final sampled suffix; finite observation only')


def evaluate_transfer(arrays, config, events, low_arrays):
    """No force changes. Same profile/time reference required for improvement/cost."""
    t=arrays['time_s']; e=arrays['gap_error_m']; v=arrays['gap_rate_m_s']; a=arrays['a_m']
    if not np.array_equal(t,low_arrays['time_s']) or not np.array_equal(a[0],low_arrays['a_m'][0]):
        raise ValueError('Control comparison needs identical time grid and initial a')
    if not np.array_equal(arrays['lambda_rad'][0],low_arrays['lambda_rad'][0]):
        raise ValueError('Initial positions differ')
    envelope=np.max(np.abs(e),axis=1);low_e=low_arrays['gap_error_m']
    switch=switching_metrics(events,config.formation.satellite_count,float(t[-1]))
    activity=event_activity(events,config.formation.satellite_count,float(t[-1]))
    extra=np.mean(low_arrays['a_m']-a,axis=1)
    signed=arrays['gap_distance_m']
    physical_signed=.5*(a[:,1:]+a[:,:-1])*np.diff(arrays['lambda_rad'],axis=1)
    reversal=np.flatnonzero(np.any(signed<=0,axis=1))
    # True same-plane Euclidean distance between ALL pairs, sampled, not collision probability.
    minimum_chord=float('inf');u=arrays['lambda_rad']
    for j in range(a.shape[1]):
        for k in range(j+1,a.shape[1]):
            square=(a[:,j]-a[:,k])**2+4*a[:,j]*a[:,k]*np.sin(.5*(u[:,j]-u[:,k]))**2
            minimum_chord=min(minimum_chord,float(np.sqrt(square.min())))
    summary=dict(max_error_m=float(envelope.max()),rms_error_m=float(np.sqrt(np.mean(e*e))),
        final_max_error_m=float(envelope[-1]),max_gap_rate_m_s=float(np.abs(v).max()),
        final_max_gap_rate_m_s=float(np.abs(v[-1]).max()),
        max_error_percent_of_100km=float(envelope.max()/config.formation.target_spacing_m*100),
        gap_rms_m=np.sqrt(np.mean(e*e,axis=0)).tolist(),
        rms_definition='sqrt(mean(e_i(t)^2)) over all saved time samples and fixed-label adjacent gaps',
        transitions=switch['n_transition_starts'],completions=switch['n_transition_completions'],
        l_to_h=sum(x['to_mode']=='H' and x['start_time']<=t[-1] for x in events),
        h_to_l=sum(x['to_mode']=='L' and x['start_time']<=t[-1] for x in events),
        completed_cycles=sum(switch['lhl_cycles_per_satellite']),
        switches_per_satellite=switch['switches_per_satellite'],cycles_per_satellite=switch['lhl_cycles_per_satellite'],
        last_switch_time_s=switch['last_switch_time_s'],
        high_stable_total_s=activity['high_stable_total_s'],high_equivalent_total_s=activity['high_equivalent_total_s'],
        high_stable_s_per_satellite=activity['high_stable_s_by_satellite'],
        high_equivalent_s_per_satellite=activity['high_equivalent_s_by_satellite'],
        max_simultaneous_stable_HIGH=int(np.sum(arrays['mode']==2,axis=1).max()),
        max_simultaneous_nonLOW=int(np.sum(arrays['mode']!=0,axis=1).max()),
        final_modes=arrays['mode'][-1].tolist(),
        extra_mean_altitude_loss_m=float(extra[-1]),
        mean_altitude_loss_m=float(np.mean(a[0]-a[-1])),
        final_altitude_m=(a[-1]-config.constants.earth_radius_m).tolist(),
        final_mean_altitude_m=float(np.mean(a[-1]-config.constants.earth_radius_m)),
        final_altitude_span_m=float(np.ptp(a[-1])),max_altitude_span_m=float(np.ptp(a,axis=1).max()),
        final_allLOW_altitude_span_m=float(np.ptp(low_arrays['a_m'][-1])),
        extra_loss_per_satellite_m=(low_arrays['a_m'][-1]-a[-1]).tolist(),
        min_signed_reference_gap_m=float(signed.min()),min_signed_current_radius_arc_m=float(physical_signed.min()),
        min_physical_pair_chord_m=minimum_chord,
        first_label_order_reversal_s=float(t[reversal[0]]) if reversal.size else None,
        labelled_precision_valid_through_s=float(t[reversal[0]]) if reversal.size else float(t[-1]),
        max_current_radius_arc_error_m=float(np.max(np.abs(physical_signed-config.formation.target_spacing_m))),
        allLOW_order_reversal_s=None)
    low_bad=np.flatnonzero(np.any(low_arrays['gap_distance_m']<=0,axis=1))
    summary['allLOW_order_reversal_s']=float(t[low_bad[0]]) if low_bad.size else None
    summary['full_horizon_precision_reduction_valid']=not reversal.size and not low_bad.size
    cut=int(low_bad[0]) if low_bad.size else len(t)
    summary['matched_pre_allLOW_reversal_improvement']=dict(
        last_valid_sample_s=float(t[cut-1]) if cut else None,
        **improvement(float(np.abs(e[:cut]).max()),float(np.abs(low_e[:cut]).max()))) if cut else None
    pre=np.arange(len(t))<(int(reversal[0]) if reversal.size else len(t))
    summary['max_error_before_order_reversal_m']=float(envelope[pre].max()) if np.any(pre) else None
    summary['precision_interpretation']='fixed-label diagnostic after reversal; not formation-keeping precision' if reversal.size else 'fixed-label adjacent precision while order retained'
    summary['error_improvement']=improvement(summary['max_error_m'],float(np.abs(low_e).max()))
    summary['RMS_improvement']=improvement(summary['rms_error_m'],float(np.sqrt(np.mean(low_e**2))))
    summary['exceedances']=[exceedance_metrics(t,envelope,x) for x in (100,500,1000)]
    settings=EvaluationSettings.from_config(config)
    summary['quiet_candidate']=quiet_candidate(t,e,v,arrays['mode'],settings)
    summary['settling']=settling_metrics(t,e,v,switch,settings,evaluation_start_s=0)
    segments=[]
    for start,stop in [(0,7),(7,30),(30,60)]:
        lo=start*86400.;hi=min(stop*86400.,t[-1])
        if hi<=lo: continue
        use=(t>=lo)&(t<=hi);last=hi==t[-1]
        count=sum(lo<=x['start_time'] and (x['start_time']<=hi if last else x['start_time']<hi) for x in events)
        left=event_activity(events,21,lo);right=event_activity(events,21,hi)
        i0=int(np.searchsorted(t,lo));i1=int(np.searchsorted(t,hi));days=(hi-lo)/86400
        segments.append(dict(start_day=start,end_day=hi/86400,transitions=count,transitions_per_day=count/days,
            high_stable_satellite_s_per_day=(right['high_stable_total_s']-left['high_stable_total_s'])/days,
            high_equivalent_satellite_s_per_day=(right['high_equivalent_total_s']-left['high_equivalent_total_s'])/days,
            error_rms_m=float(np.sqrt(np.mean(e[use]**2))),max_error_m=float(np.abs(e[use]).max()),
            extra_mean_altitude_loss_per_day_m=float((extra[i1]-extra[i0])/days),
            mean_natural_plus_control_altitude_loss_per_day_m=float(np.mean(a[i0]-a[i1])/days)))
    summary['segments']=segments
    return summary
