"""Kp-only offline diagnostics; reuse physical velocity and hybrid forecasts."""
from __future__ import annotations
import numpy as np
from .long_horizon import decompose_q
from .residual import predict_crossing,cycle_pairs


def match_all_low(config,low_config):
    """Kp only may differ: a disabled controller cannot affect physical motion."""
    left,right=config.to_dict(),low_config.to_dict()
    for item in [left,right]:
        item['controller'].pop('kp_s2_inv',None)
        item['derived'].pop('controller_kp_s2_inv')
    return left==right


def activation_state(error,rate,satellite_id,config):
    """Physical errors at a trigger; incident gaps are j-1 and j in 0-based IDs.

    Full signed vectors remain in the original run. B is the existing oriented
    open-path incidence matrix, not a new position or velocity definition.
    """
    e=np.asarray(error,dtype=float);v=np.asarray(rate,dtype=float)
    n=config.formation.satellite_count
    if e.shape!=(n-1,) or v.shape!=e.shape or not np.all(np.isfinite(e+v)):
        raise ValueError('Finite formation gap vectors required')
    if not 1<=satellite_id<=n:raise ValueError('Satellite ID outside formation')
    j=satellite_id-1;indices=[k for k in [j-1,j] if 0<=k<len(e)]
    kp,kd=config.derived.controller_kp_s2_inv,config.derived.controller_kd_s_inv
    p,d=decompose_q(e[None,:],v[None,:],kp,kd)
    return dict(max_abs_error_m=float(np.max(np.abs(e))),
        rms_error_m=float(np.sqrt(np.mean(e*e))),
        incident_max_abs_error_m=float(np.max(np.abs(e[indices]))),
        incident_gap_ids=[k+1 for k in indices],
        incident_signed_error_m=e[indices].tolist(),incident_signed_rate_m_s=v[indices].tolist(),
        B_error_m=float(p[0,j]/kp),B_rate_m_s=float(d[0,j]/kd),
        q_P_rate_m_s3=float(kp/kd*d[0,j]),
        on_margin_m_s2=float(config.controller.threshold_on_m_s2-p[0,j]-d[0,j]))


def enhance_activations(rows,arrays,config):
    """Keep sampled threshold crossing and actual HIGH command distinct."""
    t=arrays['time_s']
    for row in rows:
        for label,key in [('command','command_time_s'),('crossing','sampled_crossing_time_s')]:
            k=int(np.searchsorted(t,row[key]-1e-9))
            if k>=len(t) or abs(t[k]-row[key])>1e-8:raise ValueError('Trigger must lie on saved grid')
            row.update({label+'_'+name:value for name,value in activation_state(
                arrays['gap_error_m'][k],arrays['gap_rate_m_s'][k],row['satellite_id'],config).items()})
    return rows


def off_before_dwell(arrays,config,events):
    """First sampled off condition after HIGH start, including delayed returns.

    An early off sample does not imply it stays true at the return gate; the
    existing return_gates diagnostic checks eligibility at that separate time.
    """
    t=arrays['time_s'];off=config.controller.threshold_off_m_s2;rows=[]
    for pair in cycle_pairs(events,float(t[-1])):
        rise=pair['rise'];gate=rise['complete_time']+config.controller.dwell_time_s
        observed_end=min(pair['return_event']['command_time'],float(t[-1])) if pair['return_event'] else float(t[-1])
        k0=int(np.searchsorted(t,rise['start_time']-1e-9));k1=int(np.searchsorted(t,observed_end+1e-9))
        indices=np.flatnonzero(arrays['q_m_s2'][k0:k1,pair['satellite_id']-1]<=off)
        first=float(t[k0+indices[0]]) if indices.size else None
        rows.append(dict(satellite_id=pair['satellite_id'],local_cycle_number=pair['local_cycle_number'],
            first_off_after_HIGH_start_s=first,earliest_return_gate_s=gate,
            gate_observed=bool(gate<=t[-1]+1e-9),
            off_observed_before_gate=bool(first is not None and first<gate-1e-9),
            first_off_to_LOW_command_s=pair['return_event']['command_time']-first
                if first is not None and pair['return_event'] else None))
    return rows


def fixed_state_forecasts(error,rate,time_s,allowed,config,factors=(.3,1.,3.)):
    """Change Kp in the SAME baseline state: P offset and slope both scale."""
    result=[];kp0=1/21600.**2;kd=config.derived.controller_kd_s_inv
    for factor in factors:
        kp=kp0*factor;p,d=decompose_q(np.asarray(error)[None,:],np.asarray(rate)[None,:],kp,kd)
        candidates=[]
        for j in range(len(p[0])):
            pred=predict_crossing(p[0,j],d[0,j],time_s,allowed[j],kp,kd,
                config.controller.threshold_on_m_s2,config.controller.control_period_s)
            if pred['predicted_command_s'] is not None:
                candidates.append(dict(kp_factor=factor,satellite_id=j+1,q_P_m_s2=float(p[0,j]),
                    q_D_m_s2=float(d[0,j]),on_margin_m_s2=float(config.controller.threshold_on_m_s2-p[0,j]-d[0,j]),**pred))
        result.append(min(candidates,key=lambda x:(x['predicted_command_s'],x['satellite_id']))
            if candidates else dict(kp_factor=factor,satellite_id=None,predicted_command_s=None))
    return result
