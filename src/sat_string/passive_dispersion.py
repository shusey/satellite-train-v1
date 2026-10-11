"""Independent ALL-LOW profiles using the existing orbital_rates and RK4.

The scalar reference spacecraft stays unchanged. A positive multiplier denotes
B_j/B0, where B=Cd*A/m (not the reciprocal convention sometimes used for B).
No controller/attitude state machine or inter-profile interaction is present.
"""
from pathlib import Path
from time import perf_counter
import numpy as np
from .controller import initial_formation, gap_kinematics
from .dynamics import orbital_rates, density_at_state
from .integrator import rk4_step
from .io import prepare_output_directory, write_json


def centered_normal(count, sigma, seed):
    """Draw N(0,sigma²), then remove sample mean; do NOT rescale sample std."""
    if count<2 or sigma<0 or not np.isfinite(sigma):
        raise ValueError('Finite nonnegative sigma and at least two satellites required')
    values=np.random.default_rng(seed).normal(0.,sigma,count)
    return values-values.mean()


def centered_alternating(count, sigma):
    values=(-1.)**np.arange(count)*sigma
    return values-values.mean()


def validate_profiles(offsets_m, multipliers, count):
    offsets=np.asarray(offsets_m,float); factors=np.asarray(multipliers,float)
    if offsets.ndim==1:offsets=offsets[None,:]
    if factors.ndim==1:factors=factors[None,:]
    if offsets.shape!=factors.shape or offsets.ndim!=2 or offsets.shape[1]!=count:
        raise ValueError('Matching profiles shaped (cases,N) required')
    if not np.all(np.isfinite(offsets+factors)) or np.any(factors<=0):
        raise ValueError('Finite offsets and positive B multipliers required')
    if np.max(np.abs(offsets.mean(axis=1)))>1e-10 or np.max(np.abs(factors.mean(axis=1)-1))>1e-12:
        raise ValueError('Profile means must preserve a0 and B0')
    return offsets,factors


def passive_rhs(t,state,config,provider,factors):
    """Reuse original physical ODE; actual NRLMSIS positions for every member."""
    a,u=state
    rho=density_at_state(t,a,u,config,disturbance_enabled=False,density_provider=provider)
    da,du=orbital_rates(a,rho,np.full(a.shape,config.spacecraft.area_low_m2),
        mu_m3_s2=config.constants.mu_m3_s2,
        drag_coefficient=config.spacecraft.drag_coefficient*factors,
        mass_kg=config.spacecraft.mass_kg)
    if np.any(a<=config.constants.earth_radius_m):raise FloatingPointError('Orbit crossed Earth')
    return np.vstack((da,du))


def initial_profile_state(config,offsets_m,multipliers):
    offsets,factors=validate_profiles(offsets_m,multipliers,config.formation.satellite_count)
    a,u=initial_formation(config.formation.satellite_count,config.derived.reference_radius_m,
                          config.derived.target_spacing_rad)
    return np.vstack(((a+offsets).ravel(),np.broadcast_to(u,offsets.shape).ravel())),factors.ravel()


def theory_initial_gap_rate(config,offsets_m):
    """Fixed a_ref convention: de/dt≈-(3/2)n_ref*(delta_a[j+1]-delta_a[j])."""
    return -1.5*config.derived.reference_mean_motion_rad_s*np.diff(offsets_m,axis=-1)


def integrate_profiles(config,profiles,provider,output,*,start_s=0.,initial_state=None,progress=None):
    """Batch independent profiles, preserving original2s RK4/10s output.

    One shared vectorized density call has no coupling between profiles. New
    trajectories are queried separately: this is NOT a frozen-density replay.
    Continuations begin from saved physical state, never repeat their prefix.
    Large full states remain local .npy memory maps; compact metrics are separate.
    """
    offsets=np.array([p['delta_a_m'] for p in profiles]);factors=np.array([p['B_multiplier'] for p in profiles])
    state,flat_factors=initial_profile_state(config,offsets,factors)
    if initial_state is not None:
        supplied=np.asarray(initial_state,float)
        if supplied.shape!=state.shape or not np.all(np.isfinite(supplied)):raise ValueError('Bad continuation state')
        state=supplied.copy()
    end=config.simulation.duration_s;step=config.simulation.integration_step_s;save=config.simulation.output_step_s
    if end<=start_s or save/step!=round(save/step) or start_s/save!=round(start_s/save):
        raise ValueError('Exact step-aligned output and continuation required')
    out=prepare_output_directory(output);times=np.arange(start_s,end+save/2,save)
    if times[-1]!=end:raise ValueError('End must be output-aligned')
    np.save(out/'time_s.npy',times)
    shape=(len(times),len(profiles),config.formation.satellite_count)
    arrays={key:np.lib.format.open_memmap(out/f'{key}.npy',mode='w+',dtype='float64',shape=shape)
            for key in ('a_m','lambda_rad','density_kg_m3')}
    write_json(out/'metadata.json',dict(configuration=config.to_dict(),profiles=profiles,
        start_s=start_s,end_s=end,all_low=True,control_enabled=False,
        density='NRLMSIS2 evaluated at each changed trajectory; same historical drivers',
        integrator='Original orbital_rates + rk4_step; independent vectorized profiles',complete=False))
    t=float(start_s);index=0;began=perf_counter();next_progress=start_s+86400
    rhs=lambda tt,ss:passive_rhs(tt,ss,config,provider,flat_factors)
    while True:
        if index<len(times) and t==times[index]:
            arrays['a_m'][index]=state[0].reshape(shape[1:]);arrays['lambda_rad'][index]=state[1].reshape(shape[1:])
            arrays['density_kg_m3'][index]=provider(t,state[0],state[1]).reshape(shape[1:]);index+=1
        if t>=end:break
        state=rk4_step(rhs,t,state,min(step,end-t));t+=min(step,end-t)
        if t>=next_progress:
            for array in arrays.values():array.flush()
            if progress:progress(t,end,perf_counter()-began)
            next_progress+=86400
    for array in arrays.values():array.flush()
    write_json(out/'completion.json',dict(complete=True,profiles=len(profiles),start_s=start_s,end_s=end,
        elapsed_s=perf_counter()-began,save_count=index,all_low=True,transitions=0))
    return state


def state_metrics(t,a,u,config):
    """Unweighted all-time/all-gap RMS, existing fixed-reference definition."""
    e,v,_=gap_kinematics(a,u,reference_radius_m=config.derived.reference_radius_m,
        target_spacing_rad=config.derived.target_spacing_rad,mu_m3_s2=config.constants.mu_m3_s2)
    maximum=np.max(np.abs(e),axis=-1)
    labelled_spacing=config.formation.target_spacing_m+e
    reversal=np.flatnonzero(np.any(labelled_spacing<=0,axis=-1))
    result=dict(max_error_m=float(maximum.max()),rms_error_m=float(np.sqrt(np.mean(e**2))),
                final_max_error_m=float(maximum[-1]),max_gap_rate_m_s=float(np.abs(v).max()),
                max_error_percent_of_100km=float(maximum.max()/config.formation.target_spacing_m*100),
                first_exceedance_s={},gap_rms_m=np.sqrt(np.mean(e**2,axis=0)).tolist(),
                min_labelled_spacing_m=float(labelled_spacing.min()),
                first_label_order_reversal_s=float(t[reversal[0]]) if reversal.size else None)
    for threshold in (100,500,1000):
        ix=np.flatnonzero(maximum>threshold)
        result['first_exceedance_s'][str(threshold)]=float(t[ix[0]]) if ix.size else None
    return result
