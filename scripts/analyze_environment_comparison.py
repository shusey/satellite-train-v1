"""Offline comparison of weather-preregistered ALL-LOW runs, no physics rerun."""
import csv
import json
from pathlib import Path
import sys
import numpy as np
import yaml
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from sat_string.config import load_config
from sat_string.environment_integrals import all_low_integrals,gap_rate_from_exposure
from sat_string.nrlmsis import TrajectoryDensityCache
from sat_string.nrlmsis_analysis import current_radius_arc_error,density_statistics
from sat_string.necessity import exceedance_metrics
from sat_string.io import write_json,prepare_output_directory
from run_nrlmsis_case import make_environment

ROOT=Path(__file__).resolve().parents[1]
OUT=ROOT/'outputs/real_environment_20261011'
DOC=ROOT/'docs/real_environment_applicability'


def csv_save(path,rows):
    with path.open('w',encoding='utf-8',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def audit_and_quadrature(run,group,t,a,u,rho,cfg):
    env=ROOT/f'configs/real_environment_applicability/{group}.yaml'
    _,backend,_=make_environment(env,OUT/f'weather/{group}/drivers.csv',0)
    jumps=np.flatnonzero((t>0)&(t%10800==0))
    ix=np.unique(np.r_[np.linspace(0,len(t)-1,250,dtype=int),jumps,jumps-1])
    direct=backend(t[ix,None],a[ix],u[ix]);err=np.abs(rho[ix]/direct-1)
    interior=[]
    for anchor in np.unique(np.floor(t[ix[:-1]]/120)*120):
        if anchor+5>t[-1]:continue
        j=int(anchor/10);aa=a[j];uu=u[j];aq=.5*(a[j]+a[j+1]);uq=.5*(u[j]+u[j+1])
        cache=TrajectoryDensityCache(backend,120,cfg.constants.mu_m3_s2)
        cache(anchor,aa,uu)
        interior.extend(np.abs(cache(anchor+5,aq,uq)/backend(anchor+5,aq,uq)-1))
    interior=np.asarray(interior)
    audit=dict(max_saved_relative_error=float(err.max()),max_interior_relative_error=float(interior.max()),
               saved_density_points=int(err.size),interior_density_points=int(interior.size),passed=bool(max(err.max(),interior.max())<.001),
               correct_epoch=str(backend.geometry.epoch_utc),criterion='relative density error<0.001; actual run coordinates')
    if not audit['passed']:raise ValueError('New period failed direct/cache density audit')
    # 3h drivers and midnight DOY are discontinuous: left side must use old inputs.
    before=backend(t[jumps,None]-.001,a[jumps],u[jumps]-.001*np.sqrt(cfg.constants.mu_m3_s2/a[jumps]**3))
    dt=np.diff(t)[:,None]
    trapezoid=np.sum(.5*(rho[:-1]+rho[1:])*dt,axis=0)
    correction=np.sum(.5*(before-rho[jumps])*dt[jumps-1],axis=0)
    corrected=trapezoid+correction
    # Original 2s RK4 includes the new input in the last stage of the old interval.
    # Estimate its step-discontinuity impulse bias, without modifying that solver.
    impulse=(rho[jumps]-before)*cfg.simulation.integration_step_s/6
    bias=impulse.sum(axis=0)
    gaps=np.diff(impulse,axis=1)
    k=1.5*cfg.spacecraft.drag_coefficient*cfg.spacecraft.area_low_m2*cfg.constants.mu_m3_s2/(cfg.spacecraft.mass_kg*cfg.derived.reference_radius_m)
    conservative_spacing_bias=k*np.max(np.sum(np.abs(gaps)*(t[-1]-t[jumps,None]),axis=0))
    return audit,dict(corrected_density_quadrature_kg_s_m3=corrected.tolist(),
        endpoint_discontinuity_bias_estimate_kg_s_m3=bias.tolist(),
        spacing_bias_absolute_sum_estimate_m=float(conservative_spacing_bias),
        caution='First-order diagnostic bound-like estimate, not a rigorous solver error bound; saved10s smooth interpolation and RK rounding also contribute'),corrected,bias


def analyze(group,run):
    cfg=load_config(run/'normalized_config.yaml')
    if (np.load(run/'results.npz')['area_m2']!=cfg.spacecraft.area_low_m2).any():raise ValueError('Only ALL-LOW allowed')
    c,s=cfg.constants,cfg.spacecraft
    with np.load(run/'results.npz') as d:
        t,a,u,e,v,rho=[d[k] for k in ['time_s','a_m','lambda_rad','gap_error_m','gap_rate_m_s','density_kg_m3']]
    jr,jd=all_low_integrals(a,a[0],c.mu_m3_s2,s.drag_coefficient,s.area_low_m2,s.mass_kg)
    djr,djd=np.diff(jr,axis=1),np.diff(jd,axis=1)
    ev=np.max(np.abs(e),axis=1);rms=np.sqrt(np.mean(e**2,axis=1))
    arc=current_radius_arc_error(a,u,cfg.formation.target_spacing_m);arcm=np.max(np.abs(arc),axis=1)
    audit,q,quadrature,bias=audit_and_quadrature(run,group,t,a,u,rho,cfg)
    q['max_adjacent_quadrature_plus_RK_bias_vs_inferred_kg_s_m3']=float(np.max(np.abs(np.diff(quadrature+bias)-djr[-1])))
    predicted=gap_rate_from_exposure(jr[-1],a[0],cfg.derived.reference_radius_m,c.mu_m3_s2,s.drag_coefficient,s.area_low_m2,s.mass_kg)
    gaprows=[dict(group=group,gap_id=i+1,max_error_m=float(np.max(np.abs(e[:,i]))),rms_error_m=float(np.sqrt(np.mean(e[:,i]**2))),
      final_error_m=float(e[-1,i]),final_rate_m_s=float(v[-1,i]),final_delta_Jrho=float(djr[-1,i]),final_delta_Jdrag=float(djd[-1,i])) for i in range(e.shape[1])]
    metrics=dict(group=group,run=str(run.relative_to(ROOT)),epoch=yaml.safe_load((ROOT/f'configs/real_environment_applicability/{group}.yaml').read_text())['epoch_utc'],
      max_error_m=float(ev.max()),rms_error_m=float(np.sqrt(np.mean(e**2))),final_max_error_m=float(ev[-1]),
      max_error_percent_of_100km=float(ev.max()/cfg.formation.target_spacing_m*100),
      max_gap_speed_m_s=float(np.max(np.abs(v))),final_max_gap_speed_m_s=float(np.max(np.abs(v[-1]))),
      max_cumulative_adjacent_Jrho=float(np.max(np.abs(djr))),final_max_adjacent_Jrho=float(np.max(np.abs(djr[-1]))),
      max_cumulative_adjacent_Jdrag_m_s=float(np.max(np.abs(djd))),final_max_adjacent_Jdrag_m_s=float(np.max(np.abs(djd[-1]))),
      final_Jrho=jr[-1].tolist(),final_Jdrag=jd[-1].tolist(),mean_altitude_loss_m=float(np.mean(a[0]-a[-1])),
      rate_identity_max_abs_error_m_s=float(np.max(np.abs(predicted-v[-1]))),
      sampled_velocity_integral_final_error_m=float(np.max(np.abs(e[-1]-e[0]-np.sum(.5*(v[:-1]+v[1:])*np.diff(t)[:,None],axis=0)))),
      first_order_weighted_density_final_spacing_m=(1.5*s.drag_coefficient*s.area_low_m2*c.mu_m3_s2/(s.mass_kg*cfg.derived.reference_radius_m)*np.sum(.5*(djr[:-1]+djr[1:])*np.diff(t)[:,None],axis=0)).tolist(),
      density=density_statistics(rho),audit=audit,quadrature_check=q,
      virtual_tolerances=[exceedance_metrics(t,ev,k) for k in (100,500,1000)],
      current_radius_arc=dict(max_error_m=float(arcm.max()),rms_error_m=float(np.sqrt(np.mean(arc**2))),final_max_error_m=float(arcm[-1]),
                             virtual_tolerances=[exceedance_metrics(t,arcm,k) for k in (100,500,1000)]),
      tail=dict(last7days_max_error_m=float(ev[t>=t[-1]-7*86400].max()),end_max_error_m=float(ev[-1]),
                end14day_error_growth_m_s=float(np.polyfit(t[t>=t[-1]-14*86400]-t[-1],ev[t>=t[-1]-14*86400],1)[0])),
      integral_definition='Orbit-inferred effective integrals from exact ALL-LOW averaged ODE identities; independently checked with10s density quadrature and jump correction',
      forcing_status='Continuous NRLMSIS forcing: final gap rate is not a post-disturbance frozen residual')
    hourly=np.arange(0,len(t),360)
    csv_save(DOC/'data'/f'{group}_gaps.csv',gaprows)
    daily=np.arange(0,len(t),8640)
    csv_save(DOC/'data'/f'{group}_daily.csv',[dict(time_day=float(t[i]/86400),max_error_m=float(ev[i]),rms_error_m=float(rms[i]),
        max_gap_speed_m_s=float(np.max(np.abs(v[i]))),max_adjacent_Jrho=float(np.max(np.abs(djr[i]))),
        max_adjacent_Jdrag_m_s=float(np.max(np.abs(djd[i]))),mean_density=float(rho[i].mean())) for i in daily])
    np.savez_compressed(OUT/f'{group}_hourly.npz',time_s=t[hourly],max_error_m=ev[hourly],rms_error_m=rms[hourly],gap_error_m=e[hourly],
      gap_rate_m_s=v[hourly],adjacent_Jrho=djr[hourly],adjacent_Jdrag=djd[hourly],mean_density=rho[hourly].mean(axis=1),current_arc_error_m=arc[hourly])
    # Mechanism window is set by space weather, not by the largest spacing error.
    weather=np.genfromtxt(OUT/f'weather/{group}/drivers.csv',delimiter=',',names=True,dtype=None,encoding='utf-8')
    peak=np.argmax(weather['ap1']);center=float((np.datetime64(weather['utc'][peak])-np.datetime64(metrics['epoch']))/np.timedelta64(1,'s'))
    use=(t>=max(0,center-10800))&(t<=min(t[-1],center+10800))
    np.savez_compressed(OUT/f'{group}_mechanism.npz',time_s=t[use],adjacent_density=np.diff(rho[use],axis=1),
                        adjacent_Jrho=djr[use],adjacent_Jdrag=djd[use],gap_rate_m_s=v[use],gap_error_m=e[use])
    write_json(DOC/'data'/f'{group}_metrics.json',metrics)
    print(group,{k:metrics[k] for k in ('max_error_m','rms_error_m','final_max_error_m','max_cumulative_adjacent_Jrho','final_max_adjacent_Jrho')},flush=True)
    return metrics,cfg


def main():
    (DOC/'data').mkdir(parents=True,exist_ok=True)
    if (DOC/'data'/'comparison.json').exists():raise FileExistsError('Completed comparison exists; do not overwrite')
    entries={'quiet':ROOT/'outputs/nrlmsis2_20261010/runs/all_low_60d',
             **{g:OUT/f'runs/{g}_60d' for g in ('active','strong')}}
    metrics=[];reference=None
    for group,run in entries.items():
        cached=DOC/'data'/f'{group}_metrics.json'
        if cached.exists():
            m=json.loads(cached.read_text(encoding='utf-8'));cfg=load_config(run/'normalized_config.yaml')
            if m['run']!=str(run.relative_to(ROOT)) or not m['audit']['passed']:raise ValueError('Wrong/unvalidated cached analysis')
        else:m,cfg=analyze(group,run)
        metrics.append(m)
        if reference is None:reference=cfg.to_dict()
        elif cfg.to_dict()!=reference:raise ValueError('Physical/numerical config mismatch')
        with np.load(run/'results.npz') as d:
            expected=np.load(entries['quiet']/'results.npz')
            np.testing.assert_array_equal(d['a_m'][0],expected['a_m'][0]);np.testing.assert_array_equal(d['lambda_rad'][0],expected['lambda_rad'][0])
            np.testing.assert_array_equal(d['time_s'],expected['time_s'])
        if (run/'switch_events.csv').read_text().count('\n')!=1:raise ValueError('Unexpected transitions')
        actual=json.loads((run/'density_environment.json').read_text(encoding='utf-8'))
        intended=yaml.safe_load((ROOT/f'configs/real_environment_applicability/{group}.yaml').read_text())
        for k in ('density_model','msis_version','epoch_utc','inclination_deg','raan_deg','initial_phase_deg','geomagnetic_activity'):
            if actual['environment'][k]!=intended[k]:raise ValueError(f'Wrong environment input: {k}')
        import hashlib
        if actual['weather_sha256']!=hashlib.sha256((OUT/f'weather/{group}/drivers.csv').read_bytes()).hexdigest():raise ValueError('Wrong frozen weather input')
    write_json(DOC/'data'/'comparison.json',dict(results=metrics,input_validation='identical normalized physical/numerical config, initial state/output times; different epoch/weather only',
        control_case_calculated=False,pre_registered_protocol=json.loads((OUT/'protocol.json').read_text())))


if __name__=='__main__': main()
