"""Weather-fixed ALL-LOW dispersion screening; reuse baseline, batch new profiles."""
import argparse,csv,hashlib,json,os,sys
from pathlib import Path
import numpy as np
import yaml
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT/'src'))
from sat_string.config import load_config,config_from_dict,dump_normalized_config
from sat_string.io import write_json
from sat_string.passive_dispersion import (centered_normal,centered_alternating,
    integrate_profiles,state_metrics,theory_initial_gap_rate,initial_profile_state)
from sat_string.controller import gap_kinematics
from run_nrlmsis_case import make_environment

OUT=Path(os.environ.get('SAT_DISPERSION_OUTPUT',ROOT/'outputs/initial_dispersion_20261011'))
DOC=Path(os.environ.get('SAT_DISPERSION_DOC',ROOT/'docs/initial_dispersion_screening'));DATA=DOC/'data'
SPEC=ROOT/'configs/initial_dispersion_screening.yaml'

def write_csv(path,rows):
    with path.open('w',encoding='utf-8',newline='') as f:
        fields=list(dict.fromkeys(key for row in rows for key in row))
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader();w.writerows(rows)

def configuration(days):
    spec=yaml.safe_load(SPEC.read_text());base=load_config(ROOT/spec['baseline_run']/'normalized_config.yaml')
    raw=base.to_dict();raw['simulation']['duration_s']=days*86400
    cfg=config_from_dict(raw)
    checked=cfg.to_dict();checked['simulation']['duration_s']=base.simulation.duration_s
    if checked!=base.to_dict() or cfg.disturbance.enabled:raise ValueError('Unexpected baseline changes')
    if spec['control_enabled'] or spec['area_mode']!='all_low':raise ValueError('Screening must be ALL LOW')
    if cfg.formation.satellite_count!=21 or spec['integration_step_s']!=cfg.simulation.integration_step_s or spec['output_step_s']!=cfg.simulation.output_step_s:
        raise ValueError('Screening specification must agree with frozen baseline')
    return spec,cfg

def profile(factor,scale,seed,kind='random'):
    n=21;z_a=centered_normal(n,1,seed) if seed is not None else centered_alternating(n,1)
    z_b=centered_normal(n,1,seed+10000) if seed is not None else centered_alternating(n,1)
    off=z_a*scale if factor=='a' else np.zeros(n);eps=z_b*scale if factor=='B' else np.zeros(n)
    return dict(id=f'{factor}_{scale:g}_{kind}_{seed if seed is not None else "alt"}',factor=factor,
      nominal_sigma=scale,seed=seed,B_seed=seed+10000 if seed is not None else None,kind=kind,
      delta_a_m=off.tolist(),B_multiplier=(1+eps).tolist(),epsilon=eps.tolist(),
      realized_sigma_a_m=float(off.std()),realized_sigma_B_fraction=float(eps.std()))

def initial_profiles(spec):
    result=[profile(f,scale,seed) for f,key in [('a','sigma_a_m'),('B','sigma_B_fraction')]
            for scale in spec[key] if scale>0 for seed in spec['seeds']]
    return result+[profile('a',spec['stress_sigma_a_m'],None,'stress'),
                   profile('B',spec['stress_sigma_B_fraction'],None,'stress')]

def baseline_metrics(spec,cfg):
    cached=DATA/'baseline_prefixes.json'
    if cached.exists():return json.loads(cached.read_text())['results']
    with np.load(ROOT/spec['baseline_run']/'results.npz') as d:
        t=d['time_s'];a=d['a_m'];u=d['lambda_rad']
        rows=[dict(case_id='zero_reused',factor='zero',nominal_sigma=0,seed=None,kind='baseline',days=day,
             **state_metrics(t[t<=day*86400],a[t<=day*86400],u[t<=day*86400],cfg)) for day in (1,7,30,60)]
    write_json(cached,dict(reused=True,source=spec['baseline_run'],results=rows))
    return rows

def metrics_for_stage(path,cfg):
    if not (path/'completion.json').exists():raise ValueError('Unfinished stage')
    meta=json.loads((path/'metadata.json').read_text());t=np.load(path/'time_s.npy')
    a=np.load(path/'a_m.npy',mmap_mode='r');u=np.load(path/'lambda_rad.npy',mmap_mode='r');rows=[]
    if t[0]==0:
        for k,p in enumerate(meta['profiles']):
            for day in (1,7):
                if day*86400>t[-1]:continue
                use=t<=day*86400
                rows.append(dict(case_id=p['id'],factor=p['factor'],nominal_sigma=p['nominal_sigma'],seed=p['seed'],kind=p['kind'],
                    days=day,**state_metrics(t[use],a[use,k],u[use,k],cfg)))
    else:
        original=OUT/'initial_7d';m0=json.loads((original/'metadata.json').read_text());t0=np.load(original/'time_s.npy')
        a0=np.load(original/'a_m.npy',mmap_mode='r');u0=np.load(original/'lambda_rad.npy',mmap_mode='r')
        for k,p in enumerate(meta['profiles']):
            index=next(i for i,q in enumerate(m0['profiles']) if q['id']==p['id'])
            tt=np.r_[t0,t[1:]];aa=np.concatenate((a0[:,index],a[1:,k]));uu=np.concatenate((u0[:,index],u[1:,k]))
            rows.append(dict(case_id=p['id'],factor=p['factor'],nominal_sigma=p['nominal_sigma'],seed=p['seed'],kind=p['kind'],
                days=t[-1]/86400,**state_metrics(tt,aa,uu,cfg)))
    return rows

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--stage',choices=['initial','extend','combination','analyze'],required=True)
    args=parser.parse_args();spec,cfg=configuration(7);DATA.mkdir(parents=True,exist_ok=True)
    if args.stage=='analyze':
        rows=baseline_metrics(spec,cfg)
        for name in ('initial_7d','extended_30d','combination_7d'):
            if (OUT/name/'completion.json').exists():rows+=metrics_for_stage(OUT/name,cfg)
        write_json(DATA/'metrics.json',dict(results=rows,main_spacing='fixed-reference',RMS='sqrt(mean(e_i(t)^2)) over every stored time/gap',all_low=True))
        flattened=[{k:v for k,v in r.items() if k not in ('first_exceedance_s','gap_rms_m')}|{f'first_{key}_s':value for key,value in r['first_exceedance_s'].items()} for r in rows]
        write_csv(DATA/'metrics.csv',flattened);print([(r['case_id'],r['days'],round(r['max_error_m'],3)) for r in rows],flush=True);return
    start=0.;state=None
    if args.stage=='initial':
        profiles=initial_profiles(spec)
        if (DATA/'profiles.json').exists():raise FileExistsError('Initial profiles are frozen')
        write_json(DATA/'profiles.json',dict(spec=spec,profiles=profiles,
          limitations='User-specified provisional ranges, not measured populations; matched standardized patterns across scale levels',
          extension_rule=spec['extension_rule'],baseline_reused=True,control_reference_not_planned=True))
        write_csv(DATA/'satellite_profiles.csv',[dict(case_id=p['id'],satellite=j+1,delta_a_m=p['delta_a_m'][j],
          B_multiplier=p['B_multiplier'][j],epsilon=p['epsilon'][j]) for p in profiles for j in range(21)])
        output=OUT/'initial_7d'
    elif args.stage=='extend':
        original=OUT/'initial_7d';metrics=metrics_for_stage(original,cfg)
        frozen=json.loads((DATA/'profiles.json').read_text());selected=[];decisions=[]
        for f in ('a','B'):
            levels=sorted({p['nominal_sigma'] for p in frozen['profiles'] if p['factor']==f and p['kind']=='random'})
            for level in levels:
                candidates=[r for r in metrics if r['factor']==f and r['nominal_sigma']==level and r['kind']=='random' and r['days']==7]
                ordered=sorted(candidates,key=lambda r:r['max_error_m']);median=ordered[len(ordered)//2]
                chosen=[r for r in candidates if 100<=r['max_error_m']<=500]
                if not chosen and median['max_error_m']<100:chosen=[median]
                # The smallest B dispersion gets one long example, even if already
                # large, to compare growth with the smallest a dispersion at30d.
                if f=='B' and level==levels[0] and not chosen:chosen=[median]
                selected.extend(r['case_id'] for r in chosen)
                decisions.append(dict(factor=f,sigma=level,median_seed=median['seed'],range_7d_m=[ordered[0]['max_error_m'],ordered[-1]['max_error_m']],
                     selected_case_ids=[r['case_id'] for r in chosen],reason='100-500m: extend; below100: representative; smallest B: one growth example; other large levels stop'))
        write_json(DATA/'extension_decision.json',dict(decisions=decisions,additional_seeds=False,sigma_a_50_m=False,
             reason_no50='10m already large; no need to explore an upper bound'))
        profiles=[p for p in frozen['profiles'] if p['id'] in selected]
        if not profiles:print('No extensions needed');return
        meta=json.loads((original/'metadata.json').read_text());a=np.load(original/'a_m.npy',mmap_mode='r');u=np.load(original/'lambda_rad.npy',mmap_mode='r')
        ix=[next(i for i,q in enumerate(meta['profiles']) if q['id']==p['id']) for p in profiles]
        state=np.vstack((a[-1,ix].ravel(),u[-1,ix].ravel()));start=7*86400.;spec,cfg=configuration(30);output=OUT/'extended_30d'
    else:
        results=metrics_for_stage(OUT/'initial_7d',cfg)
        if not all(any(r['factor']==f and r['days']==7 and r['max_error_m']>=500 for r in results) for f in ('a','B')):
            print('Both factors not important: skip combination');return
        p=profile('a',1.,0,'combination');eps=centered_normal(21,.01,10000);p.update(id='combined_a1m_B1percent_seed0',factor='combined',B_multiplier=(1+eps).tolist(),epsilon=eps.tolist(),realized_sigma_B_fraction=float(eps.std()))
        profiles=[p];write_json(DATA/'combination_profile.json',dict(profiles=profiles,
            purpose='One independent A/B pattern combination; not a parameter sweep',sigma_a_m=1,sigma_B_fraction=.01));output=OUT/'combination_7d'
    env,backend,provider=make_environment(ROOT/spec['epoch_environment'],ROOT/spec['weather'],spec['cache_seconds'])
    callback=lambda t,end,wall:print(f'{args.stage}: {t/86400:g}/{end/86400:g} days; {len(profiles)} profiles; {wall:.1f}s',flush=True)
    integrate_profiles(cfg,profiles,provider,output,start_s=start,initial_state=state,progress=callback)
    dump_normalized_config(cfg,output/'normalized_config.yaml')
    write_json(output/'environment.json',dict(spec=env,cache_seconds=spec['cache_seconds'],
        weather_sha256=hashlib.sha256((ROOT/spec['weather']).read_bytes()).hexdigest(),
        msis_version='2.0',pymsis_version='0.12.0',total_mass_density_units='kg/m3',backend_points=backend.points_evaluated,
        maximum_cache_radius_offset_m=provider.maximum_radius_offset_m,maximum_cache_phase_offset_rad=provider.maximum_phase_offset_rad))
    print('Saved',output,flush=True)

if __name__=='__main__':main()
