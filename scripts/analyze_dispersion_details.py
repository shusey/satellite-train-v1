"""Compact physical diagnostics from saved dispersion runs, no new simulation."""
import json,sys
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from sat_string.controller import gap_kinematics
from sat_string.passive_dispersion import theory_initial_gap_rate
from run_dispersion_screening import ROOT,OUT,DATA,configuration,write_csv
from run_nrlmsis_case import make_environment

def group_summary(rows):
    groups=[]
    for factor in ('a','B'):
        for sigma in sorted({r['nominal_sigma'] for r in rows if r['factor']==factor and r['kind']=='random'}):
            results=[r for r in rows if r['factor']==factor and r['nominal_sigma']==sigma and r['kind']=='random']
            horizons={}
            for days in (1,7,30):
                sample=[r for r in results if r['days']==days]
                if not sample:continue
                vals=np.array([r['max_error_m'] for r in sample]);rms=[r['rms_error_m'] for r in sample]
                ordered=sorted(sample,key=lambda r:r['max_error_m']);representative=ordered[len(ordered)//2]
                horizons[str(days)]=dict(n=len(sample),median_max_error_m=float(np.median(vals)),
                    min_max_error_m=float(vals.min()),max_max_error_m=float(vals.max()),
                    median_rms_m=float(np.median(rms)),representative_seed=representative['seed'],
                    representative_case_id=representative['case_id'])
            largest=max(h['median_max_error_m'] for h in horizons.values())
            rank='A' if largest>=1000 else ('B' if largest>=300 else ('C' if largest>=70 else 'D'))
            groups.append(dict(factor=factor,sigma=sigma,horizons=horizons,provisional_rank=rank,
              rank_caution='Dynamical significance within user-specified ranges, NOT proof of real occurrence'))
    return groups

def main():
    if not (OUT/'initial_7d'/'completion.json').exists():raise ValueError('Initial stage incomplete')
    spec,cfg=configuration(7);directory=OUT/'initial_7d'
    metadata=json.loads((directory/'metadata.json').read_text());profiles=metadata['profiles']
    t=np.load(directory/'time_s.npy');a=np.load(directory/'a_m.npy',mmap_mode='r');u=np.load(directory/'lambda_rad.npy',mmap_mode='r')
    rho=np.load(directory/'density_kg_m3.npy',mmap_mode='r');rows=json.loads((DATA/'metrics.json').read_text())['results']
    groups=group_summary(rows);write_json=lambda p,d:p.write_text(json.dumps(d,indent=2,ensure_ascii=False,allow_nan=False),encoding='utf-8')
    write_json(DATA/'group_summary.json',dict(groups=groups,seeds=5,stress_excluded_from_random_statistics=True,
       matched_scales='same standardized pattern within each factor',sample_sigma_not_rescaled=True))
    theory=[]
    for k,p in enumerate(profiles):
        expected=theory_initial_gap_rate(cfg,p['delta_a_m'])
        _,actual,_=gap_kinematics(a[0,k],u[0,k],reference_radius_m=cfg.derived.reference_radius_m,
            target_spacing_rad=cfg.derived.target_spacing_rad,mu_m3_s2=cfg.constants.mu_m3_s2)
        if p['factor']=='a':
            for j in range(20):theory.append(dict(case_id=p['id'],gap=j+1,
                delta_a_initial_m=float(np.diff(a[0,k])[j]),theory_rate_m_s=float(expected[j]),actual_rate_m_s=float(actual[j]),
                error_m_s=float(actual[j]-expected[j]),relative_error=abs(float((actual[j]-expected[j])/expected[j])) if abs(expected[j])>1e-12 else None))
    write_csv(DATA/'initial_drift_theory.csv',theory)
    if max(abs(r['error_m_s']) for r in theory)>1e-6:raise ValueError('Initial theory inconsistent with physical model')
    representative=[]
    with np.load(ROOT/spec['baseline_run']/'results.npz') as base:
        base_t=base['time_s'];base_a=base['a_m'];base_u=base['lambda_rad'];use=base_t<=7*86400
        base_a=base_a[use];base_u=base_u[use]
        np.testing.assert_array_equal(t,base_t[use])
    be,bv,_=gap_kinematics(base_a,base_u,reference_radius_m=cfg.derived.reference_radius_m,
        target_spacing_rad=cfg.derived.target_spacing_rad,mu_m3_s2=cfg.constants.mu_m3_s2)
    for f,sigma in [('a',.1),('B',.003)]:
        group=next(g for g in groups if g['factor']==f and g['sigma']==sigma)
        ident=group['horizons']['7']['representative_case_id'];k=next(i for i,p in enumerate(profiles) if p['id']==ident);p=profiles[k]
        e,v,_=gap_kinematics(a[:,k],u[:,k],reference_radius_m=cfg.derived.reference_radius_m,
            target_spacing_rad=cfg.derived.target_spacing_rad,mu_m3_s2=cfg.constants.mu_m3_s2)
        da_excess=np.diff(a[:,k]-base_a,axis=1)
        B0=cfg.spacecraft.drag_coefficient*cfg.spacecraft.area_low_m2/cfg.spacecraft.mass_kg
        B=B0*np.array(p['B_multiplier']);db=np.diff(B)
        # da_dot=-B*rho*sqrt(mu*a): actual vs uniform-B reference contribution.
        adot=-B[None,:]*rho[:,k]*np.sqrt(cfg.constants.mu_m3_s2*a[:,k])
        badot=-B0*np.load(directory/'density_kg_m3.npy',mmap_mode='r')[:,k]*np.sqrt(cfg.constants.mu_m3_s2*a[:,k])
        initial_B_only_adot_excess=adot[0]-badot[0]
        summary=dict(factor=f,case_id=ident,seed=p['seed'],nominal_sigma=sigma,
            delta_B_gap_m2_kg=db.tolist(),initial_B_only_adot_excess_m_s=initial_B_only_adot_excess.tolist(),
            final_excess_radius_gap_m=da_excess[-1].tolist(),final_gap_rate_m_s=v[-1].tolist(),
            final_gap_error_m=e[-1].tolist(),final_baseline_gap_error_m=be[-1].tolist())
        representative.append(summary)
        curves=[]
        for i in range(0,len(t),360):
            for j in range(20):curves.append(dict(time_s=float(t[i]),gap=j+1,
                delta_B_m2_kg=float(db[j]),excess_radius_gap_m=float(da_excess[i,j]),
                gap_rate_m_s=float(v[i,j]),baseline_gap_rate_m_s=float(bv[i,j]),
                gap_error_m=float(e[i,j]),baseline_gap_error_m=float(be[i,j])))
        write_csv(DATA/f'representative_{f}_hourly.csv',curves)
    write_json(DATA/'mechanism_summary.json',dict(representatives=representative,
        theory_max_absolute_error_m_s=max(abs(r['error_m_s']) for r in theory),
        theory_max_relative_error=max(r['relative_error'] for r in theory if r['relative_error'] is not None),
        comparison='Excess semi-major-axis gap subtracts uniform-B NRLMSIS reference at same time',
        B_initial_causal_term='-B0*epsilon_j*rho_j(0)*sqrt(mu*a0), exactly at common initial states'))
    # A sparse same-coordinate density audit, separate from old full-history audit.
    _,backend,_=make_environment(ROOT/spec['epoch_environment'],ROOT/spec['weather'],0)
    audits=[]
    for stage in ('initial_7d','extended_30d','combination_7d'):
        path=OUT/stage
        if not (path/'completion.json').exists():continue
        st=np.load(path/'time_s.npy')
        sa=np.load(path/'a_m.npy',mmap_mode='r');su=np.load(path/'lambda_rad.npy',mmap_mode='r')
        sr=np.load(path/'density_kg_m3.npy',mmap_mode='r')
        # Six-hour samples alone hit120s cache boundaries; also check interior
        # 10/50/90s offsets, avoiding the40/80s cubic interpolation nodes.
        anchors=np.arange(0,len(st),2160)
        candidates=np.r_[np.concatenate([anchors+offset for offset in (0,1,5,9)]),len(st)-1]
        ix=np.unique(candidates[candidates<len(st)])
        direct=backend(st[ix,None,None],sa[ix],su[ix]);err=np.abs(sr[ix]/direct-1)
        audits.append(dict(stage=stage,max_relative_error=float(err.max()),points=int(err.size)))
    worst=max(r['max_relative_error'] for r in audits)
    write_json(DATA/'density_spot_audit.json',dict(stages=audits,max_relative_error=worst,points=sum(r['points'] for r in audits),
        sampling='Every6h plus10/50/90s offsets inside120s cache blocks, and final state',
        passed=bool(worst<.001),criterion='Same-coordinate relative density discrepancy <0.001',
        caveat='Sparse point check, not full integrated propagation convergence; previous120s cache error at60d about1m'))
    if worst>=.001:raise ValueError('Dispersion trajectory density cache mismatch')
    combo=OUT/'combination_7d'
    if (combo/'completion.json').exists():
        ca=np.load(combo/'a_m.npy',mmap_mode='r')[:,0];cu=np.load(combo/'lambda_rad.npy',mmap_mode='r')[:,0]
        ce,_,_=gap_kinematics(ca,cu,reference_radius_m=cfg.derived.reference_radius_m,
            target_spacing_rad=cfg.derived.target_spacing_rad,mu_m3_s2=cfg.constants.mu_m3_s2)
        single=[]
        for factor,sigma in [('a',1.),('B',.01)]:
            k=next(i for i,p in enumerate(profiles) if p['factor']==factor and p['nominal_sigma']==sigma and p['seed']==0)
            single.append(gap_kinematics(a[:,k],u[:,k],reference_radius_m=cfg.derived.reference_radius_m,
                target_spacing_rad=cfg.derived.target_spacing_rad,mu_m3_s2=cfg.constants.mu_m3_s2)[0])
        write_json(DATA/'combination_comparison.json',dict(superposition_deviation_max_m=float(np.abs(ce-single[0]-single[1]+be).max()),
            caution='One independent-pattern pair; local sign cancellations allowed; no ensemble claim'))
    print([(g['factor'],g['sigma'],g['horizons'],g['provisional_rank']) for g in groups],flush=True)

if __name__=='__main__':main()
