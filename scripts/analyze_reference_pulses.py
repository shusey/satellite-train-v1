"""Read the frozen 60d Reference runs; do not call the simulation core."""
from pathlib import Path
import argparse, csv, hashlib, json, sys
import numpy as np
import yaml

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'src'))
from sat_string.residual import cycle_pairs
from sat_string.pulse_diagnostics import (acceleration_terms, cumulative_trapezoid,
    local_control_split, node_map, reversal_flags)

SOURCE=ROOT/'outputs/reference_dispersion_20261012'
OLD=ROOT/'docs/reference_dispersion/data/60d'
FROZEN=ROOT/'docs/reference_dispersion/data/profiles.json'


def write_json(path,obj):
    path.write_text(json.dumps(obj,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')


def write_csv(path,rows):
    if not rows:return
    with path.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)


def read_csv(path):
    with path.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))


def quantiles(x):
    x=np.asarray(x,dtype=float);x=x[np.isfinite(x)]
    return dict(n=len(x),min=float(x.min()),p10=float(np.quantile(x,.1)),
        median=float(np.median(x)),p90=float(np.quantile(x,.9)),max=float(x.max())) if len(x) else dict(n=0)


def inputs_manifest():
    # Hash small inputs; size/mtime of large arrays are checked before and after.
    result={}
    for p in sorted(set(list(SOURCE.glob('batch_*d/*'))+list(OLD.glob('*'))+[FROZEN])):
        if p.is_file():
            s=p.stat();v=dict(size=s.st_size,mtime_ns=s.st_mtime_ns)
            if s.st_size<20_000_000:v['sha256']=hashlib.sha256(p.read_bytes()).hexdigest()
            result[str(p.relative_to(ROOT))]=v
    return result


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',default='outputs/reference_pulse_mechanism_20261017_v3')
    args=ap.parse_args();out=ROOT/args.output;out.mkdir(parents=True,exist_ok=False)
    before=inputs_manifest();write_json(out/'input_manifest.json',before)
    cfg=yaml.safe_load((SOURCE/'batch_30d/normalized_config.yaml').read_text())
    prof=json.loads(FROZEN.read_text(encoding='utf-8'))
    all_events=json.loads((SOURCE/'batch_60d/checkpoint.json').read_text())['events']
    c=cfg['constants'];s=cfg['spacecraft'];d=cfg['derived'];ctrl=cfg['controller']
    kp=d['controller_kp_s2_inv'];kd=d['controller_kd_s_inv'];r=d['reference_radius_m'];mu=c['mu_m3_s2']
    t=np.concatenate([np.load(SOURCE/'batch_30d/time_s.npy'),np.load(SOURCE/'batch_60d/time_s.npy')[1:]])
    assert len(t)==518401 and np.all(np.diff(t)==10)
    summary=dict(definition='observed-trajectory force attribution; not shorter-pulse counterfactual',
        minimum_equivalent_s=ctrl['slew_time_s']+ctrl['dwell_time_s'],
        minimum_cycle_s=2*ctrl['slew_time_s']+ctrl['dwell_time_s'],
        minimum_restart_s=2*(ctrl['slew_time_s']+ctrl['dwell_time_s']),cases={})
    segment_rows=[];examples=[]
    for ic,case in enumerate(['A','B','AB']):
        print(f'Reading {case}',flush=True)
        def joined(key):
            x=np.load(SOURCE/f'batch_30d/{key}.npy',mmap_mode='r')[:,ic]
            y=np.load(SOURCE/f'batch_60d/{key}.npy',mmap_mode='r')[1:,ic]
            return np.concatenate((x,y))
        a=joined('a_m');rho=joined('density_kg_m3');area=joined('area_m2');u=joined('lambda_rad')
        e=r*np.diff(u,axis=1)-cfg['formation']['target_spacing_m'];del u
        v=r*np.diff(np.sqrt(mu/a**3),axis=1)
        qp=kp*node_map(e);qd=kd*node_map(v)
        eps=np.array(prof['profiles'][ic]['epsilon'])
        terms=acceleration_terms(a,rho,area,eps,r_ref=r,mu=mu,cd=s['drag_coefficient'],mass=s['mass_kg'],area_low=s['area_low_m2'])
        integrals={k:cumulative_trapezoid(x,t) for k,x in terms.items()}
        # State-level attribution error is reported; it is not hidden in a tolerance.
        predicted=v[:1]+np.diff(sum(integrals.values()),axis=1)
        global_closure=float(np.max(abs(predicted-v)));del predicted
        pairs=cycle_pairs(all_events[ic],t[-1]);local=[];gaps=[]
        prior={ (int(z['satellite_id']),int(z['local_cycle_number'])):z for z in read_csv(OLD/f'{case}_local_cycles.csv') }
        for pair in pairs:
            rise,back=pair['rise'],pair['return_event']
            if back is None or back['complete_time']>t[-1]:continue
            j=pair['satellite_id']-1;num=pair['local_cycle_number'];old=prior[(j+1,num)]
            k0=int(np.searchsorted(t,rise['start_time']));k1=int(np.searchsorted(t,back['complete_time']))
            assert t[k0]==rise['start_time'] and t[k1]==back['complete_time']
            impulse={k:x[k1]-x[k0] for k,x in integrals.items()}
            own,neighbor=local_control_split(impulse['control'],j)
            low=np.diff(impulse['low']);mismatch=np.diff(impulse['mismatch'])
            incident=[g for g in (j-1,j) if 0<=g<v.shape[1]]
            observed=v[k1]-v[k0];closure=observed-low-mismatch-own-neighbor
            teq=float(old['high_equivalent_time_s'])
            flags=reversal_flags(v[k0,incident],v[k1,incident],own[incident])
            reduction=float(np.max(abs(v[k0,incident]))-np.max(abs(v[k1,incident])))
            p0=float(qp[k0,j]);p1=float(qp[k1,j]);d0=float(qd[k0,j]);d1=float(qd[k1,j])
            below=float(old['first_q_below_off_s']) if old['first_q_below_off_s'] else None
            # Integrate after first observed <=h_off; recrossings and dwell are retained.
            own_after_off=None
            if below is not None:
                kb=int(np.searchsorted(t,below))
                own_after_off=float((integrals['control'][k1,j]-integrals['control'][kb,j]))
            row=dict(case=case,satellite_id=j+1,cycle=num,start_s=float(t[k0]),end_s=float(t[k1]),
                high_equivalent_s=teq,minimum_pulse=bool(abs(teq-summary['minimum_equivalent_s'])<1e-8),
                return_driver=old['return_driver'],first_q_below_off_s=below,
                qP_before=p0,qP_after=p1,qD_before=d0,qD_after=d1,
                qD_own_control=float(kd*node_map(own)[j]),
                qD_neighbor_control=float(kd*node_map(neighbor)[j]),
                qD_common_LOW=float(kd*node_map(low)[j]),qD_B_mismatch=float(kd*node_map(mismatch)[j]),
                V_incident_before_m_s=float(np.max(abs(v[k0,incident]))),
                V_incident_after_m_s=float(np.max(abs(v[k1,incident]))),
                RMS_incident_before_m_s=float(np.sqrt(np.mean(v[k0,incident]**2))),
                RMS_incident_after_m_s=float(np.sqrt(np.mean(v[k1,incident]**2))),
                incident_max_velocity_reduction_m_s=reduction,
                incident_abs_growth_count=int(flags['actual_abs_growth'].sum()),
                incident_sign_reversal_count=int(flags['actual_sign_reversal'].sum()),
                own_alone_reversal_count=int(flags['own_alone_crosses_zero'].sum()),
                own_excess_radius_loss_m=float(np.trapezoid(rho[k0:k1+1,j]*(s['drag_coefficient']/s['mass_kg'])*(1+eps[j])*(area[k0:k1+1,j]-s['area_low_m2'])*np.sqrt(mu*a[k0:k1+1,j]),t[k0:k1+1])),
                own_excess_g_impulse_m_s=float(impulse['control'][j]),
                own_excess_g_impulse_after_first_off_m_s=own_after_off,
                next_own_HIGH_s=float(old['actual_next_command_s']) if old['actual_next_command_s'] else None,
                next_own_wait_s=float(old['actual_next_wait_s']) if old['actual_next_wait_s'] else None,
                closure_max_m_s=float(np.max(abs(closure[incident]))))
            # A local LOW satellite is not necessarily in an all-LOW formation.
            # Attribute its D regrowth to neighbors versus persistent LOW drag.
            follow={key:None for key in ['next_delta_qP','next_delta_qD','next_qD_common_LOW',
                'next_qD_B_mismatch','next_qD_neighbor_control','next_qD_closure']}
            if row['next_own_HIGH_s'] is not None:
                k2=int(np.searchsorted(t,row['next_own_HIGH_s']))
                increments={key:x[k2]-x[k1] for key,x in integrals.items()}
                assert abs(increments['control'][j])<1e-12,'Own area changed while LOW'
                follow.update(next_delta_qP=float(qp[k2,j]-qp[k1,j]),next_delta_qD=float(qd[k2,j]-qd[k1,j]),
                    next_qD_common_LOW=float(kd*node_map(np.diff(increments['low']))[j]),
                    next_qD_B_mismatch=float(kd*node_map(np.diff(increments['mismatch']))[j]),
                    next_qD_neighbor_control=float(kd*node_map(np.diff(increments['control']))[j]))
                follow['next_qD_closure']=follow['next_delta_qD']-sum(follow[key] for key in
                    ['next_qD_common_LOW','next_qD_B_mismatch','next_qD_neighbor_control'])
            row.update(follow)
            local.append(row)
            for g in incident:
                gaps.append(dict(case=case,satellite_id=j+1,cycle=num,gap_id=g+1,start_s=float(t[k0]),end_s=float(t[k1]),
                    e_before_m=float(e[k0,g]),e_after_m=float(e[k1,g]),v_before_m_s=float(v[k0,g]),v_after_m_s=float(v[k1,g]),
                    delta_a_before_m=float(a[k0,g+1]-a[k0,g]),delta_a_after_m=float(a[k1,g+1]-a[k1,g]),
                    dv_common_LOW_m_s=float(low[g]),dv_B_mismatch_m_s=float(mismatch[g]),
                    dv_own_control_m_s=float(own[g]),dv_neighbor_control_m_s=float(neighbor[g]),
                    closure_m_s=float(closure[g]),worsening_before=bool(e[k0,g]*v[k0,g]>0),
                    worsening_after=bool(e[k1,g]*v[k1,g]>0)))
        write_csv(out/f'{case}_pulse_cycles.csv',local);write_csv(out/f'{case}_incident_gaps.csv',gaps)
        quiet=[]
        for old in read_csv(OLD/f'{case}_quiet_cycles.csv'):
            if not old['actual_next_command_s']:continue
            t0=float(old['completion_time_s']);t1=float(old['actual_next_command_s'])
            k0=int(np.searchsorted(t,t0));k1=int(np.searchsorted(t,t1));j=int(old['actual_next_satellite_id'])-1
            dl=kd*node_map(np.diff(integrals['low'][k1]-integrals['low'][k0]))[j]
            db=kd*node_map(np.diff(integrals['mismatch'][k1]-integrals['mismatch'][k0]))[j]
            dc=kd*node_map(np.diff(integrals['control'][k1]-integrals['control'][k0]))[j]
            assert abs(dc)<1e-15, 'Quiet interval contains actuation'
            slope=kp/kd*qd[k0,j]
            pred=t0+(ctrl['threshold_on_m_s2']-qp[k0,j]-qd[k0,j])/slope if slope>0 else None
            quiet.append(dict(case=case,start_s=t0,end_s=t1,duration_s=t1-t0,trigger_satellite=j+1,
                qP_start=float(qp[k0,j]),qP_end=float(qp[k1,j]),qD_start=float(qd[k0,j]),qD_end=float(qd[k1,j]),
                delta_qP=float(qp[k1,j]-qp[k0,j]),delta_qD=float(qd[k1,j]-qd[k0,j]),
                delta_qD_common_LOW=float(dl),delta_qD_B_mismatch=float(db),
                qD_closure=float(qd[k1,j]-qd[k0,j]-dl-db),
                frozen_D_prediction_for_actual_satellite_s=pred,
                frozen_D_error_s=t1-pred if pred is not None else None,
                existing_first_arrival_valid=old['prediction_valid_on_observed_window']=='True',
                existing_first_arrival_error_s=float(old['validated_prediction_error_s']) if old['validated_prediction_error_s'] else None,
                max_velocity_start_m_s=float(np.max(abs(v[k0]))),max_velocity_end_m_s=float(np.max(abs(v[k1]))),
                worsening_fraction_start=float(np.mean(e[k0]*v[k0]>0))))
        write_csv(out/f'{case}_quiet_intervals.csv',quiet)
        # Predeclared examples: first local cycle, first completed all-LOW interval,
        # first local cycle and all-LOW interval starting on/after day 30.
        selected=[('first_pulse',local[0]),('late_pulse',next(x for x in local if x['start_s']>=30*86400))]
        qselected=[('first_quiet',quiet[0]),('late_quiet',next(x for x in quiet if x['start_s']>=30*86400))]
        for name,row in selected+qselected:
            j=int(row.get('satellite_id',row.get('trigger_satellite')))-1
            lo=max(0,row['start_s']-600);hi=min(t[-1],row['end_s']+600)
            indices=np.flatnonzero((t>=lo)&(t<=hi));trace=[]
            for k in indices:
                rec=dict(time_s=float(t[k]),satellite_id=j+1,qP=float(qp[k,j]),qD=float(qd[k,j]),q=float(qp[k,j]+qd[k,j]),area_m2=float(area[k,j]),rho_kg_m3=float(rho[k,j]))
                for g in (j-1,j):
                    if 0<=g<20:
                        rec[f'gap{g+1}_e_m']=float(e[k,g]);rec[f'gap{g+1}_v_m_s']=float(v[k,g]);rec[f'gap{g+1}_da_m']=float(a[k,g+1]-a[k,g])
                trace.append(rec)
            write_csv(out/f'{case}_{name}_trace.csv',trace);examples.append(dict(selection=name,**row))
        for lo,hi in [(0,7),(7,30),(30,60)]:
            rows=[x for x in local if lo*86400<=x['start_s']<hi*86400]
            segment_rows.append(dict(case=case,start_day=lo,end_day=hi,cycles=len(rows),
                minimum_pulse_fraction=float(np.mean([x['minimum_pulse'] for x in rows])),
                median_own_excess_radius_loss_m=float(np.median([x['own_excess_radius_loss_m'] for x in rows])),
                median_own_excess_g_impulse_m_s=float(np.median([x['own_excess_g_impulse_m_s'] for x in rows])),
                median_incident_velocity_before_m_s=float(np.median([x['V_incident_before_m_s'] for x in rows])),
                median_incident_velocity_after_m_s=float(np.median([x['V_incident_after_m_s'] for x in rows])),
                velocity_max_reduced_fraction=float(np.mean([x['incident_max_velocity_reduction_m_s']>0 for x in rows]))))
        valid=[x['existing_first_arrival_error_s'] for x in quiet if x['existing_first_arrival_valid'] and x['existing_first_arrival_error_s'] is not None]
        mismatch_abs=[abs(x['delta_qD_B_mismatch']) for x in quiet]
        low_abs=[abs(x['delta_qD_common_LOW']) for x in quiet]
        summary['cases'][case]=dict(cycles=len(local),incident_gap_records=len(gaps),quiet_intervals=len(quiet),
            minimum_pulse_count=sum(x['minimum_pulse'] for x in local),
            incident_max_velocity_reduced_fraction=float(np.mean([x['incident_max_velocity_reduction_m_s']>0 for x in local])),
            physical_gap_sign_reversal_fraction=float(np.mean([x['v_before_m_s']*x['v_after_m_s']<0 for x in gaps])),
            physical_gap_abs_speed_grew_fraction=float(np.mean([abs(x['v_after_m_s'])>abs(x['v_before_m_s']) for x in gaps])),
            own_alone_reversal_fraction=float(sum(x['own_alone_reversal_count'] for x in local)/len(gaps)),
            neighbor_control_absolute_ratio_of_sums=float(sum(abs(x['dv_neighbor_control_m_s']) for x in gaps)/sum(abs(x['dv_own_control_m_s']) for x in gaps)),
            own_excess_radius_loss_m=quantiles([x['own_excess_radius_loss_m'] for x in local]),
            own_excess_g_impulse_m_s=quantiles([x['own_excess_g_impulse_m_s'] for x in local]),
            own_impulse_after_first_off_ratio=float(sum(x['own_excess_g_impulse_after_first_off_m_s'] or 0 for x in local)/sum(x['own_excess_g_impulse_m_s'] for x in local)),
            V_incident_before_m_s=quantiles([x['V_incident_before_m_s'] for x in local]),
            V_incident_after_m_s=quantiles([x['V_incident_after_m_s'] for x in local]),
            pulse_closure_max_m_s=max(x['closure_max_m_s'] for x in local),
            cumulative_60d_closure_max_m_s=global_closure,
            quiet_delta_qD_B_abs=quantiles(mismatch_abs),quiet_delta_qD_common_LOW_abs=quantiles(low_abs),
            quiet_B_exceeds_common_LOW_fraction=float(np.mean(np.array(mismatch_abs)>np.array(low_abs))),
            quiet_delta_qP_positive_fraction=float(np.mean([x['delta_qP']>0 for x in quiet])),
            quiet_D_change_exceeds_P_change_fraction=float(np.mean([abs(x['delta_qD'])>abs(x['delta_qP']) for x in quiet])),
            validated_existing_prediction_abs_error_s=quantiles(np.abs(valid)),
            quiet_duration_s=quantiles([x['duration_s'] for x in quiet]))
        followed=[x for x in local if x['next_delta_qD'] is not None]
        summary['cases'][case]['local_LOW_to_next_HIGH']=dict(n=len(followed),
            neighbor_D_change_exceeds_natural_D_fraction=float(np.mean([abs(x['next_qD_neighbor_control'])>abs(x['next_qD_common_LOW']+x['next_qD_B_mismatch']) for x in followed])),
            D_change_exceeds_P_change_fraction=float(np.mean([abs(x['next_delta_qD'])>abs(x['next_delta_qP']) for x in followed])))
        print(json.dumps(summary['cases'][case],ensure_ascii=False),flush=True)
        del a,rho,area,e,v,qp,qd,terms,integrals
    write_csv(out/'interval_summary.csv',segment_rows);write_json(out/'examples.json',examples)
    after=inputs_manifest();assert before==after,'Source input files changed'
    summary['preservation']=dict(passed=True,files=len(before),large_file_check='size and nanosecond mtime; small inputs SHA256')
    write_json(out/'summary.json',summary)
    print(out,flush=True)


if __name__=='__main__':main()
