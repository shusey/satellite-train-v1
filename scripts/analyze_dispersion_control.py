"""Reuse OFF screening histories; evaluate new fixed-controller profiles offline."""
import argparse,csv,gc,hashlib,json,sys
from pathlib import Path
import numpy as np
import yaml
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from run_dispersion_control import ROOT,OUT,DATA,configuration,frozen_profiles
from sat_string.dispersion_control import derived_histories
from sat_string.dispersion_transfer_metrics import evaluate_transfer
from sat_string.io import write_json
from sat_string.control_cost import clean_json
from sat_string.long_horizon import decompose_q
from sat_string.dynamics import orbital_rates,mean_motion_rad_s
from sat_string.residual import residual_analysis
from sat_string.kd_analysis import enhance_residual,threshold_events,return_gates


def csv_rows(path,rows):
    if not rows:return
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with path.open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=fields);w.writeheader()
        for r in rows:w.writerow({k:json.dumps(v,ensure_ascii=False) if isinstance(v,(dict,list)) else v for k,v in r.items()})


def forecast_summary(rows):
    """Separate genuine upward extrapolation from an already-satisfied dwell gate."""
    matched=[r for r in rows if r['prediction_valid_on_observed_window'] and not r['right_censored']
             and r['prediction_target_matches_actual'] is True]
    result={}
    for status in sorted({r['prediction']['status'] for r in matched}):
        subset=[r for r in matched if r['prediction']['status']==status]
        result[status]=dict(count=len(subset),
            median_absolute_error_s=float(np.median([abs(r['prediction_error_actual_minus_predicted_s']) for r in subset])),
            median_absolute_relative_error=float(np.median([abs(r['relative_prediction_error']) for r in subset])),
            median_actual_wait_s=float(np.median([r['actual_next_wait_s'] for r in subset])),
            commands_delayed_by_gate=sum(r['prediction']['predicted_command_s']-r['prediction']['predicted_crossing_s']>10 for r in subset))
    return result


def load_member(paths,profile_id,case=None):
    result={};events=[]
    for path in paths:
        if not (path/'completion.json').exists():raise ValueError('Incomplete source')
        meta=json.loads((path/'metadata.json').read_text());ps=meta['profiles']
        j=next(i for i,p in enumerate(ps) if p['id']==profile_id and (case is None or p.get('case')==case))
        t=np.load(path/'time_s.npy');drop=1 if result else 0
        if result and result['time_s'][-1]!=t[0]:raise ValueError('Noncontiguous source')
        for key in ['time_s','a_m','lambda_rad','density_kg_m3','area_m2','mode']:
            if key=='time_s':values=t[drop:]
            elif (path/f'{key}.npy').exists():values=np.load(path/f'{key}.npy',mmap_mode='r')[drop:,j].copy()
            else:values=np.full((len(t)-drop,21),.01 if key=='area_m2' else 0,dtype=float if key=='area_m2' else 'int8')
            result[key]=np.concatenate((result[key],values)) if key in result else values
        if (path/'checkpoint.json').exists():events=json.loads((path/'checkpoint.json').read_text())['events'][j]
    return result,events


def sources(days,case,p,spec,on=True):
    if on or case=='AB':
        paths=[OUT/'batch_30d']
        if days==60:paths.append(OUT/'batch_60d')
        return paths,p['id'],case if on else 'AB_OFF'
    batch=OUT/'batch_30d'
    if (batch/'metadata.json').exists() and any(x.get('case')==case+'_OFF'
            for x in json.loads((batch/'metadata.json').read_text())['profiles']):
        paths=[batch]+([OUT/'batch_60d'] if days==60 else [])
        return paths,p['id'],case+'_OFF'
    paths=[ROOT/spec['screening_output']/'initial_7d',ROOT/spec['screening_output']/'extended_30d']
    if days==60:paths.append(OUT/'OFF_A_B_60d')
    return paths,p['source_profile'],None


def build_histories(raw,cfg):
    return derived_histories(*[raw[k] for k in ('time_s','a_m','lambda_rad','area_m2','mode','density_kg_m3')],cfg)


def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--days',type=int,choices=[30,60],required=True)
    args=parser.parse_args();spec,cfg=configuration(args.days);profiles=frozen_profiles(spec)[:3]
    provenance=ROOT/spec['baseline_run']/'density_environment.json'
    expected_env=json.loads(provenance.read_text()) if provenance.exists() else dict(environment=yaml.safe_load((ROOT/spec['environment']).read_text()))
    weather_hash=hashlib.sha256((ROOT/spec['weather']).read_bytes()).hexdigest()
    rows=[];input_checks=[];decisions=[];target=DATA/f'{args.days}d';target.mkdir(parents=True,exist_ok=False)
    for p in profiles:
        case=p['case'];paths,ident,label=sources(args.days,case,p,spec)
        raw,events=load_member(paths,ident,label);on=build_histories(raw,cfg)
        lpaths,lid,llabel=sources(args.days,case,p,spec,False)
        lraw,lev=load_member(lpaths,lid,llabel);off=build_histories(lraw,cfg)
        if lev or np.any(off['mode']!=0) or np.any(off['area_m2']!=cfg.spacecraft.area_low_m2):
            raise ValueError('OFF unexpectedly commanded or changed area/mode')
        # Mechanical equality of profile/environment/controller; only ON flag differs.
        for path,profile_id,which in [(z,ident,label) for z in paths]+[(z,lid,llabel) for z in lpaths]:
            meta=json.loads((path/'metadata.json').read_text());rawcfg=meta['configuration']
            rawcfg['simulation']['duration_s']=cfg.simulation.duration_s
            if rawcfg!=cfg.to_dict():raise ValueError('Physical/controller configuration mismatch')
            source=next(x for x in meta['profiles'] if x['id']==profile_id and (which is None or x.get('case')==which))
            if source['delta_a_m']!=p['delta_a_m'] or source['B_multiplier']!=p['B_multiplier']:raise ValueError('Profile mismatch')
            source_control=source.get('control_enabled',meta.get('control_enabled',False))
            if source_control!=(which is not None and not which.endswith('_OFF')):raise ValueError('Unexpected control flag')
            environment=json.loads((path/'environment.json').read_text())
            if (environment.get('environment',environment.get('spec'))!=expected_env['environment']
                or environment['weather_sha256']!=weather_hash or environment['cache_seconds']!=spec['cache_seconds']):
                raise ValueError('Environment/weather/cache mismatch')
        if on['time_s'][-1]!=args.days*86400:raise ValueError('Wrong horizon')
        input_checks.append(dict(case=case,on_sources=[str(z.relative_to(ROOT)) for z in paths],
            off_sources=[str(z.relative_to(ROOT)) for z in lpaths],configuration_identical=True,
            initial_a_lambda_exact=True,profiles_exact=True,difference='control_enabled only',
            environment_exact=True,weather_sha256=weather_hash,cache_seconds=spec['cache_seconds'],
            OFF_prefix_reused=any((ROOT/spec['screening_output']) in z.parents for z in lpaths),OFF_AB_new_required=case=='AB'))
        metrics=evaluate_transfer(on,cfg,events,off);lowmetrics=evaluate_transfer(off,cfg,[],off)
        metrics.update(case=case,days=args.days,control_enabled=True)
        lowmetrics.update(case=case,days=args.days,control_enabled=False)
        write_json(target/f'{case}_metrics.json',clean_json(metrics));write_json(target/f'{case}_OFF_metrics.json',clean_json(lowmetrics))
        rows.extend([lowmetrics,metrics]);csv_rows(target/f'{case}_events.csv',events)
        triggers=threshold_events(on,cfg,events);gates=return_gates(on,cfg,events)
        for r in triggers:
            k=int(np.searchsorted(on['time_s'],r['command_time_s']));j=r['satellite_id']-1
            incident=[x for x in (j-1,j) if 0<=x<20]
            r.update(max_error_at_command_m=float(np.abs(on['gap_error_m'][k]).max()),
                incident_errors_m=on['gap_error_m'][k,incident].tolist(),incident_rates_m_s=on['gap_rate_m_s'][k,incident].tolist())
        csv_rows(target/f'{case}_HIGH_triggers.csv',triggers);csv_rows(target/f'{case}_return_gates.csv',gates)
        if triggers:
            first=triggers[0];at=first['command_time_s'];j=first['satellite_id']-1
            ix=np.flatnonzero((on['time_s']>=max(0,at-3600))&(on['time_s']<=at+7200))
            pp,dd=decompose_q(on['gap_error_m'][ix],on['gap_rate_m_s'][ix],cfg.derived.controller_kp_s2_inv,cfg.derived.controller_kd_s_inv)
            np.savez_compressed(OUT/f'{case}_first_trigger_{args.days}d.npz',time_s=on['time_s'][ix],
                satellite_id=j+1,command_time_s=at,qP=pp[:,j],qD=dd[:,j],q=pp[:,j]+dd[:,j],
                gap_error_m=on['gap_error_m'][ix],gap_rate_m_s=on['gap_rate_m_s'][ix],mode=on['mode'][ix,j])
        residual,states=residual_analysis(on,cfg,events,groups=[],low_arrays=off)
        residual=enhance_residual(residual,states,cfg)
        np.savez_compressed(OUT/f'{case}_end_states_{args.days}d.npz',**states)
        for section in ['local_cycles','quiet_cycles']:csv_rows(target/f'{case}_{section}.csv',residual[section])
        diagnostics=dict(residual=residual['summary'],quiet_packets=len(residual['quiet_cycles']),
            command_contribution_counts={kind:sum(x['command_contribution_class']==kind for x in triggers)
                for kind in ('P_dominant','D_dominant','mixed')},
            quiet_onset_contribution_counts={kind:sum(x['command_contribution_class']==kind and x['quiet_packet_onset'] for x in triggers)
                for kind in ('P_dominant','D_dominant','mixed')},
            later_in_packet_contribution_counts={kind:sum(x['command_contribution_class']==kind and not x['quiet_packet_onset'] for x in triggers)
                for kind in ('P_dominant','D_dominant','mixed')},
            first_HIGH_trigger=triggers[0] if triggers else None,
            returns_dwell_limited=sum(x['observed_hold_s']==300 for x in gates),
            returns_after_dwell=sum(x['observed_hold_s'] is not None and x['observed_hold_s']>300 for x in gates),
            first_off_gate_blocked=sum(x['off_not_met_at_gate'] is True for x in gates),
            min_same_satellite_HIGH_start_interval_s=min((r['next_same_satellite_start_interval_s'] for r in residual['local_cycles']
                if r['next_same_satellite_start_interval_s'] is not None),default=None),
            cycles_at_840s=sum(r['next_same_satellite_start_interval_s']==840 for r in residual['local_cycles']),
            completed_local_Vres_median_m_s=float(np.median([r['V_res_all_m_s'] for r in residual['local_cycles']])) if residual['local_cycles'] else None,
            completed_local_incident_Vres_median_m_s=float(np.median([r['V_res_seeds_m_s'] for r in residual['local_cycles']])) if residual['local_cycles'] else None,
            own_extra_loss_per_cycle_median_m=float(np.median([r['own_extra_loss_active_increment_m'] for r in residual['local_cycles']])) if residual['local_cycles'] else None,
            initial_max_gap_rate_m_s=float(np.abs(on['gap_rate_m_s'][0]).max()),
            final_max_gap_rate_m_s=float(np.abs(on['gap_rate_m_s'][-1]).max()),
            same_satellite_HIGH_intervals_csv=f'{case}_local_cycles.csv',
            quiet_valid_origin_matched_predictions=sum(r['prediction_valid_on_observed_window'] and not r['right_censored']
                and r['prediction_target_matches_actual'] is True for r in residual['quiet_cycles']),
            quiet_forecasts_invalid_D_change=sum('D_not_constant_within_analysis_tolerance' in r['invalid_reasons'] for r in residual['quiet_cycles']),
            quiet_prediction_status_summary=forecast_summary(residual['quiet_cycles']),
            quiet_forecasts_csv=f'{case}_quiet_cycles.csv')
        write_json(target/f'{case}_diagnostics.json',clean_json(diagnostics))
        # Small postprocessing cache, local only. Plot max/RMS remain sampled envelopes.
        event_indices=np.searchsorted(on['time_s'],[x[z] for x in events for z in ('start_time','complete_time')])
        event_indices=event_indices[event_indices<len(on['time_s'])]
        ix=np.unique(np.r_[np.arange(0,len(on['time_s']),6),len(on['time_s'])-1,event_indices,
            np.argmax(np.max(np.abs(on['gap_error_m']),axis=1)),np.argmax(np.abs(on['gap_rate_m_s']),axis=0)])
        high_count=np.sum(on['mode']==2,axis=1);nonlow_count=np.sum(on['mode']!=0,axis=1)
        changes=np.unique(np.r_[0,np.flatnonzero((np.diff(high_count)!=0)|(np.diff(nonlow_count)!=0))+1,len(on['time_s'])-1])
        qp,qd=decompose_q(on['gap_error_m'][ix],on['gap_rate_m_s'][ix],cfg.derived.controller_kp_s2_inv,cfg.derived.controller_kd_s_inv)
        np.savez_compressed(OUT/f'{case}_plot_{args.days}d.npz',time_s=on['time_s'][ix],
            on_error=on['gap_error_m'][ix],off_error=off['gap_error_m'][ix],on_rate=on['gap_rate_m_s'][ix],
            a_m=on['a_m'][ix],off_a_m=off['a_m'][ix],rho=on['density_kg_m3'][ix],area=on['area_m2'][ix],mode=on['mode'][ix],qP=qp,qD=qd,
            activity_time_s=on['time_s'][changes],high_count=high_count[changes],nonLOW_count=nonlow_count[changes])
        gaprows=[dict(gap=i+1,max_error_m=float(np.abs(on['gap_error_m'][:,i]).max()),
            RMS_m=metrics['gap_rms_m'][i],final_error_m=float(on['gap_error_m'][-1,i]),final_rate_m_s=float(on['gap_rate_m_s'][-1,i])) for i in range(20)]
        csv_rows(target/f'{case}_gaps.csv',gaprows)
        # Three representative gaps, hourly, retaining signs and the persistent
        # forcing at the SAME instantaneous state. Counterfactual equal-B rates
        # isolate the direct multiplier contribution, not an alternative run.
        worst=int(np.argmax(np.max(np.abs(on['gap_error_m']),axis=0)))
        indices=np.unique(np.r_[np.arange(0,len(on['time_s']),360),len(on['time_s'])-1])
        aa=on['a_m'][indices];rho=on['density_kg_m3'][indices];area=on['area_m2'][indices]
        factors=np.asarray(p['B_multiplier']);mu=cfg.constants.mu_m3_s2;s=cfg.spacecraft
        da,n=orbital_rates(aa,rho,area,mu_m3_s2=mu,drag_coefficient=s.drag_coefficient*factors,mass_kg=s.mass_kg)
        low_da,_=orbital_rates(aa,rho,np.full_like(area,s.area_low_m2),mu_m3_s2=mu,drag_coefficient=s.drag_coefficient*factors,mass_kg=s.mass_kg)
        equal_da,_=orbital_rates(aa,rho,np.full_like(area,s.area_low_m2),mu_m3_s2=mu,drag_coefficient=s.drag_coefficient,mass_kg=s.mass_kg)
        acceleration=cfg.derived.reference_radius_m*np.diff(-1.5*n/aa*da,axis=1)
        mechanisms=[]
        for g in sorted({0,9,worst}):
            for z,k in enumerate(indices):
                mechanisms.append(dict(time_s=float(on['time_s'][k]),gap=g+1,delta_B_fraction=float(np.diff(factors)[g]),
                    radius_difference_m=float(aa[z,g+1]-aa[z,g]),
                    radius_rate_difference_m_s=float(da[z,g+1]-da[z,g]),
                    same_state_allLOW_radius_rate_difference_m_s=float(low_da[z,g+1]-low_da[z,g]),
                    same_state_equal_B_radius_rate_difference_m_s=float(equal_da[z,g+1]-equal_da[z,g]),
                    direct_B_contribution_to_radius_rate_difference_m_s=float(np.diff(low_da-equal_da,axis=1)[z,g]),
                    gap_acceleration_m_s2=float(acceleration[z,g]),gap_rate_m_s=float(on['gap_rate_m_s'][k,g]),
                    gap_error_m=float(on['gap_error_m'][k,g]),both_satellites_stable_LOW=bool(np.all(on['mode'][k,g:g+2]==0))))
        csv_rows(target/f'{case}_mechanism_hourly.csv',mechanisms)
        csv_rows(target/f'{case}_satellites.csv',[dict(satellite_id=j+1,delta_a_m=p['delta_a_m'][j],epsilon=p['epsilon'][j],
            switches=metrics['switches_per_satellite'][j],cycles=metrics['cycles_per_satellite'][j],
            final_mode=metrics['final_modes'][j],final_altitude_m=metrics['final_altitude_m'][j],
            extra_altitude_loss_m=metrics['extra_loss_per_satellite_m'][j],high_stable_s=metrics['high_stable_s_per_satellite'][j],
            high_equivalent_s=metrics['high_equivalent_s_per_satellite'][j]) for j in range(21)])
        safe=metrics['first_label_order_reversal_s'] is None and metrics['min_signed_reference_gap_m']>0
        active=metrics['segments'][-1]['transitions']>0
        # Candidate only. A research review must judge abnormal cost/failure;
        # there is no invented mission-specific error/cost cutoff here.
        extend=safe and active
        decisions.append(dict(case=case,extend_to_60d=extend,safe=safe,finite=True,
            latter_segment_transitions=metrics['segments'][-1]['transitions'],max_error_m=metrics['max_error_m'],
            reason='positive signed gaps, finite states, ongoing activity; cost and precision require review' if extend else 'stop: check safety/activity'))
        print(case,args.days,metrics['max_error_m'],metrics['rms_error_m'],metrics['transitions'],metrics['extra_mean_altitude_loss_m'],flush=True)
        del raw,on,lraw,off,states,residual;gc.collect()
    write_json(target/'comparison.json',clean_json(dict(cases=rows,input_checks=input_checks)))
    smallkeys=['case','days','control_enabled','max_error_m','rms_error_m','final_max_error_m','max_gap_rate_m_s','transitions','completed_cycles','last_switch_time_s','extra_mean_altitude_loss_m','high_stable_total_s','high_equivalent_total_s','first_label_order_reversal_s']
    csv_rows(target/'comparison.csv',[{k:r[k] for k in smallkeys} for r in rows])
    if args.days==30:write_json(DATA/'extension_decision.json',dict(extend_to_60d=False,reviewed=False,
        all_cases_eligible=all(x['extend_to_60d'] for x in decisions),cases=decisions,
        interpretation='normal order and meaningful ongoing activity; no tuning or asymptotic stability inference'))


if __name__=='__main__':main()
