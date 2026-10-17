"""Paired seven-day A-ref D300/D150 comparison, preserving existing runs."""
import argparse,hashlib,json,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'));sys.path.insert(0,str(ROOT/'src'))
from sat_string.config import load_config,config_from_dict
from sat_string.dispersion_control import derived_histories
from sat_string.dispersion_transfer_metrics import evaluate_transfer
from sat_string.residual import residual_analysis,cycle_pairs
from sat_string.long_horizon import decompose_q
from sat_string.pulse_diagnostics import acceleration_terms,cumulative_trapezoid,local_control_split,node_map
from sat_string.dwell_intervention_metrics import distribution,followup_statistics
from sat_string.regression import compare_baseline
from sat_string.io import write_json
from sat_string.control_cost import clean_json
from analyze_reference_pulses import write_csv
from run_reference_a_dwell import SOURCE,PROFILE

DEFAULT=ROOT/'outputs/reference_a_dwell_20261017'


def snapshot():
    paths=list(SOURCE.glob('*'))+[PROFILE,ROOT/'docs/real_environment_applicability/data/active_drivers.csv']
    result={}
    for p in paths:
        if not p.is_file():continue
        stat=p.stat();record=dict(size=stat.st_size,mtime_ns=stat.st_mtime_ns)
        if stat.st_size<20_000_000:record['sha256']=hashlib.sha256(p.read_bytes()).hexdigest()
        result[str(p.relative_to(ROOT))]=record
    return result


def load_arrays(folder, member, config, days=7):
    t=np.load(folder/'time_s.npy',mmap_mode='r');stop=int(np.searchsorted(t,days*86400,side='right'))
    t=np.array(t[:stop]);assert t[-1]==days*86400
    a,u,area,mode,rho=[np.asarray(np.load(folder/f'{key}.npy',mmap_mode='r')[:stop,member]) for key in
        ['a_m','lambda_rad','area_m2','mode','density_kg_m3']]
    return derived_histories(t,a,u,area,mode,rho,config)


def events_for(folder,member,end_s):
    return [x for x in json.loads((folder/'checkpoint.json').read_text())['events'][member] if x['command_time']<=end_s]


def diagnose(arrays,config,events,low,out,label):
    t=arrays['time_s'];a=arrays['a_m'];rho=arrays['density_kg_m3'];area=arrays['area_m2']
    e=arrays['gap_error_m'];v=arrays['gap_rate_m_s'];end=float(t[-1]);n=a.shape[1]
    kp=config.derived.controller_kp_s2_inv;kd=config.derived.controller_kd_s_inv
    p,d=decompose_q(e,v,kp,kd);c=config.controller;s=config.spacecraft
    terms=acceleration_terms(a,rho,area,np.zeros(n),r_ref=config.derived.reference_radius_m,
        mu=config.constants.mu_m3_s2,cd=s.drag_coefficient,mass=s.mass_kg,area_low=s.area_low_m2)
    integ={key:cumulative_trapezoid(value,t) for key,value in terms.items()}
    metric=evaluate_transfer(arrays,config,events,low)
    residual,states=residual_analysis(arrays,config,events,groups=[],low_arrays=low)
    np.savez_compressed(out/f'{label}_endpoint_states.npz',**states)
    old={(r['satellite_id'],r['local_cycle_number']):r for r in residual['local_cycles']}
    allpairs=cycle_pairs(events,end);byid={(x['satellite_id'],x['local_cycle_number']):x for x in allpairs}
    rows=[];edges=[];stage_rows=[];stages=[];closures=[];qclosures=[];incomplete=[]
    for pair in allpairs:
        rise=pair['rise'];back=pair['return_event'];j=pair['satellite_id']-1;num=pair['local_cycle_number']
        if back is None or back['complete_time']>end:
            incomplete.append(dict(satellite_id=j+1,cycle=num,start_s=rise['start_time'],completed_by_end=False));continue
        source=old[(j+1,num)]
        k0=int(np.searchsorted(t,rise['start_time']));k1=int(np.searchsorted(t,back['complete_time']))
        assert t[k0]==rise['start_time'] and t[k1]==back['complete_time']
        impulse={key:x[k1]-x[k0] for key,x in integ.items()}
        own,neighbor=local_control_split(impulse['control'],j);natural=np.diff(impulse['low'])
        delta_v=v[k1]-v[k0];closure=delta_v-own-neighbor-natural
        incident=[g for g in [j-1,j] if 0<=g<n-1]
        nextpair=byid.get((j+1,num+1));nxt=nextpair['rise'] if nextpair else None
        next_t=nxt['start_time'] if nxt else None
        follow={key:None for key in ['next_delta_qP','next_delta_qD','next_qD_own','next_qD_neighbor','next_qD_natural',
            'next_qD_closure','q_at_LOW_dwell_release','first_q_above_on_while_LOW_s','trigger_minus_first_threshold_s']}
        if nxt is not None:
            k2=int(np.searchsorted(t,next_t));inc={key:x[k2]-x[k1] for key,x in integ.items()}
            fown,fneighbor=local_control_split(inc['control'],j)
            follow.update(next_delta_qP=float(p[k2,j]-p[k1,j]),next_delta_qD=float(d[k2,j]-d[k1,j]),
                next_qD_own=float(kd*node_map(fown)[j]),next_qD_neighbor=float(kd*node_map(fneighbor)[j]),
                next_qD_natural=float(kd*node_map(np.diff(inc['low']))[j]))
            follow['next_qD_closure']=follow['next_delta_qD']-sum(follow[key] for key in ['next_qD_own','next_qD_neighbor','next_qD_natural'])
            assert abs(follow['next_qD_own'])<1e-15
            qclosures.append(abs(follow['next_qD_closure']))
            kr=int(np.searchsorted(t,t[k1]+c.dwell_time_s))
            follow['q_at_LOW_dwell_release']=float(p[kr,j]+d[kr,j])
            hit=np.flatnonzero(p[k1:k2+1,j]+d[k1:k2+1,j]>=c.threshold_on_m_s2)
            if len(hit):
                first=float(t[k1+hit[0]]);follow['first_q_above_on_while_LOW_s']=first
                follow['trigger_minus_first_threshold_s']=float(next_t-first)
        row=dict(label=label,satellite_id=j+1,cycle=num,start_s=float(t[k0]),end_s=float(t[k1]),
            high_complete_s=rise['complete_time'],low_command_s=back['command_time'],
            high_equivalent_s=source['high_equivalent_time_s'],high_stable_s=source['high_stable_time_s'],
            minimum_pulse=bool(abs(source['high_equivalent_time_s']-(c.slew_time_s+c.dwell_time_s))<1e-9),
            return_driver=source['return_driver'],first_q_below_off_s=source['first_q_below_off_s'],
            pulse_g_impulse_m_s=float(impulse['control'][j]),
            V_incident_before_m_s=float(np.max(abs(v[k0,incident]))),V_incident_after_m_s=float(np.max(abs(v[k1,incident]))),
            e_incident_max_before_m=float(np.max(abs(e[k0,incident]))),e_incident_max_after_m=float(np.max(abs(e[k1,incident]))),
            qP_before=float(p[k0,j]),qD_before=float(d[k0,j]),qP_after=float(p[k1,j]),qD_after=float(d[k1,j]),
            next_HIGH_s=next_t,next_wait_s=float(next_t-t[k1]) if next_t is not None else None,
            next_HIGH_start_interval_s=float(next_t-t[k0]) if next_t is not None else None,
            right_censored=nxt is None,observed_followup_s=float((next_t if next_t is not None else end)-t[k1]),
            LOW_dwell_release_s=float(t[k1]+c.dwell_time_s),**follow)
        rows.append(row);closures.append(float(np.max(abs(closure[incident]))))
        # The physical state is continuous across a zero-delay command. The saved
        # attitude mode is post-command, while e/v equal the just-before state.
        stage_indices=[]
        for stage,at in [('HIGH_start',rise['start_time']),('HIGH_complete',rise['complete_time']),
                         ('LOW_command',back['command_time']),('LOW_complete',back['complete_time'])]:
            k=int(np.searchsorted(t,at));assert t[k]==at;stage_indices.append(k)
            for g in incident:
                stage_rows.append(dict(label=label,satellite_id=j+1,cycle=num,stage=stage,time_s=float(at),gap_id=g+1,
                    error_m=float(e[k,g]),rate_m_s=float(v[k,g]),e_times_v_m2_s=float(e[k,g]*v[k,g]),
                    own_qP=float(p[k,j]),own_qD=float(d[k,j]),own_q=float(p[k,j]+d[k,j]),
                    own_area_m2=float(area[k,j]),post_event_mode=int(arrays['mode'][k,j])))
        stages.append(dict(satellite_id=j+1,cycle=num,time_s=t[stage_indices],e=e[stage_indices],v=v[stage_indices],qP=p[stage_indices],qD=d[stage_indices]))
        for g in incident:
            edges.append(dict(label=label,satellite_id=j+1,cycle=num,gap_id=g+1,
                error_before_m=float(e[k0,g]),error_after_m=float(e[k1,g]),velocity_before_m_s=float(v[k0,g]),velocity_after_m_s=float(v[k1,g]),
                abs_speed_decreased=bool(abs(v[k1,g])<abs(v[k0,g])),abs_speed_increased=bool(abs(v[k1,g])>abs(v[k0,g])),
                sign_reversed=bool(v[k0,g]*v[k1,g]<0),worsening_before=bool(e[k0,g]*v[k0,g]>0),worsening_after=bool(e[k1,g]*v[k1,g]>0),
                e_v_sign_before=int(np.sign(e[k0,g]*v[k0,g])),e_v_sign_after=int(np.sign(e[k1,g]*v[k1,g])),
                own_impulse_m_s=float(own[g]),neighbor_impulse_m_s=float(neighbor[g]),natural_impulse_m_s=float(natural[g]),closure_m_s=float(closure[g])))
    write_csv(out/f'{label}_cycles.csv',rows);write_csv(out/f'{label}_incident_gaps.csv',edges)
    write_csv(out/f'{label}_cycle_stages.csv',stage_rows);write_csv(out/f'{label}_incomplete_cycles.csv',incomplete)
    np.savez_compressed(out/f'{label}_four_stage_states.npz',
        satellite_id=[x['satellite_id'] for x in stages],cycle=[x['cycle'] for x in stages],
        **{key:np.asarray([x[key] for x in stages]) for key in ['time_s','e','v','qP','qD']})
    # Keep scalar/nested summaries in JSON; state arrays remain in NPZ.
    write_csv(out/f'{label}_quiet_intervals.csv',residual['quiet_cycles'])
    write_json(out/f'{label}_residual_summary.json',residual['summary'])
    quiet=residual['quiet_cycles'];completedquiet=[x for x in quiet if not x['right_censored']]
    matched=[x for x in completedquiet if x['prediction_valid_on_observed_window'] and x['prediction_target_matches_actual']]
    followed=[x for x in rows if x['next_HIGH_s'] is not None]
    low_extra=low['a_m']-a
    stat=dict(label=label,dwell_s=c.dwell_time_s,days=end/86400,formation=metric,
        completed_cycles=len(rows),incomplete_cycles=incomplete,
        theoretical_minimum_HIGH_equivalent_s=c.slew_time_s+c.dwell_time_s,
        theoretical_minimum_cycle_s=2*c.slew_time_s+c.dwell_time_s,
        theoretical_minimum_HIGH_restart_s=2*(c.slew_time_s+c.dwell_time_s),
        minimum_pulse_count=sum(x['minimum_pulse'] for x in rows),minimum_pulse_fraction=float(np.mean([x['minimum_pulse'] for x in rows])),
        high_equivalent_s_per_cycle=distribution(x['high_equivalent_s'] for x in rows),
        pulse_velocity_contribution_m_s=distribution(x['pulse_g_impulse_m_s'] for x in rows),
        residual_abs_incident_gap_velocity_m_s=distribution(abs(x['velocity_after_m_s']) for x in edges),
        residual_max_incident_velocity_per_cycle_m_s=distribution(x['V_incident_after_m_s'] for x in rows),
        initial_max_incident_velocity_per_cycle_m_s=distribution(x['V_incident_before_m_s'] for x in rows),
        residual_worsening_fraction=float(np.mean([x['worsening_after'] for x in edges])),
        incident_velocity_decreased_fraction=float(np.mean([x['abs_speed_decreased'] for x in edges])),
        incident_velocity_increased_fraction=float(np.mean([x['abs_speed_increased'] for x in edges])),
        incident_velocity_sign_reversed_fraction=float(np.mean([x['sign_reversed'] for x in edges])),
        ev_sign_change_fraction=float(np.mean([x['e_v_sign_before']!=x['e_v_sign_after'] for x in edges])),
        re_trigger=followup_statistics(rows,end),
        same_satellite_HIGH_start_interval_s=distribution(x['next_HIGH_start_interval_s'] for x in followed),
        next_HIGH_at_LOW_dwell_release_count=sum(abs(x['next_wait_s']-c.dwell_time_s)<1e-8 for x in followed),
        same_HIGH_interval_at_minimum_count=sum(abs(x['next_HIGH_start_interval_s']-2*(c.slew_time_s+c.dwell_time_s))<1e-8 for x in followed),
        LOW_release_already_above_on_fraction=float(np.mean([x['q_at_LOW_dwell_release']>=c.threshold_on_m_s2 for x in followed])),
        neighbor_D_exceeds_natural_fraction=float(np.mean([abs(x['next_qD_neighbor'])>abs(x['next_qD_natural']) for x in followed])),
        neighbor_D_change_abs_m_s2=distribution(abs(x['next_qD_neighbor']) for x in followed),
        natural_D_change_abs_m_s2=distribution(abs(x['next_qD_natural']) for x in followed),
        dwell_limited_returns=residual['summary']['dwell_limited_returns'],
        threshold_limited_returns=residual['summary']['threshold_limited_returns'],
        quiet=dict(total_intervals=len(quiet),observed_next_HIGH_intervals=len(completedquiet),right_censored_intervals=sum(x['right_censored'] for x in quiet),
            duration_observed_s=distribution(x['actual_next_wait_s'] for x in completedquiet),
            total_observed_allLOW_time_s=float(sum(x['observation_window_end_s']-x['completion_time_s'] for x in quiet)),
            right_censored_duration_lower_bound_s=[end-x['completion_time_s'] for x in quiet if x['right_censored']],
            matched_valid_forecasts=len(matched),matched_forecast_abs_error_s=distribution(abs(x['validated_prediction_error_s']) for x in matched)),
        closure=dict(pulse_max_m_s=max(closures),local_LOW_followup_max_qD_m_s2=max(qclosures,default=0)),
        cycle_counting_note='completed LOW return <= tmax; pending cycles excluded from residual distributions; adjacent-gap samples overlap')
    write_json(out/f'{label}_metrics.json',clean_json(stat))
    write_csv(out/f'{label}_satellites.csv',[dict(satellite_id=j+1,cycles=metric['cycles_per_satellite'][j],
        transitions=metric['switches_per_satellite'][j],HIGH_equivalent_s=metric['high_equivalent_s_per_satellite'][j],
        altitude_loss_m=float(a[0,j]-a[-1,j]),extra_altitude_loss_vs_same_allLOW_m=float(low_extra[-1,j]),
        final_altitude_m=metric['final_altitude_m'][j]) for j in range(n)])
    write_csv(out/f'{label}_events.csv',events)
    saveix=np.arange(0,len(t),6) # 60s plot table only, metrics always original 10s.
    np.savez_compressed(out/f'{label}_plot_series.npz',time_s=t[saveix],error_m=e[saveix])
    write_csv(out/f'{label}_history_60s.csv',[dict(time_s=float(t[k]),max_error_m=float(np.max(abs(e[k]))),
        instantaneous_RMS_m=float(np.sqrt(np.mean(e[k]**2))),max_gap_rate_m_s=float(np.max(abs(v[k]))),
        stable_HIGH_count=int(np.sum(arrays['mode'][k]==2)),nonLOW_count=int(np.sum(arrays['mode'][k]!=0)),
        extra_mean_altitude_loss_m=float(np.mean(low_extra[k]))) for k in saveix])
    return stat,rows,edges


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=DEFAULT)
    args=ap.parse_args();out=args.root/'analysis';out.mkdir(parents=True,exist_ok=False)
    before=snapshot();write_json(out/'input_manifest.json',before)
    base=load_config(SOURCE/'normalized_config.yaml');raw=base.to_dict();raw['simulation']['duration_s']=7*86400.;cfg300=config_from_dict(raw)
    cfg150=load_config(args.root/'D150_7d/normalized_config.yaml')
    old=cfg300.to_dict();new=cfg150.to_dict();new['controller']['dwell_time_s']=old['controller']['dwell_time_s'];assert old==new
    p1=json.loads((SOURCE/'metadata.json').read_text())['profiles'][0]
    p2=json.loads((args.root/'D150_7d/metadata.json').read_text())['profiles'][0];assert p1==p2
    source_env=json.loads((SOURCE/'environment.json').read_text());new_env=json.loads((args.root/'D150_7d/environment.json').read_text())
    for key in ['environment','weather_sha256','frozen_profiles_sha256','cache_seconds','pymsis_version','msis_version']:
        assert source_env[key]==new_env[key],key
    low=load_arrays(SOURCE,3,cfg300)
    data={label:load_arrays(folder,member,cfg) for label,folder,member,cfg in
        [('D300',SOURCE,0,cfg300),('D150',args.root/'D150_7d',0,cfg150)]}
    for key in ['a_m','lambda_rad','density_kg_m3','area_m2','mode']:
        np.testing.assert_array_equal(data['D300'][key][0],data['D150'][key][0])
    common={}
    for key in ['a_m','lambda_rad','density_kg_m3','area_m2','gap_error_m','gap_rate_m_s']:
        x=data['D300'][key][:28];y=data['D150'][key][:28] # through t=270s inclusive
        common[key]=float(np.max(abs(x-y)))
        assert np.allclose(x,y,rtol=0,atol=1e-10 if key!='density_kg_m3' else 1e-22),key
    stats={};cycles={}
    for label,cfg,folder,member in [('D300',cfg300,SOURCE,0),('D150',cfg150,args.root/'D150_7d',0)]:
        print('Analyze',label,flush=True)
        events=events_for(folder,member,7*86400)
        stats[label],cycles[label],_=diagnose(data[label],cfg,events,low,out,label)
    paired=[]
    for row in cycles['D300']:
        if row['start_s']!=0:continue
        other=next(x for x in cycles['D150'] if x['start_s']==0 and x['satellite_id']==row['satellite_id'])
        paired.append(dict(satellite_id=row['satellite_id'],D300_end_s=row['end_s'],D150_end_s=other['end_s'],
            D300_pulse_m_s=row['pulse_g_impulse_m_s'],D150_pulse_m_s=other['pulse_g_impulse_m_s'],
            D150_D300_pulse_ratio=other['pulse_g_impulse_m_s']/row['pulse_g_impulse_m_s'],
            D300_residual_max_m_s=row['V_incident_after_m_s'],D150_residual_max_m_s=other['V_incident_after_m_s']))
    write_csv(out/'paired_initial_pulses.csv',paired)
    ref=ROOT/'tests/fixtures/regression_24h';suite=ROOT/'outputs/reference_dispersion_20261012/zero_regression_24h/suite'
    regression=compare_baseline(suite,ref);assert regression['passed'];regression['execution']='recompared saved previously regenerated suite; no new baseline simulation'
    write_json(out/'regression_24h.json',regression)
    assert before==snapshot()
    report=dict(cases=stats,input_check=dict(passed=True,changes=['controller.dwell_time_s: 300 -> 150'],
        common_0_to_270s_max_abs_difference=common,profile_sha256=new_env['frozen_profiles_sha256'],
        weather_sha256=new_env['weather_sha256'],old_inputs_unchanged=len(before)),
        matched_initial_pulse_count=len(paired),
        aggregate_pulse_median_ratio=stats['D150']['pulse_velocity_contribution_m_s']['median']/stats['D300']['pulse_velocity_contribution_m_s']['median'],
        theory_minimum_equivalent_ratio=270/420,
        interpretation='Total effect of simultaneous HIGH and LOW dwell intervention; pulse-size and refractory-gate pathways are not independently identified.')
    write_json(out/'comparison.json',clean_json(report))
    print(json.dumps(clean_json(report),ensure_ascii=False),flush=True)


if __name__=='__main__':main()
