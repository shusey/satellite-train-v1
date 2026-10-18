"""One common postprocessor for the four 30-day factorial cells; no integration."""
import argparse
import csv
import gc
import hashlib
import json
from pathlib import Path
import sys
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'));sys.path.insert(0,str(ROOT/'src'))
from run_reference_dwell_factorial import OUTPUT,SOURCE,PROFILE
from analyze_reference_a_dwell import diagnose,load_arrays,events_for
from analyze_reference_pulses import write_csv
from sat_string.config import load_config,config_from_dict,dump_normalized_config
from sat_string.dwell_horizon import joined_saved_parts
from sat_string.dispersion_control import derived_histories
from sat_string.dwell_intervention_metrics import distribution
from sat_string.dwell_factorial import CASE_ORDER,contrasts
from sat_string.long_horizon import decompose_q
from sat_string.control_cost import clean_json
from sat_string.io import write_json

D150_PARTS=[ROOT/'outputs/reference_a_dwell_20261017/D150_7d',
            ROOT/'outputs/reference_a_dwell_60d_20261017/D150_30d']
KEYS=['a_m','lambda_rad','area_m2','mode','density_kg_m3']


def read_csv(path):
    with path.open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))


def equal_environment(left,right):
    keys=['environment','weather_sha256','frozen_profiles_sha256','cache_seconds','pymsis_version','msis_version']
    if any(left[k]!=right[k] for k in keys):raise ValueError('Environment changed')


def analyze(label):
    out=OUTPUT/'analysis'/label;out.mkdir(parents=True,exist_ok=False)
    base=load_config(SOURCE/'normalized_config.yaml')
    folder=SOURCE if label=='H300-L300' else D150_PARTS[-1] if label=='H150-L150' else OUTPUT/label
    cfg=load_config(folder/'normalized_config.yaml')
    raw=cfg.to_dict();raw['simulation']['duration_s']=30*86400.;cfg=config_from_dict(raw)
    b,c=base.to_dict(),cfg.to_dict()
    for x in [b,c]:
        x['simulation'].pop('duration_s')
        for k in ['dwell_time_s','dwell_high_time_s','dwell_low_time_s']:x['controller'].pop(k,None)
    assert b==c
    high,low=cfg.controller.effective_dwell_high_s,cfg.controller.effective_dwell_low_s
    assert label==f'H{high:g}-L{low:g}'
    baseline_meta=json.loads((SOURCE/'metadata.json').read_text())
    meta=json.loads((folder/'metadata.json').read_text())
    assert meta['profiles'][0]==baseline_meta['profiles'][0]
    equal_environment(json.loads((SOURCE/'environment.json').read_text()),json.loads((folder/'environment.json').read_text()))
    off=load_arrays(SOURCE,3,base,days=30)
    checks=[]
    if label=='H150-L150':
        for part in D150_PARTS:
            equal_environment(json.loads((SOURCE/'environment.json').read_text()),json.loads((part/'environment.json').read_text()))
        t,x,checks=joined_saved_parts(D150_PARTS,0,KEYS)
        arrays=derived_histories(t,x['a_m'],x['lambda_rad'],x['area_m2'],x['mode'],x['density_kg_m3'],cfg)
        del x
    else:arrays=load_arrays(folder,0,cfg,days=30)
    for key in KEYS:
        np.testing.assert_array_equal(arrays[key][0],load_arrays_initial(SOURCE,key))
    events=events_for(folder,0,30*86400.)
    print('Analyze '+label,flush=True)
    stat,cycles,edges=diagnose(arrays,cfg,events,off,out,label)
    t=arrays['time_s'];end=float(t[-1]);qP,qD=decompose_q(arrays['gap_error_m'],arrays['gap_rate_m_s'],
        cfg.derived.controller_kp_s2_inv,cfg.derived.controller_kd_s_inv)
    quiet=read_csv(out/f'{label}_quiet_intervals.csv')
    actual_rows=[]
    for row in quiet:
        start=float(row['completion_time_s']);stop=float(row['observation_window_end_s'])
        j=int(row['actual_next_satellite_id'])-1 if row['actual_next_satellite_id'] else None
        k0=int(np.searchsorted(t,start));k1=int(np.searchsorted(t,stop))
        actual_rows.append(dict(start_s=start,end_s=stop,duration_s=stop-start,
            right_censored=row['right_censored']=='True',next_satellite_id=None if j is None else j+1,
            start_qP=None if j is None else qP[k0,j],start_qD=None if j is None else qD[k0,j],
            trigger_qP=None if j is None else qP[k1,j],trigger_qD=None if j is None else qD[k1,j],
            predicted_satellite_id=row['predicted_satellite_id'],
            forecast_valid=row['prediction_valid_on_observed_window']=='True',
            target_matches=row['prediction_target_matches_actual']=='True',
            timing_error_s=float(row['prediction_error_actual_minus_predicted_s']) if row['prediction_error_actual_minus_predicted_s'] else None))
    write_csv(out/'quiet_diagnostics.csv',actual_rows)
    durations=[x['duration_s'] for x in actual_rows]
    all_low_fraction=float(np.sum(np.diff(t)*np.all(arrays['mode'][:-1]==0,axis=1))/end)
    # Initial all-LOW time, if any, is separate from response-end forecasts.
    first_high=min((e['start_time'] for e in events if e['to_mode']=='H'),default=end)
    assert np.isclose((sum(durations)+first_high)/end,all_low_fraction,atol=1e-12,rtol=0)
    valid=[r for r in actual_rows if r['forecast_valid'] and not r['right_censored']]
    matched=[r for r in valid if r['target_matches']]
    stat['all_LOW']=dict(count=len(actual_rows)+(first_high>0),fraction=all_low_fraction,
        initial_interval_s=first_high,duration_with_censored_lower_bound_s=distribution(durations+([first_high] if first_high else [])),
        observed_duration_s=distribution(r['duration_s'] for r in actual_rows if not r['right_censored']),
        right_censored=sum(r['right_censored'] for r in actual_rows),
        valid_observed_forecasts=len(valid),target_matches=len(matched),
        matched_abs_timing_error_s=distribution(abs(r['timing_error_s']) for r in matched),
        next_satellite_counts={str(j):sum(r['next_satellite_id']==j for r in actual_rows) for j in range(1,22)})
    followed=[r for r in cycles if r['next_HIGH_s'] is not None]
    qualifying=[r for r in followed if abs(r['next_qD_neighbor'])>abs(r['next_qD_natural'])]
    stat['neighbor_qualifying_intervals']=len(qualifying)
    stat['neighbor_qualifying_intervals_per_day']=len(qualifying)/30
    stat['theoretical_minimum_LOW_complete_to_HIGH_s']=low
    stat['re_trigger']['observed_wait_count']=len(followed)
    for seg in stat['formation']['segments']:
        lo,hi=seg['start_day']*86400,seg['end_day']*86400
        count=sum(lo<=r['end_s'] and (r['end_s']<=hi if hi==end else r['end_s']<hi) for r in cycles)
        seg['completed_cycles']=int(count);seg['cycles_per_day']=count/((hi-lo)/86400)
        seg['HIGH_equivalent_satellite_hours_per_day']=seg['high_equivalent_satellite_s_per_day']/3600
    pulse_sum=sum(r['high_equivalent_s'] for r in cycles)
    assert stat['formation']['completed_cycles']==len(cycles)
    assert min(r['high_equivalent_s'] for r in cycles)==stat['theoretical_minimum_HIGH_equivalent_s']
    # Allow incomplete final pulses explicitly; event integral remains authoritative.
    stat['completed_cycle_HIGH_equivalent_sum_s']=pulse_sum
    stat['incomplete_pulse_HIGH_equivalent_s']=stat['formation']['high_equivalent_total_s']-pulse_sum
    assert stat['incomplete_pulse_HIGH_equivalent_s']>=-1e-8
    assert stat['closure']['pulse_max_m_s']<1e-5
    assert stat['closure']['local_LOW_followup_max_qD_m_s2']<1e-9
    stat['verification']=dict(inputs_match_except_resolved_dwell=True,profile_exact=True,environment_exact=True,
        continuation_boundaries=checks,profile_sha256=hashlib.sha256(PROFILE.read_bytes()).hexdigest(),
        no_order_reversal=stat['formation']['first_label_order_reversal_s'] is None)
    write_json(out/f'{label}_metrics.json',clean_json(stat))
    dump_normalized_config(cfg,out/'evaluation_config.yaml')
    del arrays,off,qP,qD;gc.collect()
    print(json.dumps({k:stat['formation'][k] for k in ['max_error_m','rms_error_m','completed_cycles','extra_mean_altitude_loss_m']}),flush=True)


def load_arrays_initial(folder,key):
    return np.load(folder/f'{key}.npy',mmap_mode='r')[0,0]


def summarize():
    stats={label:json.loads((OUTPUT/'analysis'/label/f'{label}_metrics.json').read_text()) for label in CASE_ORDER}
    metrics={
        'max_error_m':lambda s:s['formation']['max_error_m'],
        'RMS_m':lambda s:s['formation']['rms_error_m'],
        'cycles':lambda s:s['completed_cycles'],
        'HIGH_equivalent_satellite_hours':lambda s:s['formation']['high_equivalent_total_s']/3600,
        'extra_altitude_loss_m':lambda s:s['formation']['extra_mean_altitude_loss_m'],
        'residual_incident_abs_velocity_median_mm_s':lambda s:s['residual_abs_incident_gap_velocity_m_s']['median']*1000,
        'residual_cycle_max_velocity_median_mm_s':lambda s:s['residual_max_incident_velocity_per_cycle_m_s']['median']*1000}
    rows=[]
    for label,s in stats.items():
        m=s['formation'];rows.append(dict(case=label,**{k:fn(s) for k,fn in metrics.items()},
            final_max_error_m=m['final_max_error_m'],transitions=m['transitions'],
            exceed_100m_fraction=m['exceedances'][0]['exceedance_fraction'],exceed_500m=m['exceedances'][1]['ever_exceeded'],
            allLOW_time_fraction=s['all_LOW']['fraction'],
            theoretical_min_HIGH_equivalent_s=s['theoretical_minimum_HIGH_equivalent_s'],
            theoretical_min_HIGH_restart_s=s['theoretical_minimum_HIGH_restart_s'],
            observed_retrigger_median_h=s['re_trigger']['observed_wait_s']['median']/3600))
    contrasts_rows=[dict(metric=k,**contrasts([fn(stats[label]) for label in CASE_ORDER])) for k,fn in metrics.items()]
    direct=[]
    for factor,pairs in [('HIGH',[(0,1),(2,3)]),('LOW',[(0,2),(1,3)])]:
        for i,j in pairs:
            for name,fn in metrics.items():
                before,after=fn(stats[CASE_ORDER[i]]),fn(stats[CASE_ORDER[j]])
                direct.append(dict(factor=factor,baseline=CASE_ORDER[i],intervention=CASE_ORDER[j],metric=name,
                    before=before,after=after,difference=after-before,change_percent=100*(after/before-1) if before else None))
    write_csv(OUTPUT/'analysis/comparison.csv',rows)
    write_csv(OUTPUT/'analysis/contrasts.csv',contrasts_rows)
    write_csv(OUTPUT/'analysis/direct_comparisons.csv',direct)
    write_json(OUTPUT/'analysis/comparison.json',dict(case_order=list(CASE_ORDER),cases=stats,
        contrasts=contrasts_rows,direct_comparisons=direct,
        inference='Deterministic single-profile descriptive contrasts; no statistical ANOVA or permanent stability claim'))
    print(json.dumps(rows,indent=2),flush=True)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--case',choices=CASE_ORDER);ap.add_argument('--summarize',action='store_true')
    args=ap.parse_args()
    if args.summarize:summarize()
    elif args.case:analyze(args.case)
    else:ap.error('Choose one case or --summarize')
