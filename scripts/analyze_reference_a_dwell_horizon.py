"""Extend only A-ref D150; reuse archived D300 diagnostics and trajectories."""
import argparse,csv,json,gc,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'));sys.path.insert(0,str(ROOT/'src'))
from sat_string.config import load_config,config_from_dict
from sat_string.dispersion_control import derived_histories
from sat_string.controller import gap_kinematics
from sat_string.dwell_horizon import joined_saved_parts
from sat_string.dwell_intervention_metrics import distribution,followup_statistics
from sat_string.regression import compare_baseline
from sat_string.io import write_json
from sat_string.control_cost import clean_json
from analyze_reference_a_dwell import diagnose
from analyze_reference_pulses import read_csv,write_csv

ORIGINAL=ROOT/'outputs/reference_dispersion_20261012'
SEVEN=ROOT/'outputs/reference_a_dwell_20261017'
LONG=ROOT/'outputs/reference_a_dwell_60d_20261017'
PULSE=ROOT/'outputs/reference_pulse_mechanism_20261017_v3'


def archive_D300(days,out):
    """Reuse old metrics/whole-cycle tables; censor using the NEW observation end."""
    end=days*86400.;folder=ROOT/f'docs/reference_dispersion/data/{days}d'
    metric=json.loads((folder/'A_metrics.json').read_text())
    rows=[]
    for x in read_csv(PULSE/'A_pulse_cycles.csv'):
        if float(x['end_s'])>end:continue
        nxt=float(x['next_own_HIGH_s']) if x['next_own_HIGH_s'] else None
        if nxt is not None and nxt>end:nxt=None
        rows.append(dict(satellite_id=int(x['satellite_id']),cycle=int(x['cycle']),start_s=float(x['start_s']),end_s=float(x['end_s']),
            next_HIGH_s=nxt,next_wait_s=nxt-float(x['end_s']) if nxt is not None else None,
            next_HIGH_start_interval_s=nxt-float(x['start_s']) if nxt is not None else None,
            pulse_g_impulse_m_s=float(x['own_excess_g_impulse_m_s']),
            V_incident_after_m_s=float(x['V_incident_after_m_s']),high_equivalent_s=float(x['high_equivalent_s']),
            minimum_pulse=x['minimum_pulse']=='True'))
    edges=[x for x in read_csv(PULSE/'A_incident_gaps.csv') if float(x['end_s'])<=end]
    quiet=read_csv(folder/'A_quiet_cycles.csv');uncensored=[x for x in quiet if x['actual_next_command_s']]
    valid=[x for x in uncensored if x['prediction_valid_on_observed_window']=='True' and x['prediction_target_matches_actual']=='True']
    stat=dict(label='D300',dwell_s=300,days=days,formation=metric,completed_cycles=len(rows),
        minimum_pulse_count=sum(x['minimum_pulse'] for x in rows),minimum_pulse_fraction=float(np.mean([x['minimum_pulse'] for x in rows])),
        pulse_velocity_contribution_m_s=distribution(x['pulse_g_impulse_m_s'] for x in rows),
        high_equivalent_s_per_cycle=distribution(x['high_equivalent_s'] for x in rows),
        residual_abs_incident_gap_velocity_m_s=distribution(abs(float(x['v_after_m_s'])) for x in edges),
        residual_max_incident_velocity_per_cycle_m_s=distribution(x['V_incident_after_m_s'] for x in rows),
        re_trigger=followup_statistics(rows,end),
        same_satellite_HIGH_start_interval_s=distribution(x['next_HIGH_start_interval_s'] for x in rows if x['next_HIGH_s'] is not None),
        quiet=dict(total_intervals=len(quiet),observed_next_HIGH_intervals=len(uncensored),
            right_censored_intervals=len(quiet)-len(uncensored),
            duration_observed_s=distribution(float(x['actual_next_wait_s']) for x in uncensored),
            total_observed_allLOW_time_s=sum(float(x['observation_window_end_s'])-float(x['completion_time_s']) for x in quiet),
            matched_valid_forecasts=len(valid),matched_forecast_abs_error_s=distribution(abs(float(x['validated_prediction_error_s'])) for x in valid)),
        reuse_provenance=[str(folder.relative_to(ROOT)),str(PULSE.relative_to(ROOT))])
    assert len(rows)==metric['completed_cycles']
    write_json(out/'D300_metrics.json',stat);write_csv(out/'D300_cycles.csv',rows)
    write_csv(out/'D300_quiet_intervals.csv',quiet)
    checkpoint=json.loads((ORIGINAL/f'batch_{days}d/checkpoint.json').read_text())
    write_csv(out/'D300_events.csv',[x for x in checkpoint['events'][0] if x['command_time']<=end])
    cache=np.load(ORIGINAL/f'A_plot_{days}d.npz')
    t,error,on_a,off_a=[cache[k] for k in ['time_s','on_error','a_m','off_a_m']]
    ix=np.arange(0,len(t),6)
    write_csv(out/'D300_history_60s.csv',[dict(time_s=float(t[k]),max_error_m=float(abs(error[k]).max()),
        instantaneous_RMS_m=float(np.sqrt(np.mean(error[k]**2))),
        extra_mean_altitude_loss_m=float(np.mean(off_a[k]-on_a[k]))) for k in ix])
    return stat


def preservation_manifest():
    paths=[]
    for folder in [ORIGINAL/'batch_30d',ORIGINAL/'batch_60d',SEVEN/'D150_7d']:
        paths.extend(p for p in folder.iterdir() if p.is_file())
    result={}
    import hashlib
    for p in paths:
        s=p.stat();result[str(p.relative_to(ROOT))]=dict(size=s.st_size,mtime_ns=s.st_mtime_ns)
        if s.st_size<20_000_000:result[str(p.relative_to(ROOT))]['sha256']=hashlib.sha256(p.read_bytes()).hexdigest()
    return result


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--days',type=int,choices=[30,60],required=True);args=ap.parse_args()
    days=args.days;out=LONG/f'analysis_{days}d';out.mkdir(parents=True,exist_ok=False)
    before=preservation_manifest();write_json(out/'input_manifest.json',before)
    paths=[SEVEN/'D150_7d',LONG/'D150_30d']+([LONG/'D150_60d'] if days==60 else [])
    for p in paths:
        if not json.loads((p/'completion.json').read_text())['complete']:raise ValueError('Incomplete part')
    cfg=load_config(paths[-1]/'normalized_config.yaml')
    assert cfg.controller.dwell_time_s==150 and cfg.simulation.duration_s==days*86400.
    events=json.loads((paths[-1]/'checkpoint.json').read_text())['events'][0]
    old_events=json.loads((paths[0]/'checkpoint.json').read_text())['events'][0]
    assert events[:len(old_events)]==old_events
    # Only a/lambda are needed for the reused all-LOW accuracy/altitude reference.
    offpaths=[ORIGINAL/'batch_30d']+([ORIGINAL/'batch_60d'] if days==60 else [])
    t,off,offchecks=joined_saved_parts(offpaths,3,['a_m','lambda_rad'])
    e,v,gap=gap_kinematics(off['a_m'],off['lambda_rad'],reference_radius_m=cfg.derived.reference_radius_m,
        target_spacing_rad=cfg.derived.target_spacing_rad,mu_m3_s2=cfg.constants.mu_m3_s2)
    low=dict(time_s=t,a_m=off['a_m'],lambda_rad=off['lambda_rad'],gap_error_m=e,gap_distance_m=gap)
    del v,off
    t,raw,checks=joined_saved_parts(paths,0,['a_m','lambda_rad','area_m2','mode','density_kg_m3'])
    assert np.array_equal(t,low['time_s'])
    arrays=derived_histories(t,raw['a_m'],raw['lambda_rad'],raw['area_m2'],raw['mode'],raw['density_kg_m3'],cfg)
    print(f'Analyze D150 at {days}d',flush=True)
    stat,rows,_=diagnose(arrays,cfg,events,low,out,'D150')
    # Interval cycles are assigned at LOW completion, avoiding start/end ambiguity.
    for interval in stat['formation']['segments']:
        left=interval['start_day']*86400;right=interval['end_day']*86400
        count=int(sum(left<=x['end_s'] and (x['end_s']<=right if right==t[-1] else x['end_s']<right) for x in rows))
        interval['completed_cycles']=count;interval['cycles_per_day']=count/((right-left)/86400)
        interval['HIGH_equivalent_satellite_hours_per_day']=interval['high_equivalent_satellite_s_per_day']/3600
    stat=clean_json(stat);write_json(out/'D150_metrics.json',stat)
    del low,arrays,raw,e,gap;gc.collect()
    base=archive_D300(days,out)
    reg=compare_baseline(ORIGINAL/'zero_regression_24h/suite',ROOT/'tests/fixtures/regression_24h')
    assert reg['passed'];reg['execution']='recompare saved previously regenerated 24h suite; no new baseline simulation'
    write_json(out/'regression_24h.json',reg)
    assert preservation_manifest()==before
    report=dict(days=days,cases=dict(D300=base,D150=stat),verification=dict(passed=True,
        old_input_files_unchanged=len(before),continuation_boundaries=checks+offchecks,
        seven_day_prefix_events_exact=True,regression_24h_passed=True),
        provenance=dict(D150_parts=[str(p.relative_to(ROOT)) for p in paths],D300='reused saved metrics, pulse diagnostics and plot histories'))
    write_json(out/'comparison.json',clean_json(report))
    key=stat['formation'];print(json.dumps(dict(days=days,max_error_m=key['max_error_m'],RMS_m=key['rms_error_m'],
        cycles=stat['completed_cycles'],extra_height_m=key['extra_mean_altitude_loss_m'],min_gap_m=key['min_signed_reference_gap_m']),ensure_ascii=False),flush=True)


if __name__=='__main__':main()
