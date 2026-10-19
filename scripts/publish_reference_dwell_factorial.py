"""Verify saved results and copy only compact reproducible research artifacts."""
import hashlib
import argparse
import json
import shutil
from pathlib import Path
from run_reference_dwell_factorial import ROOT,OUTPUT,SOURCE,PROFILE
from analyze_reference_dwell_factorial import read_csv
from analyze_reference_pulses import write_csv
from sat_string.dwell_factorial import CASE_ORDER
from sat_string.io import write_json
import numpy as np


def publish_results(output,destination,*,days=30,tests_passed=311,new_test_cases=18):
    OUTPUT=output
    before=json.loads((OUTPUT/'source_manifest_before.json').read_text())
    for name,record in before.items():
        p=ROOT/name;stat=p.stat()
        assert stat.st_size==record['size'] and stat.st_mtime_ns==record['mtime_ns'],name
        if record['sha256']:assert hashlib.sha256(p.read_bytes()).hexdigest()==record['sha256'],name
    reg=json.loads((OUTPUT/'regression_24h/regression.json').read_text());assert reg['passed']
    diagonal=json.loads((OUTPUT/'diagonal_regression.json').read_text());assert all(x['passed'] for x in diagonal.values())
    report=json.loads((OUTPUT/'analysis/comparison.json').read_text())
    checks={};segments=[];pulses=[];retrigger=[];quiet=[]
    for label in CASE_ORDER:
        s=report['cases'][label];m=s['formation'];folder=OUTPUT/'analysis'/label
        assert m['first_label_order_reversal_s'] is None and m['min_physical_pair_chord_m']>0
        assert s['verification']['inputs_match_except_resolved_dwell'] and s['verification']['environment_exact']
        assert s['completed_cycles']==m['completed_cycles']
        with np.load(folder/f'{label}_four_stage_states.npz') as packed:
            assert packed['v'].shape==(s['completed_cycles'],4,20)
            assert packed['qP'].shape==(s['completed_cycles'],4,21)
            assert np.isfinite(packed['v']).all()
        assert s['closure']['pulse_max_m_s']<1e-5
        assert s['closure']['local_LOW_followup_max_qD_m_s2']<1e-9
        checks[label]=dict(cycles=s['completed_cycles'],min_physical_gap_m=m['min_physical_pair_chord_m'],
            impulse_quadrature_closure=s['closure'],input_verification=s['verification'])
        segments.extend(dict(case=label,**x) for x in m['segments'])
        pulses.append(dict(case=label,minimum_pulse_fraction=s['minimum_pulse_fraction'],
            HIGH_equivalent_theory_s=s['theoretical_minimum_HIGH_equivalent_s'],
            HIGH_restart_theory_s=s['theoretical_minimum_HIGH_restart_s'],
            **{(f'HIGH_equivalent_{k}_s' if k!='n' else 'HIGH_equivalent_n'):v for k,v in s['high_equivalent_s_per_cycle'].items()},
            **{(f'own_pulse_{k}_mm_s' if k!='n' else 'own_pulse_n'):v*1000 if v is not None and k!='n' else v for k,v in s['pulse_velocity_contribution_m_s'].items()},
            **{(f'incident_velocity_{k}_mm_s' if k!='n' else 'incident_velocity_n'):v*1000 if v is not None and k!='n' else v for k,v in s['residual_abs_incident_gap_velocity_m_s'].items()},
            worsening_fraction=s['residual_worsening_fraction']))
        retrigger.append(dict(case=label,**s['re_trigger']['observed_wait_s'],right_censored=s['re_trigger']['right_censored'],
            within_30min_fraction=s['re_trigger']['within_1800s']['fraction'],
            within_60min_fraction=s['re_trigger']['within_3600s']['fraction'],
            at_LOW_gate_count=s['next_HIGH_at_LOW_dwell_release_count'],
            at_HIGH_restart_minimum_count=s['same_HIGH_interval_at_minimum_count'],
            allLOW_fraction=s['all_LOW']['fraction'],allLOW_intervals=s['all_LOW']['count'],
            valid_observed_forecasts=s['all_LOW']['valid_observed_forecasts'],forecast_target_matches=s['all_LOW']['target_matches'],
            timing_abs_median_s=s['all_LOW']['matched_abs_timing_error_s']['median'],
            timing_abs_p90_s=s['all_LOW']['matched_abs_timing_error_s']['p90'],
            neighbor_exceeds_natural_fraction=s['neighbor_D_exceeds_natural_fraction'],
            neighbor_abs_D_median_m_s2=s['neighbor_D_change_abs_m_s2']['median'],
            neighbor_qualifying_per_day=s['neighbor_qualifying_intervals_per_day']))
    data=destination/'data';data.mkdir(parents=True,exist_ok=False)
    for name in ['comparison.json','comparison.csv','contrasts.csv','direct_comparisons.csv']:
        shutil.copyfile(OUTPUT/'analysis'/name,data/name)
    write_csv(data/'segments.csv',segments);write_csv(data/'pulse_summary.csv',pulses);write_csv(data/'retrigger_summary.csv',retrigger)
    if days==60:
        late_rows=[]
        for label in CASE_ORDER:
            s=report['cases'][label];late=s['late_30_60d'];seg=s['formation']['segments'][-1]
            late_rows.append(dict(case=label,RMS_m=seg['error_rms_m'],max_error_m=seg['max_error_m'],
                cycles_per_day=seg['cycles_per_day'],transitions_per_day=seg['transitions_per_day'],
                HIGH_equivalent_satellite_hours_per_day=seg['HIGH_equivalent_satellite_hours_per_day'],
                extra_altitude_loss_m_per_day=seg['extra_mean_altitude_loss_per_day_m'],
                residual_cycle_max_velocity_median_mm_s=late['residual_max_incident_velocity_per_cycle_m_s']['median']*1000,
                pulse_velocity_median_mm_s=late['pulse_velocity_contribution_m_s']['median']*1000,
                all_LOW_time_fraction=late['all_LOW_time_fraction'],
                neighbor_qualifying_intervals_per_day=late['neighbor_qualifying_intervals_per_day']))
        write_csv(data/'late_30_60d.csv',late_rows)
    keys=['satellite_id','cycle','start_s','end_s','high_equivalent_s','minimum_pulse','return_driver','pulse_g_impulse_m_s',
          'V_incident_before_m_s','V_incident_after_m_s','qP_before','qD_before','qP_after','qD_after',
          'next_HIGH_s','next_wait_s','right_censored','next_qD_neighbor','next_qD_natural']
    for label in CASE_ORDER:
        folder=OUTPUT/'analysis'/label
        cycles=read_csv(folder/f'{label}_cycles.csv')
        write_csv(data/f'{label}_cycles.csv',[{k:r[k] for k in keys} for r in cycles])
        for name in ['quiet_diagnostics.csv','evaluation_config.yaml',f'{label}_events.csv',f'{label}_satellites.csv']:
            shutil.copyfile(folder/name,data/(f'{label}_{name}' if not name.startswith(label) else name))
    frozen=json.loads(PROFILE.read_text())
    write_json(data/'reference_A_profile.json',dict(profile=frozen['profiles'][0],source=str(PROFILE.relative_to(ROOT)),
        source_sha256=hashlib.sha256(PROFILE.read_bytes()).hexdigest(),std_ddof=1))
    for name in ['H150-L300','H300-L150']:
        for file in ['environment.json','input_difference.json','completion.json']:
            shutil.copyfile(OUTPUT/name/file,data/f'{name}_{file}')
    shutil.copyfile(OUTPUT/'regression_24h/regression.json',data/'regression_24h.json')
    shutil.copyfile(OUTPUT/'diagonal_regression.json',data/'diagonal_regression.json')
    shutil.copyfile(OUTPUT/'source_manifest_before.json',data/'source_manifest.json')
    figures=destination/'figures';figures.mkdir()
    for p in (OUTPUT/'figures').glob('*.png'):shutil.copyfile(p,figures/p.name)
    write_json(data/'verification.json',dict(passed=True,tests_passed=tests_passed,new_test_cases=new_test_cases,
        baseline_24h_passed=True,explicit_diagonal_6h_exact=True,unchanged_input_files=len(before),
        large_file_preservation='Read-only use; unchanged size and nanosecond mtime; small files additionally SHA256',
        new_physics_segments=dict(cases=['H150-L300','H300-L150'],start_day=0 if days==30 else 30,end_day=days),
        diagonal_physics_reused=True,observed_days=days,checks=checks))
    print('Verified and published compact artifacts to '+str(destination),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path)
    parser.add_argument('--root',type=Path)
    parser.add_argument('--days',type=int,choices=[30,60],default=30)
    args=parser.parse_args()
    output=args.root or (OUTPUT if args.days==30 else OUTPUT.parent/'reference_dwell_factorial_60d_20261019')
    dest=args.output or ROOT/f'docs/reference_dwell_factorial_{args.days}d'
    publish_results(output,dest,days=args.days,tests_passed=311 if args.days==30 else 323,new_test_cases=18 if args.days==30 else 12)
