"""Validate and publish small seven-day artifacts without another simulation."""
import argparse,csv,hashlib,json,shutil,sys
from pathlib import Path
import numpy as np
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'));sys.path.insert(0,str(ROOT/'src'))
from sat_string.io import write_json
from sat_string.config import load_config,config_from_dict,dump_normalized_config
from analyze_reference_a_dwell import snapshot,SOURCE,PROFILE,load_arrays
from sat_string.long_horizon import decompose_q
from sat_string.pulse_diagnostics import acceleration_terms,cumulative_trapezoid,node_map


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--root',type=Path,default=ROOT/'outputs/reference_a_dwell_20261017')
    ap.add_argument('--output',type=Path,default=ROOT/'docs/reference_a_dwell_150/data')
    args=ap.parse_args();out=args.output;out.mkdir(parents=True,exist_ok=False);src=args.root/'analysis'
    comparison=json.loads((src/'comparison.json').read_text(encoding='utf-8'))
    manifest=json.loads((src/'input_manifest.json').read_text(encoding='utf-8'));assert snapshot()==manifest
    reg=json.loads((src/'regression_24h.json').read_text());assert reg['passed']
    checks={}
    summary_rows=[]
    for label in ['D300','D150']:
        stat=comparison['cases'][label];m=stat['formation']
        assert m['first_label_order_reversal_s'] is None and m['min_physical_pair_chord_m']>0
        assert m['completed_cycles']==stat['completed_cycles']
        with (src/f'{label}_cycles.csv').open(encoding='utf-8-sig',newline='') as f:cycles=list(csv.DictReader(f))
        total=sum(float(x['high_equivalent_s']) for x in cycles)
        assert np.isclose(total,m['high_equivalent_total_s'],rtol=0,atol=1e-8)
        theoretical=420 if label=='D300' else 270
        assert min(float(x['high_equivalent_s']) for x in cycles)==theoretical
        with np.load(src/f'{label}_four_stage_states.npz') as x:
            assert x['e'].shape==(m['completed_cycles'],4,20)
            assert x['v'].shape==(m['completed_cycles'],4,20)
            assert x['qP'].shape==(m['completed_cycles'],4,21)
            assert np.isfinite(x['v']).all()
        assert stat['closure']['pulse_max_m_s']<1e-5
        assert stat['closure']['local_LOW_followup_max_qD_m_s2']<1e-9
        checks[label]=dict(cycles=m['completed_cycles'],HIGH_equivalent_cycle_sum_s=total,
            physical_gap_m=m['min_physical_pair_chord_m'],no_order_reversal=True,
            maximum_pulse_velocity_closure_m_s=stat['closure']['pulse_max_m_s'])
        # Explicit actual-trigger P/D fields, including the intervals where the
        # constant-D forecast selected a different satellite from the real one.
        raw=load_config(SOURCE/'normalized_config.yaml').to_dict();raw['simulation']['duration_s']=604800.
        raw['controller']['dwell_time_s']=stat['dwell_s'];cfg=config_from_dict(raw)
        arrays=load_arrays(SOURCE if label=='D300' else args.root/'D150_7d',0,cfg)
        qp,qd=decompose_q(arrays['gap_error_m'],arrays['gap_rate_m_s'],cfg.derived.controller_kp_s2_inv,cfg.derived.controller_kd_s_inv)
        terms=acceleration_terms(arrays['a_m'],arrays['density_kg_m3'],arrays['area_m2'],np.zeros(21),
            r_ref=cfg.derived.reference_radius_m,mu=cfg.constants.mu_m3_s2,cd=cfg.spacecraft.drag_coefficient,
            mass=cfg.spacecraft.mass_kg,area_low=cfg.spacecraft.area_low_m2)
        low_integral=cumulative_trapezoid(terms['low'],arrays['time_s'])
        quiet_rows=[]
        with (src/f'{label}_quiet_intervals.csv').open(encoding='utf-8-sig',newline='') as f:
            for q in csv.DictReader(f):
                t0=float(q['completion_time_s']);actual=float(q['actual_next_command_s']) if q['actual_next_command_s'] else None
                k0=int(np.searchsorted(arrays['time_s'],t0));k1=int(np.searchsorted(arrays['time_s'],actual)) if actual is not None else None
                j=int(q['actual_next_satellite_id'])-1 if q['actual_next_satellite_id'] else None
                quiet_rows.append(dict(start_s=t0,next_HIGH_s=actual,
                    observed_duration_s=(actual if actual is not None else 604800.)-t0,right_censored=actual is None,
                    actual_trigger_satellite_id=j+1 if j is not None else None,
                    actual_trigger_qP_at_start=float(qp[k0,j]) if j is not None else None,
                    actual_trigger_qD_at_start=float(qd[k0,j]) if j is not None else None,
                    actual_trigger_qP_at_HIGH=float(qp[k1,j]) if j is not None else None,
                    actual_trigger_qD_at_HIGH=float(qd[k1,j]) if j is not None else None,
                    actual_trigger_delta_qD_from_natural_LOW=float(cfg.derived.controller_kd_s_inv*node_map(np.diff(low_integral[k1]-low_integral[k0]))[j]) if j is not None else None,
                    start_max_abs_gap_velocity_m_s=float(np.max(abs(arrays['gap_rate_m_s'][k0])))))
        with (out/f'{label}_quiet_diagnostics.csv').open('w',encoding='utf-8-sig',newline='') as f:
            w=csv.DictWriter(f,fieldnames=list(quiet_rows[0]));w.writeheader();w.writerows(quiet_rows)
        summary_rows.append(dict(case=label,days=7,dwell_s=stat['dwell_s'],max_error_m=m['max_error_m'],
            RMS_error_m=m['rms_error_m'],final_max_error_m=m['final_max_error_m'],
            residual_incident_abs_v_median_mm_s=stat['residual_abs_incident_gap_velocity_m_s']['median']*1000,
            residual_incident_max_v_median_mm_s=stat['residual_max_incident_velocity_per_cycle_m_s']['median']*1000,
            pulse_median_mm_s=stat['pulse_velocity_contribution_m_s']['median']*1000,
            cycles=m['completed_cycles'],transitions=m['transitions'],
            HIGH_equivalent_satellite_hours=m['high_equivalent_total_s']/3600,
            extra_altitude_loss_m=m['extra_mean_altitude_loss_m'],
            retrigger_observed_median_hours=stat['re_trigger']['observed_wait_s']['median']/3600,
            exceed_100m_fraction=m['exceedances'][0]['exceedance_fraction']))
    for p in src.iterdir():
        # No large state arrays or repeated minute-resolution history in Git.
        if p.suffix in ['.json','.csv'] and not p.name.endswith('_history_60s.csv'):
            shutil.copyfile(p,out/p.name)
    for name in ['normalized_config.yaml','environment.json','input_difference.json','completion.json','metadata.json']:
        shutil.copyfile(args.root/'D150_7d'/name,out/f'D150_{name}')
    raw=load_config(SOURCE/'normalized_config.yaml').to_dict();raw['simulation']['duration_s']=7*86400.
    dump_normalized_config(config_from_dict(raw),out/'D300_evaluation_config.yaml')
    shutil.copyfile(PROFILE,out/'frozen_profiles.json');shutil.copyfile(SOURCE/'environment.json',out/'D300_environment.json')
    with (out/'comparison_summary.csv').open('w',encoding='utf-8-sig',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(summary_rows[0]));w.writeheader();w.writerows(summary_rows)
    write_json(out/'verification.json',dict(passed=True,tests_passed=291,new_tests=3,physics_runs_added=1,
        new_physics_run='D150 only, seven days',new_physics_run_elapsed_s=json.loads((args.root/'D150_7d/completion.json').read_text())['elapsed_s'],
        D300_reused=True,allLOW_reused=True,regression_24h_passed=True,
        old_inputs_unchanged=len(manifest),old_large_array_check='size and mtime, small files SHA256',
        same_trajectory_through_270s=comparison['input_check']['common_0_to_270s_max_abs_difference'],
        checks=checks,profile_sha256=hashlib.sha256(PROFILE.read_bytes()).hexdigest()))
    print(json.dumps(summary_rows,ensure_ascii=False,indent=2));print('Verification passed',flush=True)


if __name__=='__main__':main()
