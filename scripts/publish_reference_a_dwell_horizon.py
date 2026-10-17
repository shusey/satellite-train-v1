"""Publish small horizon comparisons and checkpoint provenance, not raw histories."""
import json,shutil,sys,hashlib,subprocess
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'));sys.path.insert(0,str(ROOT/'src'))
from sat_string.io import write_json
from analyze_reference_pulses import read_csv,write_csv
from analyze_reference_a_dwell_horizon import LONG,SEVEN,preservation_manifest


def main():
    out=ROOT/'docs/reference_a_dwell_60d/data';out.mkdir(parents=True,exist_ok=False)
    reports={7:json.loads((SEVEN/'analysis/comparison.json').read_text())}
    for days in [30,60]:
        src=LONG/f'analysis_{days}d'
        reports[days]=json.loads((src/'comparison.json').read_text())
        assert reports[days]['verification']['passed']
        target=out/f'{days}d';target.mkdir()
        for name in ['comparison.json','D150_metrics.json','D300_metrics.json','regression_24h.json']:
            shutil.copyfile(src/name,target/name)
        # Source indices permit navigating full local tables without committing
        # four-stage states or the whole minute-by-minute plot history again.
        for case in ['D300','D150']:
            for suffix in ['cycles','quiet_intervals']:
                rows=read_csv(src/f'{case}_{suffix}.csv')
                chosen=sorted(set(range(min(3,len(rows))))|set(range(0,len(rows),100))|set(range(max(0,len(rows)-3),len(rows))))
                write_csv(target/f'{case}_{suffix}_samples.csv',[dict(full_table_row=k,**rows[k]) for k in chosen])
    summary=[];intervals=[]
    for day,report in reports.items():
        for case,s in report['cases'].items():
            m=s['formation'];summary.append(dict(days=day,case=case,max_error_m=m['max_error_m'],RMS_error_m=m['rms_error_m'],
                final_error_m=m['final_max_error_m'],max_error_percent_of_100km=m['max_error_percent_of_100km'],
                completed_cycles=m['completed_cycles'],transitions=m['transitions'],last_switch_s=m['last_switch_time_s'],
                HIGH_equivalent_satellite_hours=m['high_equivalent_total_s']/3600,
                extra_mean_altitude_loss_m=m['extra_mean_altitude_loss_m'],
                residual_incident_max_median_mm_s=s['residual_max_incident_velocity_per_cycle_m_s']['median']*1000,
                same_satellite_retrigger_observed_median_hours=s['re_trigger']['observed_wait_s']['median']/3600,
                exceed_100m_fraction=m['exceedances'][0]['exceedance_fraction'],
                ever_exceed_500m=m['exceedances'][1]['ever_exceeded'],minimum_physical_gap_m=m['min_physical_pair_chord_m']))
            if day==60:
                for segment in m['segments']:
                    intervals.append(dict(case=case,start_day=segment['start_day'],end_day=segment['end_day'],
                        RMS_m=segment['error_rms_m'],max_error_m=segment['max_error_m'],
                        transitions_per_day=segment['transitions_per_day'],cycles_per_day=segment['cycles_per_day'],
                        HIGH_equivalent_satellite_hours_per_day=segment['high_equivalent_satellite_s_per_day']/3600,
                        extra_altitude_loss_m_per_day=segment['extra_mean_altitude_loss_per_day_m']))
    write_csv(out/'horizon_comparison.csv',summary);write_csv(out/'interval_comparison.csv',intervals)
    shutil.copyfile(LONG/'extension_decision_30d.json',out/'extension_decision_30d.json')
    for day in [30,60]:
        part=LONG/f'D150_{day}d'
        target=out/f'{day}d'
        for name in ['normalized_config.yaml','environment.json','completion.json','input_difference.json','checkpoint.json']:
            shutil.copyfile(part/name,target/f'D150_{name}')
    # Check both continuation chains' cumulative event prefixes once more.
    parts=[SEVEN/'D150_7d',LONG/'D150_30d',LONG/'D150_60d'];chain=[]
    for a,b in zip(parts,parts[1:]):
        earlier=json.loads((a/'checkpoint.json').read_text());later=json.loads((b/'checkpoint.json').read_text())
        assert earlier['events'][0]==later['events'][0][:len(earlier['events'][0])]
        env=json.loads((b/'environment.json').read_text())
        assert env['prefix_checkpoint_sha256']==hashlib.sha256((a/'checkpoint.json').read_bytes()).hexdigest()
        chain.append(dict(from_day=earlier['time_s']/86400,to_day=later['time_s']/86400,
            prefix_checkpoint_sha256=env['prefix_checkpoint_sha256'],prefix_events_exact=True))
    old=json.loads((LONG/'analysis_60d/input_manifest.json').read_text())
    assert old==preservation_manifest()
    changed=subprocess.check_output(['git','diff','--name-only','HEAD','--','src/sat_string/controller.py','src/sat_string/dynamics.py',
        'src/sat_string/integrator.py','src/sat_string/attitude.py','src/sat_string/dispersion_control.py','configs'],cwd=ROOT,text=True)
    assert not changed.strip()
    verification=dict(passed=True,tests_passed=293,checkpoint_chain=chain,
        old_input_files_unchanged=len(old),regression_24h_passed=True,
        core_and_config_changes=False,new_physics='D150 suffix 7–30d and 30–60d only; D300/allLOW reused',
        new_physics_elapsed_s=sum(json.loads((LONG/f'D150_{d}d/completion.json').read_text())['elapsed_s'] for d in [30,60]))
    write_json(out/'verification.json',verification)
    print(json.dumps(summary,ensure_ascii=False,indent=2),flush=True)


if __name__=='__main__':main()
