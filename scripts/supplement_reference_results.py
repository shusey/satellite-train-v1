"""Publish common-window, low/reference and signed trigger diagnostics offline."""
import argparse,csv,json
import numpy as np
import yaml
from run_reference_dispersion import ROOT,SPEC,DATA
from analyze_dispersion_control import csv_rows
from sat_string.io import write_json


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--days',type=int,choices=[30,60],required=True)
    days=ap.parse_args().days;target=DATA/f'{days}d'
    out=ROOT/yaml.safe_load(SPEC.read_text())['output']
    comparisons=[];old_new=[];segments=[];triggers=[];corrections=[]
    profile_meta=json.loads((DATA/'profile_metadata.json').read_text())
    for c in ['A','B','AB']:
        on=json.loads((target/f'{c}_metrics.json').read_text());off=json.loads((target/f'{c}_OFF_metrics.json').read_text())
        comparisons.append(dict(case=c,days=days,ON_max_error_m=on['max_error_m'],ON_RMS_m=on['rms_error_m'],
            ON_final_error_m=on['final_max_error_m'],ON_max_gap_rate_m_s=on['max_gap_rate_m_s'],
            OFF_full_fixed_label_max_diagnostic_m=off['max_error_m'],
            OFF_max_before_order_reversal_m=off['max_error_before_order_reversal_m'],
            OFF_reversal_day=off['first_label_order_reversal_s']/86400 if off['first_label_order_reversal_s'] is not None else None,
            common_window_max_error=on['matched_pre_allLOW_reversal_improvement'],
            common_window_RMS=on['matched_pre_allLOW_reversal_RMS_improvement'],
            transitions=on['transitions'],completed_cycles=on['completed_cycles'],
            extra_mean_altitude_loss_m=on['extra_mean_altitude_loss_m']))
        for s in on['segments']:segments.append(dict(case=c,**s))
        low=json.loads((ROOT/f'docs/dispersion_control_transfer/data/{days}d/{c}_metrics.json').read_text())
        shape=c!='AB'
        old_new.append(dict(case=c,old_profiles_same_shape=shape,
            old_a_sample_std_m=profile_meta['old_low'][['A','B','AB'].index(c)]['a_sample_std_m'],
            old_B_sample_std=profile_meta['old_low'][['A','B','AB'].index(c)]['B_sample_std'],
            old_control_max_error_m=low['max_error_m'],reference_control_max_error_m=on['max_error_m'],
            old_cycles=low['completed_cycles'],reference_cycles=on['completed_cycles'],
            old_transitions=low['transitions'],reference_transitions=on['transitions'],
            old_extra_altitude_loss_m=low['extra_mean_altitude_loss_m'],reference_extra_altitude_loss_m=on['extra_mean_altitude_loss_m']))
        def read_csv(name):
            with (target/name).open(encoding='utf-8-sig',newline='') as f:return list(csv.DictReader(f))
        local=read_csv(f'{c}_local_cycles.csv');quiet=read_csv(f'{c}_quiet_cycles.csv') if (target/f'{c}_quiet_cycles.csv').exists() else []
        # Full endpoint arrays/wide cycle tables stay local. Publish a compact,
        # explicitly sampled table as a navigation aid, never for recomputing totals.
        for name,records in [('local',local),('quiet',quiet)]:
            if records:
                ix=sorted(set(range(min(21,len(records))))|set(range(0,len(records),100))|
                          set(range(max(0,len(records)-21),len(records))))
                csv_rows(target/f'{c}_{name}_cycle_samples.csv',
                    [dict(full_table_row=i,sampling_note='first21 + every100 + last21; totals use full events',**records[i]) for i in ix])
        first=local[0] if local else None
        corrections.append(dict(case=c,first_local_LOW_return_s=float(first['low_complete_s']) if first else None,
            first_local_return_max_gap_rate_m_s=float(first['V_res_all_m_s']) if first else None,
            first_local_next_HIGH_s=float(first['actual_next_command_s']) if first and first['actual_next_command_s'] else None,
            first_allLOW_boundary_s=float(quiet[0]['completion_time_s']) if quiet else None,
            longest_observed_allLOW_interval_s=max((float(r['observation_window_end_s'])-float(r['completion_time_s']) for r in quiet),default=None),
            latter_segment_transitions=on['segments'][-1]['transitions'],
            quiet_candidate=on['quiet_candidate'],last_completion_s=on['last_switch_time_s']))
        cache=dict(np.load(out/f'{c}_plot_{days}d.npz'))
        checkpoint=json.loads((out/f'batch_{days}d/checkpoint.json').read_text())
        events=checkpoint['events'][['A','B','AB'].index(c)]
        high=[e for e in events if e['to_mode']=='H']
        chosen=high[:1]+next(([e] for e in high if e['command_time']>=days*.5*86400),[])
        for event in chosen:
            at=event['command_time'];j=event['satellite_id']-1
            k=int(np.searchsorted(cache['time_s'],at));assert cache['time_s'][k]==at
            previous=[e for e in events if e['satellite_id']==j+1 and e['complete_time']<=at]
            last=max((e['complete_time'] for e in previous),default=None)
            incident=[g for g in (j-1,j) if 0<=g<20]
            qp=float(cache['qP'][k,j]);qd=float(cache['qD'][k,j])
            triggers.append(dict(case=c,command_time_s=at,satellite_id=j+1,
                qP_m_s2=qp,qD_m_s2=qd,q_m_s2=qp+qd,pre_command_mode='LOW',
                post_command_mode=int(cache['mode'][k,j]),last_completion_s=last,
                dwell_release_s=last+300 if last is not None else 0,
                pre_command_dwell_remaining_s=max(0,last+300-at) if last is not None else 0,
                incident_gap_ids=[g+1 for g in incident],error_m=cache['on_error'][k,incident].tolist(),
                rate_m_s=cache['on_rate'][k,incident].tolist()))
            ix=np.flatnonzero((cache['time_s']>=max(0,at-3600))&(cache['time_s']<=at+7200))
            csv_rows(target/f'{c}_q_satellite{j+1}_at{int(at)}s.csv',[
                dict(time_s=float(cache['time_s'][i]),satellite_id=j+1,qP=float(cache['qP'][i,j]),
                    qD=float(cache['qD'][i,j]),q=float(cache['qP'][i,j]+cache['qD'][i,j]),
                    mode=int(cache['mode'][i,j]),area_m2=float(cache['area'][i,j])) for i in ix])
    csv_rows(target/'common_window_comparison.csv',comparisons)
    csv_rows(target/'low_reference_comparison.csv',old_new)
    csv_rows(target/'interval_metrics.csv',segments)
    csv_rows(target/'representative_triggers.csv',triggers)
    write_json(target/'reference_summary.json',dict(comparison=comparisons,low_reference=old_new,
        representative_triggers=triggers,segments=segments,correction_diagnostics=corrections))


if __name__=='__main__':main()
