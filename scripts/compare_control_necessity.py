#!/usr/bin/env python3
"""Compare four saved 60-day runs; never call the simulator or overwrite output."""
import argparse
import csv
import json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from sat_string.config import load_config
from sat_string.evaluation import switching_metrics
from sat_string.control_cost import event_activity, clean_json
from sat_string.io import load_metadata, load_switch_events, prepare_output_directory, write_json
from sat_string.necessity import exceedance_metrics, improvement, validate_comparison_config


def csv_rows(path, rows):
    with path.open('w', encoding='utf-8-sig', newline='') as stream:
        w = csv.DictWriter(stream, fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ['all-low', 'kp03', 'baseline', 'kp3']:
        p.add_argument('--'+name, type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    args = p.parse_args()
    cases = [('all_low', '全LOW（無制御）', args.all_low, 1., False),
             ('kp03', '制御：Kp＝0.3倍', args.kp03, .3, True),
             ('baseline', '制御：Kp＝1倍', args.baseline, 1., True),
             ('kp3', '制御：Kp＝3倍', args.kp3, 3., True)]
    reference = load_config(args.all_low / 'normalized_config.yaml').to_dict()
    frozen = load_config('tests/fixtures/regression_24h/config.yaml').to_dict()
    frozen['simulation']['duration_s'] = 60*86400
    validate_comparison_config(reference, frozen, 1.)
    out = prepare_output_directory(args.output)
    rows, thresholds, gaps, checkpoints, checks = [], [], [], [], []
    initial_keys = ['a_m', 'lambda_rad', 'area_m2', 'mode', 'density_kg_m3']
    for key, label, source, factor, controlled in cases:
        cfg = load_config(source / 'normalized_config.yaml').to_dict()
        validate_comparison_config(cfg, reference, factor)
        meta = load_metadata(source)
        if meta['configuration'] != cfg or not meta['case']['disturbance_enabled'] or meta['case']['control_enabled'] != controlled:
            raise ValueError('Metadata/configuration or case flags mismatch')
        events = load_switch_events(source)
        with np.load(source / 'results.npz', allow_pickle=False) as data:
            t, e, a = data['time_s'], data['gap_error_m'], data['a_m']
            final_rates = data['gap_rate_m_s'][-1].copy()
            initials = {k: data[k][0].copy() for k in initial_keys}
            if key == 'all_low':
                ref_t, ref_initial, low_final = t.copy(), initials, a[-1].copy()
                if events or not np.all(data['area_m2'] == cfg['spacecraft']['area_low_m2']) or not np.all(data['mode'] == 0):
                    raise ValueError('All LOW saved run contains actuation')
            if not np.array_equal(t, ref_t) or any(not np.array_equal(initials[k], ref_initial[k]) for k in initial_keys):
                raise ValueError('Time grid or actual initial state mismatch')
            if t[0] != 0 or t[-1] != 60*86400 or e.shape != (len(t), 20) or not np.all(np.isfinite(e)):
                raise ValueError('Invalid horizon, shape or nonfinite errors')
            envelope = np.max(np.abs(e), axis=1)
            instant_rms = np.sqrt(np.mean(e*e, axis=1))
            switching = switching_metrics(events, 21, t[-1])
            exposure = event_activity(events, 21, t[-1])
            row = dict(case=key, label=label, source_directory=source.as_posix(),
                max_error_m=float(envelope.max()), rms_error_m=float(np.sqrt(np.mean(e*e))),
                final_max_error_m=float(envelope[-1]),
                transitions=switching['n_transition_starts'],
                completed_cycles=sum(switching['lhl_cycles_per_satellite']),
                last_switch_time_s=switching['last_switch_time_s'],
                high_stable_total_s=exposure['high_stable_total_s'],
                high_equivalent_total_s=exposure['high_equivalent_total_s'],
                initial_mean_altitude_loss_m=float(np.mean(a[0]-a[-1])),
                extra_altitude_loss_m=float(np.mean(low_final-a[-1])),
                minimum_spacing_m=float(100000+e.min()),
                maximum_spacing_m=float(100000+e.max()))
            for metric in ['max_error', 'rms_error', 'final_max_error']:
                row[metric+'_percent_of_target'] = row[metric+'_m']/100000*100
            # Check against the existing formal sample-RMS metrics, not older time-weighted RMS.
            saved = json.loads((source/'metrics.json').read_text(encoding='utf-8'))
            for new, old in [('max_error_m','max_abs_error_m'), ('rms_error_m','rms_error_m'), ('final_max_error_m','final_max_abs_error_m')]:
                if not np.isclose(row[new], saved['formation'][old], rtol=1e-12, atol=1e-9):
                    raise ValueError('Saved formation metrics mismatch')
            rows.append(row)
            for tol in [50, 100, 200, 500, 1000]:
                thresholds.append(dict(case=key, **exceedance_metrics(t, envelope, tol)))
            for i in range(20):
                gaps.append(dict(case=key, gap_number=i+1, max_abs_error_m=float(np.abs(e[:,i]).max()),
                    rms_error_m=float(np.sqrt(np.mean(e[:,i]**2))), final_error_m=float(e[-1,i]),
                    final_rate_m_s=float(final_rates[i])))
            for day in [1, 3, 7, 14, 21, 30, 45, 60]:
                j=int(np.searchsorted(t, day*86400))
                checkpoints.append(dict(case=key, day=day, instantaneous_max_error_m=float(envelope[j]),
                    instantaneous_rms_error_m=float(instant_rms[j])))
            # Only a small plotting cache; full-resolution science stays in source NPZ.
            ix=np.unique(np.r_[np.arange(0,len(t),max(1,len(t)//20000)),len(t)-1,envelope.argmax()])
            np.savez_compressed(out/(key+'_plot.npz'), time_s=t[ix], max_error_m=envelope[ix],
                rms_error_m=instant_rms[ix], gap_error_m=e[ix][:,[0,8,9,10,11,19]])
            write_json(out/(key+'_switching.json'), clean_json(dict(switching=switching, exposure=exposure)))
            write_json(out/(key+'_input.json'), dict(configuration=cfg, case=meta['case']))
            checks.append(dict(case=key, configuration=True, flags=True, initial_states=True,
                time_grid=True, saved_metrics=True, all_low_no_actuation=True if not controlled else None))
        print(label, row['max_error_m'], row['rms_error_m'], flush=True)
    reductions=[]
    for row in rows[1:]:
        for metric in ['max_error_m','rms_error_m']:
            reductions.append(dict(case=row['case'], metric=metric, **improvement(row[metric],rows[0][metric])))
    write_json(out/'comparison.json',dict(cases=rows))
    write_json(out/'input_checks.json',dict(passed=True, cases=checks,
        initial_array_keys=initial_keys, allowed_differences=['Kp override/derived Kp', 'control_enabled', 'case name'],
        new_physics_runs=0, optional_no_disturbance_60d='not present; no new run'))
    for name, data in [('comparison',rows),('exceedance',thresholds),('gap_metrics',gaps),
                       ('checkpoints',checkpoints),('reductions',reductions)]: csv_rows(out/(name+'.csv'),data)


if __name__ == '__main__': main()
