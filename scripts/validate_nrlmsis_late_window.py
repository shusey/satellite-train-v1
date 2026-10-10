#!/usr/bin/env python3
"""One-day 10/120 s convergence at the first maximum-Ap day in the frozen period."""
import argparse
from dataclasses import replace
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from sat_string.config import load_config
from sat_string.continuation import continue_case
from sat_string.io import prepare_output_directory,write_json
from sat_string.simulator import CaseDefinition
from run_nrlmsis_case import make_environment


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run',type=Path,required=True);p.add_argument('--weather',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);args=p.parse_args();out=prepare_output_directory(args.output)
    config=load_config(args.run/'normalized_config.yaml')
    _,backend,_=make_environment(Path('configs/nrlmsis2_environment.yaml'),args.weather,0)
    weather=backend.weather
    peak=int(np.argmax(weather.values[:,3]))  # ap1 is column3; first maximum, not chosen from response.
    start=float(np.floor(weather.seconds[peak]/86400)*86400);end=start+86400
    with np.load(args.run/'results.npz') as d:
        k=int(np.searchsorted(d['time_s'],start));snapshot={key:d[key][k:k+1] for key in d.files}
    source=replace(config,simulation=replace(config.simulation,duration_s=start))
    destination=replace(config,simulation=replace(config.simulation,duration_s=end))
    records=[];reference=None
    for step in [10,120]:
        _,backend,provider=make_environment(Path('configs/nrlmsis2_environment.yaml'),args.weather,step)
        result=continue_case(snapshot,source,destination,CaseDefinition('nrlmsis2_all_low',False,False),[],density_provider=provider)
        e=result.gap_error_m;v=result.gap_rate_m_s
        if reference is None:reference=(e.copy(),v.copy())
        row=dict(cache_s=step,max_abs_error_m=float(np.abs(e).max()),rms_error_m=float(np.sqrt(np.mean(e*e))),
            max_error_history_difference_m=float(np.max(np.abs(e-reference[0]))),
            max_gap_rate_history_difference_m_s=float(np.max(np.abs(v-reference[1]))))
        records.append(row)
        np.savez_compressed(out/f'window_{step}s.npz',**result.arrays())
        print(row,flush=True)
    report=dict(start_s=start,end_s=end,first_maximum_ap1=float(weather.values[peak,3]),
        selection='first maximum observed ap1 in frozen weather, independent of formation response',
        same_initial_physical_state=True,control_enabled=False,records=records,
        passed=records[-1]['max_error_history_difference_m']<=.05,
        maximum_acceptable_error_history_difference_m=.05)
    write_json(out/'late_convergence.json',report)
    if not report['passed']:raise RuntimeError('Late-window convergence failed')


if __name__=='__main__':main()
