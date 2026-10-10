#!/usr/bin/env python3
"""Short validation and cache convergence before any 60-day run."""
import argparse
from pathlib import Path
import sys
import subprocess
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from sat_string.io import prepare_output_directory,write_json
from run_nrlmsis_case import make_environment


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--weather',type=Path,required=True);p.add_argument('--case',choices=['all_low','baseline'],required=True)
    p.add_argument('--direct-run',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();out=prepare_output_directory(args.output)
    _,backend,_=make_environment(Path('configs/nrlmsis2_environment.yaml'),args.weather,0)
    reference=np.load(args.direct_run/'results.npz');rt=reference['time_s'];re=reference['gap_error_m']
    report=[]
    for step in [10,30,60,120]:
        run=out/f'{args.case}_{step}s'
        subprocess.run([sys.executable,'scripts/run_nrlmsis_case.py','--weather',str(args.weather),'--case',args.case,
            '--duration-s','86400','--cache-seconds',str(step),'--output',str(run)],check=True)
        with np.load(run/'results.npz') as d:
            e=d['gap_error_m'];t=d['time_s'];rho=d['density_kg_m3']
            ix=np.unique(np.r_[np.arange(0,len(t),31),len(t)-1])
            a,u=d['a_m'][ix],d['lambda_rad'][ix]
            direct=backend.evaluate(t[ix,None],a,u)
            relative=np.abs(rho[ix]/direct-1)
            short=e[:len(rt)]
            record=dict(cache_s=step,max_error_m=float(np.abs(e).max()),rms_error_m=float(np.sqrt(np.mean(e*e))),
                short_max_error_m=float(np.abs(short).max()),short_rms_error_m=float(np.sqrt(np.mean(short*short))),
                short_max_history_difference_m=float(np.max(np.abs(short-re))),
                direct_density_max_relative_error=float(relative.max()),direct_density_rms_relative_error=float(np.sqrt(np.mean(relative**2))),
                run=str(run))
            report.append(record)
    ref=report[0]
    direct_max=float(np.abs(re).max());direct_rms=float(np.sqrt(np.mean(re*re)))
    for row in report:
        row['passed']=(abs(row['max_error_m']-ref['max_error_m'])<=max(.05,.01*ref['max_error_m'])
            and abs(row['rms_error_m']-ref['rms_error_m'])<=max(.05,.01*ref['rms_error_m'])
            and abs(row['short_max_error_m']-direct_max)<=max(.05,.01*direct_max)
            and abs(row['short_rms_error_m']-direct_rms)<=max(.05,.01*direct_rms)
            and row['direct_density_max_relative_error']<.001)
    selected=max((r['cache_s'] for r in report if r['passed']),default=None)
    write_json(out/'convergence.json',dict(case=args.case,records=report,selected_cache_s=selected,
        direct_duration_s=float(rt[-1]),duration_s=86400,direct_reference=str(args.direct_run),
        selection_basis='largest tested interval satisfying pre-registered density and max/RMS tolerances'))
    print('Selected',args.case,selected,flush=True)
    if selected is None:raise RuntimeError('No validated cache interval')


if __name__=='__main__':main()
