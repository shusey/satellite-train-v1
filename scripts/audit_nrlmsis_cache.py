#!/usr/bin/env python3
"""Direct density audits at saved positions and inside cache blocks; no physics rerun."""
import argparse
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'src'))
from sat_string.config import load_config
from sat_string.io import write_json
from sat_string.nrlmsis import TrajectoryDensityCache
from run_nrlmsis_case import make_environment


def audit(run, weather, step, probes=250):
    _, backend, _ = make_environment(Path('configs/nrlmsis2_environment.yaml'), weather, 0)
    mu = load_config(run/'normalized_config.yaml').constants.mu_m3_s2
    with np.load(run/'results.npz') as d:
        t, a, u, rho = (d[k] for k in ['time_s','a_m','lambda_rad','density_kg_m3'])
    jumps=np.arange(10800,t[-1]+1,10800)
    ix = np.unique(np.r_[np.linspace(0,len(t)-1,probes,dtype=int),
        np.searchsorted(t,jumps),np.searchsorted(t,jumps)-1])
    direct = backend(t[ix,None],a[ix],u[ix])
    saved = np.abs(rho[ix]/direct-1)
    # Query mid-output instants, including first/last blocks before daily jumps.
    anchors = np.unique(np.r_[np.floor(t[ix[:-1]]/step)*step,jumps-step,jumps])
    errors = []
    for anchor in anchors:
        query = anchor + min(5., step/2)
        if query > t[-1]: continue
        def state(at, values):
            k=min(int(np.searchsorted(t,at,side='right'))-1,len(t)-2)
            fraction=(at-t[k])/(t[k+1]-t[k])
            return values[k]+fraction*(values[k+1]-values[k])
        aa, uu = state(anchor,a), state(anchor,u)
        aq, uq = state(query,a), state(query,u)
        cache = TrajectoryDensityCache(backend,step,mu)
        cache(anchor,aa,uu)
        errors.extend(np.abs(cache(query,aq,uq)/backend(query,aq,uq)-1))
    errors = np.asarray(errors)
    return dict(cache_s=step, saved_positions_max_relative_error=float(saved.max()),
        saved_positions_rms_relative_error=float(np.sqrt(np.mean(saved**2))),
        interior_max_relative_error=float(errors.max()), interior_rms_relative_error=float(np.sqrt(np.mean(errors**2))),
        interior_probe_points=len(errors), saved_probe_points=int(direct.size),
        passed=bool(max(saved.max(),errors.max())<.001),
        interior_state_convention='Same linearly interpolated saved state for direct/cache; no claim of exact physical state between outputs',
        maximum_acceptable_relative_density_error=.001)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--run',type=Path,required=True);p.add_argument('--weather',type=Path,required=True)
    p.add_argument('--cache-seconds',type=int,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    if args.output.exists(): raise FileExistsError(args.output)
    result=audit(args.run,args.weather,args.cache_seconds)
    write_json(args.output,result);print(result)
    if not result['passed']:raise RuntimeError('Direct density audit failed')


if __name__=='__main__':main()
