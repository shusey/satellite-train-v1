"""Storm-time density/interpolation integral test, NO new physical run.

Prescribed400km circular tracks at the first maximum observed3h Ap. Tests
10/30/60/120s cache against direct2s density along exactly the same tracks.
This does not claim a60day coupled trajectory convergence proof.
"""
import argparse
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from sat_string.config import load_config
from sat_string.controller import initial_formation
from sat_string.nrlmsis import TrajectoryDensityCache
from sat_string.io import write_json
from run_nrlmsis_case import make_environment


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--environment',type=Path,required=True);p.add_argument('--weather',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True);args=p.parse_args()
    if args.output.exists():raise FileExistsError(args.output)
    cfg=load_config('tests/fixtures/regression_24h/config.yaml')
    spec,backend,_=make_environment(args.environment,args.weather,0)
    a,u0=initial_formation(21,cfg.derived.reference_radius_m,cfg.derived.target_spacing_rad)
    peak=int(np.argmax(backend.weather.values[:,3]));start=float(backend.weather.seconds[peak])
    t=start+np.arange(0,16680+1,2.);u=u0+cfg.derived.reference_mean_motion_rad_s*t[:,None]
    direct=backend(t[:,None],np.broadcast_to(a,u.shape),u)
    def cumul(values):return np.r_[np.zeros((1,21)),np.cumsum((values[:-1]+values[1:])*np.diff(t)[:,None]/2,axis=0)]
    ref=cumul(direct);records=[]
    coefficient=1.5*cfg.spacecraft.drag_coefficient*cfg.spacecraft.area_low_m2*cfg.constants.mu_m3_s2/(cfg.spacecraft.mass_kg*cfg.derived.reference_radius_m)
    for step in (10,30,60,120):
        cache=TrajectoryDensityCache(backend,step,cfg.constants.mu_m3_s2)
        sampled=np.array([cache(at,a,phase) for at,phase in zip(t,u)])
        jr_error=np.diff(cumul(sampled)-ref,axis=1)
        velocity_error=coefficient*jr_error
        spacing_error=np.sum((velocity_error[:-1]+velocity_error[1:])*np.diff(t)[:,None]/2,axis=0)
        records.append(dict(cache_s=step,max_relative_density_error=float(np.abs(sampled/direct-1).max()),
          max_adjacent_integral_error_kg_s_m3=float(np.abs(jr_error).max()),max_implied_gap_rate_error_m_s=float(np.abs(velocity_error).max()),
          max_implied_end_spacing_error_m=float(np.abs(spacing_error).max())))
    result=dict(epoch=spec['epoch_utc'],start_s=start,end_s=float(t[-1]),ap1=float(backend.weather.values[peak,3]),
       reference='direct NRLMSIS2.0 at2s on prescribed400km circular positions; no coupled physics rerun',
       records=records,passed=records[-1]['max_relative_density_error']<.001 and records[-1]['max_implied_end_spacing_error_m']<.05,
       limitations='Fixed-a density-field quadrature only; actual60d-state density audited separately;3h jump quadrature shared by reference/candidate')
    write_json(args.output,result);print(result,flush=True)
    if not result['passed']:raise RuntimeError('Density integral convergence failed')


if __name__=='__main__':main()
