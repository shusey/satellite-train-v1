"""Direct replay at ALL saved positions, not another physical simulation.

Measures cumulative differential cache-density error, rather than accepting
pointwise common-mode density accuracy alone. Frozen-trajectory correction
does not include dynamical/density feedback and is not a rigorous error bound.
"""
import argparse
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from sat_string.config import load_config
from sat_string.io import write_json
from run_nrlmsis_case import make_environment


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ('run','environment','weather','output'):p.add_argument('--'+name,type=Path,required=True)
    args=p.parse_args()
    if args.output.exists():raise FileExistsError(args.output)
    cfg=load_config(args.run/'normalized_config.yaml');_,backend,_=make_environment(args.environment,args.weather,0)
    with np.load(args.run/'results.npz') as d:
        t,a,u,rho,e,v=[d[k] for k in ('time_s','a_m','lambda_rad','density_kg_m3','gap_error_m','gap_rate_m_s')]
    delta=np.empty_like(rho);max_relative=0.;max_adjacent=0.
    for start in range(0,len(t),5000):
        stop=min(len(t),start+5000)
        direct=backend(t[start:stop,None],a[start:stop],u[start:stop])
        delta[start:stop]=rho[start:stop]-direct
        max_relative=max(max_relative,float(np.abs(rho[start:stop]/direct-1).max()))
        max_adjacent=max(max_adjacent,float(np.abs(np.diff(delta[start:stop],axis=1)).max()))
    cumulative=np.r_[np.zeros((1,a.shape[1])),np.cumsum(.5*(delta[:-1]+delta[1:])*np.diff(t)[:,None],axis=0)]
    k=cfg.spacecraft.drag_coefficient*cfg.spacecraft.area_low_m2*np.sqrt(cfg.constants.mu_m3_s2)/(2*cfg.spacecraft.mass_kg)
    # If cached density exceeds direct density, direct-forcing radius is larger.
    corrected_radius=(np.sqrt(a)+k*cumulative)**2
    corrected_rate=cfg.derived.reference_radius_m*np.diff(np.sqrt(cfg.constants.mu_m3_s2/corrected_radius**3),axis=1)
    rate_difference=v-corrected_rate
    spacing_difference=np.r_[np.zeros((1,v.shape[1])),np.cumsum(.5*(rate_difference[:-1]+rate_difference[1:])*np.diff(t)[:,None],axis=0)]
    max_spacing=float(np.abs(spacing_difference).max());allowance=max(.05,.01*float(np.abs(e).max()))
    result=dict(all_saved_points=int(rho.size),saved_step_s=float(np.diff(t).max()),
      max_relative_density_error=max_relative,max_adjacent_density_error_kg_m3=max_adjacent,
      max_cumulative_adjacent_density_error_kg_s_m3=float(np.abs(np.diff(cumulative,axis=1)).max()),
      final_adjacent_density_error_kg_s_m3=np.diff(cumulative[-1]).tolist(),
      max_estimated_gap_rate_error_m_s=float(np.abs(rate_difference).max()),max_estimated_spacing_error_m=max_spacing,
      acceptance_estimated_spacing_error_m=allowance,passed=max_relative<.001 and max_spacing<allowance,
      method='Direct NRLMSIS2.0 replay of saved10s positions, same drivers; trapezoid cache-minus-direct density, frozen-trajectory radius correction',
      limitations='No new coupled physics; no rigorous error bound; common-mode and phase/radius feedback not repropagated; original2s RK jumps considered in separate quadrature diagnostic')
    write_json(args.output,result);print(result,flush=True)
    if not result['passed']:raise RuntimeError('Integrated differential cache audit needs investigation')


if __name__=='__main__':main()
