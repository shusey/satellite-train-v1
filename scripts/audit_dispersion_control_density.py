"""Short offline density/cache audit on already saved controlled trajectories."""
import argparse,json
import numpy as np
from run_dispersion_control import ROOT,OUT,DATA,configuration
from run_nrlmsis_case import make_environment
from sat_string.io import write_json


def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--days',type=int,choices=[30,60],required=True)
    args=p.parse_args();spec,cfg=configuration(args.days)
    env,backend,provider=make_environment(ROOT/spec['environment'],ROOT/spec['weather'],spec['cache_seconds'])
    times=[];radii=[];phases=[];densities=[]
    paths=[OUT/f'batch_{day}d' for day in (30,60) if day<=args.days]
    if args.days==60 and (OUT/'OFF_A_B_60d').exists():paths.append(OUT/'OFF_A_B_60d')
    for path in paths:
        t=np.load(path/'time_s.npy')
        a=np.load(path/'a_m.npy',mmap_mode='r');u=np.load(path/'lambda_rad.npy',mmap_mode='r');rho=np.load(path/'density_kg_m3.npy',mmap_mode='r')
        start=t[0];end=t[-1]
        # Interior offsets test interpolation, not merely exact cache nodes.
        selected=np.ravel([base+off for base in np.arange(start,end,21600) for off in (10,50,90)])
        selected=selected[(selected>=start)&(selected<=end)]
        ix=np.searchsorted(t,selected)
        times.append(np.broadcast_to(selected[:,None,None],a[ix].shape).ravel())
        radii.append(a[ix].ravel());phases.append(u[ix].ravel());densities.append(rho[ix].ravel())
    tt=np.concatenate(times);aa=np.concatenate(radii);uu=np.concatenate(phases);saved=np.concatenate(densities)
    direct=backend.evaluate(tt,aa,uu)
    relative=np.abs(saved/direct-1)
    report=dict(points=len(tt),interior_offsets_s=[10,50,90],selection='every6h; all saved profiles and21 satellites',
        audited_sources=[str(p.relative_to(OUT)) for p in paths],
        maximum_relative_density_error=float(relative.max()),RMS_relative_density_error=float(np.sqrt(np.mean(relative**2))),
        maximum_saved_density_kg_m3=float(saved.max()),minimum_saved_density_kg_m3=float(saved.min()),
        msis_version='2.0',pymsis_version='0.12.0',density_units='kg/m3',
        # Existing convergence studies already justify120s. This guards new trajectories only.
        passed=bool(np.all(np.isfinite(saved)) and np.all(saved>0) and relative.max()<.001),
        audit_relative_limit=.001,interpretation='spot audit, not proof over every RK stage')
    write_json(DATA/f'density_audit_{args.days}d.json',report);print(report,flush=True)
    return 0 if report['passed'] else 1


if __name__=='__main__':raise SystemExit(main())
