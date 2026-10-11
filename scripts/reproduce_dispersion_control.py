"""Fresh-checkout reproduction only: three frozen ON profiles and matched OFF.

Do not use when the published local runs already exist: reuse them. This helper
is not executed as an extra experiment in the study. It avoids rebuilding the
old37-profile screening solely to obtain two reference trajectories.
"""
import argparse,hashlib,json
from pathlib import Path
from run_dispersion_control import ROOT,configuration,frozen_profiles
from run_nrlmsis_case import make_environment
from sat_string.dispersion_control import integrate_control_profiles
from sat_string.config import dump_normalized_config
from sat_string.io import write_json


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--days',type=int,choices=[30,60],required=True)
    parser.add_argument('--output-root',type=Path,required=True)
    parser.add_argument('--prefix',type=Path)
    args=parser.parse_args();spec,cfg=configuration(args.days)
    on=frozen_profiles(spec)[:3]
    profiles=on+[dict(p,case=p['case']+'_OFF',control_enabled=False) for p in on]
    env,backend,provider=make_environment(ROOT/spec['environment'],ROOT/spec['weather'],spec['cache_seconds'])
    weather_hash=hashlib.sha256((ROOT/spec['weather']).read_bytes()).hexdigest()
    checkpoint=None
    if args.days==60 and not args.prefix:raise ValueError('Review30d first;60d needs checkpoint')
    if args.prefix:
        if not (args.prefix/'completion.json').exists():raise ValueError('Prefix incomplete')
        checkpoint=json.loads((args.prefix/'checkpoint.json').read_text())
        previous=json.loads((args.prefix/'environment.json').read_text())
        if previous['environment']!=env or previous['weather_sha256']!=weather_hash or previous['cache_seconds']!=spec['cache_seconds']:
            raise ValueError('Changed space weather/geometry/cache')
    target=args.output_root/f'batch_{args.days}d'
    progress=lambda t,end,wall,counts,gap:print(f'{t/86400:g}/{end/86400:g}d; {wall:.1f}s; commands={counts}; signed minimum={gap:g}m',flush=True)
    integrate_control_profiles(cfg,profiles,provider,target,checkpoint=checkpoint,progress=progress)
    dump_normalized_config(cfg,target/'normalized_config.yaml')
    write_json(target/'environment.json',dict(environment=env,weather_sha256=weather_hash,cache_seconds=spec['cache_seconds'],
        pymsis_version='0.12.0',msis_version='2.0',density_units='kg/m3',backend_points=backend.points_evaluated,
        maximum_cache_radius_offset_m=provider.maximum_radius_offset_m,maximum_cache_phase_offset_rad=provider.maximum_phase_offset_rad,
        prefix=str(args.prefix) if args.prefix else None,fresh_checkout_reproduction=True))


if __name__=='__main__':main()
