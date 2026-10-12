"""Frozen reference profiles with unchanged controller, dynamics and cache method."""
import argparse, hashlib, json
from pathlib import Path
import yaml
from run_dispersion_control import ROOT, configuration
from run_nrlmsis_case import make_environment
from sat_string.reference_condition import make_reference, load_frozen
from sat_string.dispersion_control import integrate_control_profiles
from sat_string.config import dump_normalized_config
from sat_string.io import write_json

SPEC=ROOT/'configs/reference_dispersion.yaml'
DOC=ROOT/'docs/reference_dispersion'; DATA=DOC/'data'


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--days',type=int,choices=[30,60],required=True)
    ap.add_argument('--prefix',type=Path)
    args=ap.parse_args(); spec=yaml.safe_load(SPEC.read_text())
    if (spec['sigma_a_m'],spec['sigma_B_fraction'],spec['std_ddof']) != (.5,.005,1):
        raise ValueError('This runner is the fixed 0.5m / 0.5percent sample-std reference')
    _,cfg=configuration(args.days)
    frozen=ROOT/spec['frozen_profiles']; DATA.mkdir(parents=True,exist_ok=True)
    if not frozen.exists():
        original=ROOT/spec['source_profiles']
        obj=make_reference(json.loads(original.read_text()))
        obj.update(source_file=spec['source_profiles'],source_sha256=hashlib.sha256(original.read_bytes()).hexdigest())
        write_json(frozen,obj)
    obj=load_frozen(frozen); profiles=obj['profiles']
    checkpoint=None
    if args.days==60:
        if not args.prefix: raise ValueError('60d requires a 30d checkpoint')
        decision=json.loads((DATA/'extension_decision.json').read_text())
        if not decision.get('reviewed') or not decision['extend_to_60d']:
            raise ValueError('Reviewed 30d safety/go decision required')
    if args.prefix:
        if not (args.prefix/'completion.json').exists(): raise ValueError('Prefix incomplete')
        checkpoint=json.loads((args.prefix/'checkpoint.json').read_text())
    env,backend,provider=make_environment(ROOT/spec['environment'],ROOT/spec['weather'],spec['cache_seconds'])
    whash=hashlib.sha256((ROOT/spec['weather']).read_bytes()).hexdigest()
    old=json.loads((ROOT/spec['baseline_run']/'density_environment.json').read_text())
    if old['environment']!=env or old['weather_sha256']!=whash or old['cache_seconds']!=spec['cache_seconds']:
        raise ValueError('Frozen active environment differs')
    if args.prefix:
        prior=json.loads((args.prefix/'environment.json').read_text())
        if prior['environment']!=env or prior['weather_sha256']!=whash or prior['cache_seconds']!=spec['cache_seconds']:
            raise ValueError('Continuation environment differs')
    output=ROOT/spec['output']/f'batch_{args.days}d'
    callback=lambda t,end,wall,counts,gap:print(
        f'{t/86400:g}/{end/86400:g}d; elapsed={wall:.1f}s; event counts={counts}; minimum signed gap={gap:.3f}m',flush=True)
    integrate_control_profiles(cfg,profiles,provider,output,checkpoint=checkpoint,progress=callback)
    dump_normalized_config(cfg,output/'normalized_config.yaml')
    write_json(output/'environment.json',dict(environment=env,weather_sha256=whash,
        cache_seconds=spec['cache_seconds'],pymsis_version='0.12.0',msis_version='2.0',
        density_units='kg/m3',backend_points=backend.points_evaluated,
        maximum_cache_radius_offset_m=provider.maximum_radius_offset_m,
        maximum_cache_phase_offset_rad=provider.maximum_phase_offset_rad,
        frozen_profiles_sha256=hashlib.sha256(frozen.read_bytes()).hexdigest(),
        prefix=str(args.prefix) if args.prefix else None))
    print('Saved',output,flush=True)


if __name__=='__main__': main()
