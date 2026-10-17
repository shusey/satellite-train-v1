"""Fixed A-ref D150 intervention and checkpoint extension; D300 is reused."""
import argparse, hashlib, json, sys
from pathlib import Path
import numpy as np
import yaml
ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/'scripts'));sys.path.insert(0,str(ROOT/'src'))
from run_nrlmsis_case import make_environment
from sat_string.config import load_config,config_from_dict,dump_normalized_config
from sat_string.reference_condition import load_frozen
from sat_string.dispersion_control import integrate_control_profiles
from sat_string.io import write_json

SOURCE=ROOT/'outputs/reference_dispersion_20261012/batch_30d'
SPEC=ROOT/'configs/reference_dispersion.yaml'
PROFILE=ROOT/'docs/reference_dispersion/data/profiles.json'


def intervention_config(base, dwell=150, days=7):
    if dwell != 150 or days not in (7,30,60):raise ValueError('Fixed D150; supported horizons are 7, 30 and 60 days')
    raw=base.to_dict();raw['controller']['dwell_time_s']=float(dwell)
    raw['simulation']['duration_s']=float(days)*86400
    result=config_from_dict(raw)
    expected=base.to_dict();observed=result.to_dict()
    for obj in (expected,observed):
        obj['controller'].pop('dwell_time_s');obj['simulation'].pop('duration_s')
    if expected != observed:raise ValueError('Unexpected intervention input difference')
    return result


def main():
    ap=argparse.ArgumentParser();ap.add_argument('--output',type=Path,default=ROOT/'outputs/reference_a_dwell_20261017/D150_7d')
    ap.add_argument('--days',type=int,choices=[7,30,60],default=7)
    ap.add_argument('--prefix',type=Path)
    args=ap.parse_args()
    if args.output.exists():raise FileExistsError('Use a new output; no duplicate simulation or overwriting')
    spec=yaml.safe_load(SPEC.read_text());base=load_config(SOURCE/'normalized_config.yaml')
    cfg=intervention_config(base,days=args.days)
    if args.days>7 and args.prefix is None:raise ValueError('Extensions require a saved checkpoint; never recalculate the prefix')
    if args.days==7 and args.prefix is not None:raise ValueError('Seven-day experiment starts at epoch')
    frozen=load_frozen(PROFILE);profile=frozen['profiles'][0]
    assert profile['case']=='A' and profile['realized_sigma_a_m']==.5
    assert not np.any(profile['epsilon'])
    meta=json.loads((SOURCE/'metadata.json').read_text())
    assert meta['profiles'][0]==profile and meta['configuration']==base.to_dict()
    env,backend,provider=make_environment(ROOT/spec['environment'],ROOT/spec['weather'],spec['cache_seconds'])
    previous=json.loads((SOURCE/'environment.json').read_text())
    whash=hashlib.sha256((ROOT/spec['weather']).read_bytes()).hexdigest()
    phash=hashlib.sha256(PROFILE.read_bytes()).hexdigest()
    assert previous['environment']==env and previous['weather_sha256']==whash
    assert previous['cache_seconds']==spec['cache_seconds'] and previous['frozen_profiles_sha256']==phash
    initial=np.load(SOURCE/'a_m.npy',mmap_mode='r')[0,0]
    expected=base.derived.reference_radius_m+np.array(profile['delta_a_m'])
    np.testing.assert_array_equal(initial,expected)
    checkpoint=None;prefix_hash=None
    if args.prefix:
        done=json.loads((args.prefix/'completion.json').read_text())
        if not done['complete']:raise ValueError('Incomplete checkpoint')
        old_env=json.loads((args.prefix/'environment.json').read_text())
        for key,value in [('environment',env),('weather_sha256',whash),('frozen_profiles_sha256',phash),('cache_seconds',spec['cache_seconds'])]:
            if old_env[key]!=value:raise ValueError('Continuation environment changed: '+key)
        payload=(args.prefix/'checkpoint.json').read_bytes();prefix_hash=hashlib.sha256(payload).hexdigest()
        checkpoint=json.loads(payload)
        expected_day=7 if args.days==30 else 30
        if checkpoint['time_s']!=expected_day*86400:raise ValueError('Expected preceding 7d or 30d checkpoint')
        state=np.asarray(checkpoint['state']);gaps=cfg.derived.reference_radius_m*np.diff(state[1])
        if not np.isfinite(state).all() or np.min(gaps)<=0 or np.min(state[0])<=cfg.constants.earth_radius_m:
            raise ValueError('Cannot extend an invalid or reversed formation')
    callback=lambda t,end,wall,counts,gap:print(f'{t/86400:g}/{end/86400:g}d; {wall:.1f}s; events={counts}; min gap={gap:.3f}m',flush=True)
    integrate_control_profiles(cfg,[profile],provider,args.output,checkpoint=checkpoint,progress=callback)
    dump_normalized_config(cfg,args.output/'normalized_config.yaml')
    write_json(args.output/'environment.json',dict(environment=env,weather_sha256=whash,
        frozen_profiles_sha256=phash,cache_seconds=spec['cache_seconds'],pymsis_version=__import__('pymsis').__version__,
        msis_version='2.0',density_units='kg/m3',backend_points=backend.points_evaluated,
        D300_reused=str(SOURCE.relative_to(ROOT)),D300_member=0,
        prefix=str(args.prefix) if args.prefix else None,prefix_checkpoint_sha256=prefix_hash,
        maximum_cache_radius_offset_m=provider.maximum_radius_offset_m,
        maximum_cache_phase_offset_rad=provider.maximum_phase_offset_rad))
    write_json(args.output/'input_difference.json',dict(passed=True,
        physical_input_changes={'controller.dwell_time_s':{'D300':300,'D150':150}},
        observation_window={'D300_saved_days':60 if args.days==60 else 30,'compared_days':args.days},
        same_profile=True,same_initial_state=True,same_weather=True,same_environment=True,same_cache_method=True,
        note='The same density field and evaluation method are used; densities need not match after trajectories diverge.'))


if __name__=='__main__':main()
