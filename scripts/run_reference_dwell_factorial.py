"""Run only the two missing A-ref factorial cells, for exactly 30 days."""
import argparse
import hashlib
import json
import sys
from pathlib import Path
import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'scripts'))
sys.path.insert(0, str(ROOT/'src'))
from run_reference_a_dwell import SOURCE, PROFILE, SPEC
from run_nrlmsis_case import make_environment
from sat_string.config import load_config, dump_normalized_config
from sat_string.reference_condition import load_frozen
from sat_string.dispersion_control import integrate_control_profiles
from sat_string.dwell_factorial import factorial_config
from sat_string.io import write_json

OUTPUT = ROOT/'outputs/reference_dwell_factorial_20261018'


def source_manifest():
    folders=[SOURCE, ROOT/'outputs/reference_a_dwell_20261017/D150_7d',
             ROOT/'outputs/reference_a_dwell_60d_20261017/D150_30d']
    paths=[p for folder in folders for p in folder.rglob('*') if p.is_file()]
    paths += [PROFILE,ROOT/'docs/real_environment_applicability/data/active_drivers.csv']
    return {str(p.relative_to(ROOT)):dict(size=p.stat().st_size,mtime_ns=p.stat().st_mtime_ns,
        sha256=hashlib.sha256(p.read_bytes()).hexdigest() if p.stat().st_size<20_000_000 else None)
        for p in paths}


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--case', choices=['H150-L300', 'H300-L150'], required=True)
    ap.add_argument('--root', type=Path, default=OUTPUT)
    args = ap.parse_args()
    folder = args.root/args.case
    if folder.exists():
        raise FileExistsError('Never overwrite or repeat an existing run')
    args.root.mkdir(parents=True,exist_ok=True)
    manifest=args.root/'source_manifest_before.json'
    if not manifest.exists():write_json(manifest,source_manifest())
    high, low = (150, 300) if args.case == 'H150-L300' else (300, 150)
    base = load_config(SOURCE/'normalized_config.yaml')
    cfg = factorial_config(base, high, low)
    profile = load_frozen(PROFILE)['profiles'][0]
    saved = json.loads((SOURCE/'metadata.json').read_text())
    assert saved['profiles'][0] == profile and saved['configuration'] == base.to_dict()
    assert not np.any(profile['epsilon'])
    assert np.std(profile['delta_a_m'], ddof=1) == .5
    np.testing.assert_array_equal(np.load(SOURCE/'a_m.npy', mmap_mode='r')[0, 0],
                                  cfg.derived.reference_radius_m+np.array(profile['delta_a_m']))
    spec = yaml.safe_load(SPEC.read_text())
    env, backend, provider = make_environment(ROOT/spec['environment'], ROOT/spec['weather'], spec['cache_seconds'])
    previous = json.loads((SOURCE/'environment.json').read_text())
    hashes = dict(weather_sha256=hashlib.sha256((ROOT/spec['weather']).read_bytes()).hexdigest(),
                  frozen_profiles_sha256=hashlib.sha256(PROFILE.read_bytes()).hexdigest())
    assert previous['environment'] == env and previous['cache_seconds'] == spec['cache_seconds']
    assert all(previous[k] == v for k, v in hashes.items())
    callback = lambda t,end,wall,counts,gap: print(
        f'{args.case} {t/86400:g}/{end/86400:g}d; {wall:.1f}s; events={counts}; min gap={gap:.3f}m', flush=True)
    integrate_control_profiles(cfg, [profile], provider, folder, progress=callback)
    dump_normalized_config(cfg, folder/'normalized_config.yaml')
    write_json(folder/'environment.json', dict(environment=env, **hashes,
        cache_seconds=spec['cache_seconds'], pymsis_version=__import__('pymsis').__version__,
        msis_version='2.0', density_units='kg/m3', backend_points=backend.points_evaluated,
        maximum_cache_radius_offset_m=provider.maximum_radius_offset_m,
        maximum_cache_phase_offset_rad=provider.maximum_phase_offset_rad))
    write_json(folder/'input_difference.json', dict(passed=True,
        changes={'controller.dwell_high_time_s':high, 'controller.dwell_low_time_s':low},
        same_initial_state=True, same_profile=True, same_environment=True, same_weather=True,
        same_cache_method=True, compared_days=30, no_automatic_extension=True))


if __name__ == '__main__':
    main()
