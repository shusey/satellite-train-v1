"""Extend the two completed cross cells from day 30 to day 60, never rerun a prefix."""
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import yaml
from run_reference_dwell_factorial import ROOT,OUTPUT as PREFIX_ROOT,SOURCE,PROFILE,SPEC,source_manifest
from run_nrlmsis_case import make_environment
from sat_string.config import load_config,dump_normalized_config
from sat_string.dwell_factorial import extension_config
from sat_string.reference_condition import load_frozen
from sat_string.dispersion_control import integrate_control_profiles
from sat_string.continuation import restore_attitude
from sat_string.io import write_json

OUTPUT=ROOT/'outputs/reference_dwell_factorial_60d_20261019'


def extended_manifest():
    result=source_manifest()
    folders=[PREFIX_ROOT/'H150-L300',PREFIX_ROOT/'H300-L150',SOURCE.parent/'batch_60d',
             ROOT/'outputs/reference_a_dwell_60d_20261017/D150_60d']
    for folder in folders:
        for p in folder.rglob('*'):
            if p.is_file():
                st=p.stat();result[str(p.relative_to(ROOT))]=dict(size=st.st_size,mtime_ns=st.st_mtime_ns,
                    sha256=hashlib.sha256(p.read_bytes()).hexdigest() if st.st_size<20_000_000 else None)
    return result


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--case',choices=['H150-L300','H300-L150'],required=True)
    args=ap.parse_args();prefix=PREFIX_ROOT/args.case;out=OUTPUT/args.case
    if out.exists():raise FileExistsError('Never overwrite or duplicate a saved suffix')
    OUTPUT.mkdir(exist_ok=True)
    if not (OUTPUT/'source_manifest_before.json').exists():write_json(OUTPUT/'source_manifest_before.json',extended_manifest())
    done=json.loads((prefix/'completion.json').read_text());assert done['complete'] and done['end_s']==30*86400.
    payload=(prefix/'checkpoint.json').read_bytes();checkpoint=json.loads(payload)
    profile=load_frozen(PROFILE)['profiles'][0];cfg=extension_config(checkpoint,profile)
    assert load_config(prefix/'normalized_config.yaml').to_dict()==checkpoint['configuration']
    assert json.loads((prefix/'metadata.json').read_text())['profiles']==[profile]
    state=np.asarray(checkpoint['state'])
    for k,key in enumerate(['a_m','lambda_rad']):
        np.testing.assert_array_equal(state[k],np.load(prefix/f'{key}.npy',mmap_mode='r')[-1,0])
    restored=restore_attitude(cfg,checkpoint['events'][0],checkpoint['time_s'])
    np.testing.assert_array_equal(restored.mode,np.load(prefix/'mode.npy',mmap_mode='r')[-1,0])
    np.testing.assert_array_equal(restored.areas(checkpoint['time_s']),np.load(prefix/'area_m2.npy',mmap_mode='r')[-1,0])
    spec=yaml.safe_load(SPEC.read_text())
    env,backend,provider=make_environment(ROOT/spec['environment'],ROOT/spec['weather'],spec['cache_seconds'])
    previous=json.loads((prefix/'environment.json').read_text())
    whash=hashlib.sha256((ROOT/spec['weather']).read_bytes()).hexdigest()
    phash=hashlib.sha256(PROFILE.read_bytes()).hexdigest()
    assert previous['environment']==env and previous['weather_sha256']==whash
    assert previous['frozen_profiles_sha256']==phash and previous['cache_seconds']==spec['cache_seconds']
    assert previous['pymsis_version']==__import__('pymsis').__version__ and previous['msis_version']=='2.0'
    callback=lambda t,end,wall,counts,gap:print(f'{args.case} {t/86400:g}/{end/86400:g}d; {wall:.1f}s; events={counts}; min gap={gap:.3f}m',flush=True)
    integrate_control_profiles(cfg,[profile],provider,out,checkpoint=checkpoint,progress=callback)
    dump_normalized_config(cfg,out/'normalized_config.yaml')
    write_json(out/'environment.json',dict(environment=env,weather_sha256=whash,frozen_profiles_sha256=phash,
        cache_seconds=spec['cache_seconds'],pymsis_version=__import__('pymsis').__version__,msis_version='2.0',
        density_units='kg/m3',backend_points=backend.points_evaluated,
        prefix=str(prefix.relative_to(ROOT)),prefix_checkpoint_sha256=hashlib.sha256(payload).hexdigest(),
        maximum_cache_radius_offset_m=provider.maximum_radius_offset_m,
        maximum_cache_phase_offset_rad=provider.maximum_phase_offset_rad))
    write_json(out/'input_difference.json',dict(passed=True,changed_only='simulation.duration_s: 30d -> 60d',
        exact_checkpoint_state=True,exact_checkpoint_actuator=True,same_profile=True,same_environment=True,
        same_weather=True,same_density_method=True,prefix_reused=True))


if __name__=='__main__':main()
