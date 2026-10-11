"""Extend ONLY the reused A/B all-LOW30d checkpoints; no prefix recalculation."""
import hashlib,json
from pathlib import Path
import numpy as np
from run_dispersion_control import ROOT,OUT,DATA,configuration,frozen_profiles
from run_nrlmsis_case import make_environment
from sat_string.passive_dispersion import integrate_profiles
from sat_string.config import dump_normalized_config
from sat_string.io import write_json


def main():
    spec,cfg=configuration(60)
    decision=json.loads((DATA/'extension_decision.json').read_text())
    if not decision['extend_to_60d']:raise ValueError('No go decision')
    selected=frozen_profiles(spec)[:2]
    original=ROOT/spec['screening_output']/'extended_30d'
    meta=json.loads((original/'metadata.json').read_text());old=meta['profiles']
    ix=[next(i for i,p in enumerate(old) if p['id']==q['source_profile']) for q in selected]
    profiles=[old[i] for i in ix]
    a=np.load(original/'a_m.npy',mmap_mode='r');u=np.load(original/'lambda_rad.npy',mmap_mode='r')
    state=np.vstack((a[-1,ix].ravel(),u[-1,ix].ravel()))
    env,backend,provider=make_environment(ROOT/spec['environment'],ROOT/spec['weather'],spec['cache_seconds'])
    previous=json.loads((original/'environment.json').read_text())
    weather_hash=hashlib.sha256((ROOT/spec['weather']).read_bytes()).hexdigest()
    if env!=previous['spec'] or weather_hash!=previous['weather_sha256']:raise ValueError('Frozen OFF environment mismatch')
    out=OUT/'OFF_A_B_60d'
    callback=lambda t,end,wall:print(f'OFF A/B suffix: {t/86400:g}/{end/86400:g}d; {wall:.1f}s',flush=True)
    integrate_profiles(cfg,profiles,provider,out,start_s=30*86400,initial_state=state,progress=callback)
    dump_normalized_config(cfg,out/'normalized_config.yaml')
    write_json(out/'environment.json',dict(environment=env,weather_sha256=weather_hash,cache_seconds=spec['cache_seconds'],
        pymsis_version='0.12.0',msis_version='2.0',density_units='kg/m3',backend_points=backend.points_evaluated,
        maximum_cache_radius_offset_m=provider.maximum_radius_offset_m,maximum_cache_phase_offset_rad=provider.maximum_phase_offset_rad,
        prefix=str(original)))


if __name__=='__main__':main()
