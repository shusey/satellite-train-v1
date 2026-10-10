#!/usr/bin/env python3
"""One NRLMSIS case using original RK4/control/actuator primitives."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
from time import perf_counter
import yaml
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from sat_string.config import load_config,config_from_dict
from sat_string.continuation import continue_case
from sat_string.io import prepare_output_directory,load_case_arrays,load_switch_events,load_metadata,save_case,write_json
from sat_string.evaluation import compute_evaluation_metrics
from sat_string.nrlmsis import OrbitGeometry,SpaceWeather,NRLMSIS2,TrajectoryDensityCache
from sat_string.simulator import simulate_case,CaseDefinition
from sat_string.control_cost import clean_json


def make_environment(path,weather,cache_seconds):
    spec=yaml.safe_load(path.read_text(encoding='utf-8'))
    if spec['density_model']!='nrlmsis2' or str(spec['msis_version'])!='2.0' or spec['geomagnetic_activity']!=-1:
        raise ValueError('Explicit NRLMSIS 2.0 and storm-time Ap required')
    geometry=OrbitGeometry(**{k:spec[k] for k in ['epoch_utc','inclination_deg','raan_deg','initial_phase_deg']})
    backend=NRLMSIS2(geometry,SpaceWeather(weather,geometry.epoch_utc))
    mu=load_config('tests/fixtures/regression_24h/config.yaml').constants.mu_m3_s2
    return spec,backend,(TrajectoryDensityCache(backend,cache_seconds,mu) if cache_seconds else backend)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--environment',type=Path,default=Path('configs/nrlmsis2_environment.yaml'))
    p.add_argument('--weather',type=Path,required=True);p.add_argument('--case',choices=['all_low','baseline'],required=True)
    p.add_argument('--duration-s',type=float,default=5184000);p.add_argument('--cache-seconds',type=int,choices=[0,10,30,60,120],required=True)
    p.add_argument('--prefix',type=Path);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();spec,backend,provider=make_environment(args.environment,args.weather,args.cache_seconds)
    raw=load_config('tests/fixtures/regression_24h/config.yaml').to_dict()
    raw['atmosphere']['density_model']='nrlmsis2';raw['disturbance']['enabled']=False;raw['simulation']['duration_s']=args.duration_s
    config=config_from_dict(raw);case=CaseDefinition('nrlmsis2_'+args.case,args.case=='baseline',False)
    out=prepare_output_directory(args.output);start=perf_counter()
    def progress(t,end):print(f'{args.case}: {t/86400:.3f}/{end/86400:g} days; {perf_counter()-start:.1f} seconds',flush=True)
    if args.prefix:
        saved_case=load_metadata(args.prefix)['case']
        if saved_case!=dict(name=case.name,control_enabled=case.control_enabled,disturbance_enabled=case.disturbance_enabled):
            raise ValueError('Cannot change case flags in a continuation')
        provenance=json.loads((args.prefix/'density_environment.json').read_text(encoding='utf-8'))
        if provenance['environment']!=spec or provenance['cache_seconds']!=args.cache_seconds or provenance['weather_sha256']!=hashlib.sha256(args.weather.read_bytes()).hexdigest():
            raise ValueError('Cannot change environment/cache/weather in a continuation')
        source_config=load_config(args.prefix/'normalized_config.yaml')
        result=continue_case(load_case_arrays(args.prefix),source_config,config,case,load_switch_events(args.prefix),
            density_provider=provider,progress_callback=progress)
    else:result=simulate_case(config,case,density_provider=provider,progress_callback=progress)
    save_case(result,config,out)
    metrics=compute_evaluation_metrics(result.arrays(),config,result.switch_events,evaluation_start_s=0)
    write_json(out/'metrics.json',clean_json(metrics))
    write_json(out/'density_environment.json',dict(environment=spec,cache_seconds=args.cache_seconds,
        weather_sha256=hashlib.sha256(args.weather.read_bytes()).hexdigest(),parameter_version='2.0',pymsis_version='0.12.0',
        total_mass_density_units='kg/m3',elapsed_s=perf_counter()-start,backend_points=backend.points_evaluated,
        prefix=str(args.prefix) if args.prefix else None,cache_blocks=getattr(provider,'blocks',None),
        max_cache_radius_offset_m=getattr(provider,'maximum_radius_offset_m',None),
        max_cache_phase_offset_rad=getattr(provider,'maximum_phase_offset_rad',None),
        msis20_parameter_sha256=hashlib.sha256((Path(__import__('pymsis').__file__).parent/'msis2.0.parm').read_bytes()).hexdigest()))
    print('Saved',out,metrics['formation']['max_abs_error_m'],metrics['formation']['rms_error_m'],flush=True)


if __name__=='__main__':main()
