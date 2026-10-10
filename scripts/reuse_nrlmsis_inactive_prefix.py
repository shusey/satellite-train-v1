#!/usr/bin/env python3
"""Prepare controlled case from proven inactive all-LOW history; never infer future inactivity."""
import argparse
from dataclasses import replace
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from sat_string.config import load_config
from sat_string.io import load_case_arrays,load_switch_events,load_metadata,save_case,prepare_output_directory,write_json
from sat_string.simulator import SimulationResult,CaseDefinition
from sat_string.evaluation import compute_evaluation_metrics
from sat_string.nrlmsis_analysis import inactive_prefix_proof
from sat_string.control_cost import clean_json


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--all-low',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();config=load_config(args.all_low/'normalized_config.yaml')
    meta=load_metadata(args.all_low)
    if meta['case']['control_enabled'] or meta['case']['disturbance_enabled'] or config.atmosphere.density_model!='nrlmsis2':
        raise ValueError('Expected uncontaminated NRLMSIS all-LOW reference')
    if load_switch_events(args.all_low):raise ValueError('All-LOW reference contains transitions')
    arrays=load_case_arrays(args.all_low)
    env=json.loads((args.all_low/'density_environment.json').read_text(encoding='utf-8'))
    proof=inactive_prefix_proof(arrays,config,env['cache_seconds'])
    if proof['reusable_until_s']<=0:raise ValueError('No reusable prefix; run controlled case directly')
    prefix={k:v[arrays['time_s']<=proof['reusable_until_s']] for k,v in arrays.items()}
    config=replace(config,simulation=replace(config.simulation,duration_s=proof['reusable_until_s']))
    result=SimulationResult(case=CaseDefinition('nrlmsis2_baseline',True,False),switch_events=[],**prefix)
    out=prepare_output_directory(args.output);save_case(result,config,out)
    env.update(inactive_prefix_source=str(args.all_low),computation_method='all-LOW trajectory proven controller-inactive at every tick')
    write_json(out/'density_environment.json',env);write_json(out/'inactivity_proof.json',proof)
    write_json(out/'metrics.json',clean_json(compute_evaluation_metrics(prefix,config,evaluation_start_s=0)))
    print(proof,flush=True)
    print('Complete equivalent case' if proof['first_potential_high_time_s'] is None else 'Continue real controlled case from this prefix',flush=True)


if __name__=='__main__':main()
