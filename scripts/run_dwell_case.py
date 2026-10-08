#!/usr/bin/env python3
"""Change only dwell from the initial LOW state, reusing a proven pre-switch prefix."""
import argparse
from pathlib import Path
import sys
from time import perf_counter
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from sat_string.config import load_config,config_from_dict
from sat_string.continuation import continue_case
from sat_string.io import load_case_arrays,load_switch_events,load_metadata,prepare_output_directory,save_case,write_json
from sat_string.metrics import compute_case_metrics,save_metrics
from sat_string.simulator import CaseDefinition


def dwell_config(reference,dwell_s,days):
    data=reference.to_dict()
    data["controller"]["dwell_time_s"]=float(dwell_s)
    data["simulation"]["duration_s"]=float(days)*86400
    return config_from_dict(data)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--dwell",type=float,choices=[150.,600.],required=True)
    p.add_argument("--days",type=float,default=30.)
    p.add_argument("--source",type=Path,default=Path("outputs/long_horizon_q_20261008/regression_24h/suite/controlled_disturbed"))
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    source_cfg=load_config(args.source/"normalized_config.yaml")
    frozen=load_config("tests/fixtures/regression_24h/config.yaml")
    if source_cfg!=frozen:raise ValueError("Source must match the frozen baseline configuration")
    case=CaseDefinition(**load_metadata(args.source)["case"])
    if not case.control_enabled or not case.disturbance_enabled:raise ValueError("Controlled disturbed source required")
    events=load_switch_events(args.source)
    cut=min(e["command_time"] for e in events)-source_cfg.controller.control_period_s
    arrays=load_case_arrays(args.source)
    prefix={k:v[arrays["time_s"]<=cut].copy() for k,v in arrays.items()}
    # Initial LOW has next_allowed_time=0. With no commands/completions yet,
    # neither HIGH nor LOW dwell can have influenced ANY prefix state.
    if any(e["command_time"]<=cut for e in events) or np.any(prefix["mode"]!=0) or not np.all(prefix["area_m2"]==source_cfg.spacecraft.area_low_m2):
        raise ValueError("Prefix is not the untouched initial LOW regime")
    config=dwell_config(source_cfg,args.dwell,args.days)
    prefix_cfg=dwell_config(source_cfg,args.dwell,cut/86400)
    output=prepare_output_directory(args.output);start=perf_counter()
    result=continue_case(prefix,prefix_cfg,config,case,[],progress_callback=lambda t,end:print(
        f"dwell={args.dwell:g}: day {t/86400:.2f}/{end/86400:g} ({perf_counter()-start:.1f} s)",flush=True))
    save_case(result,config,output)
    save_metrics(compute_case_metrics(result.arrays(),config,result.switch_events),output)
    write_json(output/"dwell_initialization.json",dict(source_directory=str(args.source.resolve()),
        reused_through_s=cut,configuration_changes=["controller.dwell_time_s","simulation.duration_s"],
        justification="No commands, all initial LOW, next_allowed_time=0; dwell has not acted on the shared prefix",
        elapsed_s=perf_counter()-start))
    print(f"Saved {output}",flush=True)


if __name__=="__main__":main()
