#!/usr/bin/env python3
"""Change only Kd and duration; initial-state simulation uses the existing core."""
import argparse
from pathlib import Path
import sys
from time import perf_counter
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from sat_string.config import load_config,config_from_dict
from sat_string.io import prepare_output_directory,save_case,write_json
from sat_string.metrics import compute_case_metrics,save_metrics
from sat_string.simulator import simulate_case,STANDARD_CASES


def kd_config(reference,multiplier,days):
    data=reference.to_dict()
    data["controller"]["kd_s_inv"]=float(multiplier)*reference.derived.controller_kd_s_inv
    data["simulation"]["duration_s"]=float(days)*86400
    return config_from_dict(data)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--kd-factor",type=float,choices=[.3,3.],required=True)
    p.add_argument("--days",type=float,default=30.)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    frozen=load_config("tests/fixtures/regression_24h/config.yaml")
    config=kd_config(frozen,args.kd_factor,args.days)
    output=prepare_output_directory(args.output);start=perf_counter()
    result=simulate_case(config,STANDARD_CASES[0],progress_callback=lambda t,end:print(
        f"Kd factor {args.kd_factor:g}: day {t/86400:.2f}/{end/86400:g} ({perf_counter()-start:.1f} s)",flush=True))
    save_case(result,config,output)
    save_metrics(compute_case_metrics(result.arrays(),config,result.switch_events),output)
    write_json(output/"kd_initialization.json",dict(kd_factor=args.kd_factor,
        Kp_s2_inv=config.derived.controller_kp_s2_inv,Kd_s_inv=config.derived.controller_kd_s_inv,
        source_configuration="tests/fixtures/regression_24h/config.yaml",
        configuration_changes=["controller.kd_s_inv and derived Kd","simulation.duration_s"],
        prefix_reuse=False,reason="Kd changes trigger times; no controlled pre-switch prefix borrowed",
        elapsed_s=perf_counter()-start))
    print(f"Saved {output}",flush=True)


if __name__=="__main__":main()
