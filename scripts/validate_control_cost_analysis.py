#!/usr/bin/env python3
"""Reproduction and step-refinement checks using the unchanged simulation core."""
import argparse
import json
from pathlib import Path
import re
import subprocess
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
import numpy as np
from sat_string.config import config_from_dict,load_config
from sat_string.control_cost import evaluate_cost
from sat_string.io import load_case_arrays,load_switch_events,save_case,write_json
from sat_string.simulator import STANDARD_CASES,simulate_case


def differences(left,right):
    return {key:dict(bitwise_equal=bool(np.array_equal(left[key],right[key])),
                      max_abs_difference=float(np.max(np.abs(left[key].astype(float)-right[key].astype(float)))))
            for key in left}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--before",type=Path,required=True)
    p.add_argument("--legacy",type=Path,default=Path("runs/experiment_13_1.5"))
    p.add_argument("--refine",action="append",default=[])
    args=p.parse_args()
    out=args.output
    report={}
    core=["config.py","controller.py","dynamics.py","simulator.py","integrator.py",
          "atmosphere.py","attitude.py","disturbance.py"]
    report["unchanged_core_vs_source_commit"]={
        name:subprocess.check_output(["git","diff","9449a78","--",f"src/sat_string/{name}"],text=True)==""
        for name in core}
    report["before_after"]={}
    report["saved_legacy"]={}
    for case in STANDARD_CASES:
        name=case.name
        path=out/"baseline"/name
        report["before_after"][name]=differences(load_case_arrays(args.before/name),load_case_arrays(path))
        report["before_after"][name]["switch_events_identical"]=load_switch_events(args.before/name)==load_switch_events(path)
        old=json.loads((args.legacy/name/"metrics.json").read_text(encoding="utf-8-sig"))
        new=json.loads((path/"metrics.json").read_text(encoding="utf-8-sig"))
        log=(args.legacy/name/"run.log").read_text(encoding="utf-8-sig")
        starts=re.findall(r"t=([\d.e+-]+) s satellite=(\d+) command (L|H)->(L|H); start=([\d.e+-]+) complete=([\d.e+-]+) s",log)
        events=[dict(command_time=float(t),satellite_id=int(i),from_mode=f,to_mode=to,start_time=float(s),complete_time=float(c))
                for t,i,f,to,s,c in starts]
        report["saved_legacy"][name]=dict(all_metrics_equal=old==new,event_log_equal=events==load_switch_events(path))
    report["refinement"]=json.loads((out/"validation.json").read_text(encoding="utf-8-sig")).get("refinement",{}) if (out/"validation.json").exists() else {}
    targets=args.refine or ["baseline/controlled_disturbed"]
    low=load_case_arrays(out/"baseline"/"all_low_drag_baseline")
    for relative in targets:
        path=out/relative
        cfg=load_config(path/"normalized_config.yaml")
        coarse=load_case_arrays(path)
        data=cfg.to_dict()
        # Repository standard_suite already uses 1 s, used here only as numerical refinement.
        data["simulation"]["integration_step_s"]=cfg.simulation.integration_step_s/2
        finecfg=config_from_dict(data)
        destination=out/"validation_runs"/(path.name+"_half_step")
        destination.mkdir(parents=True,exist_ok=False)
        result=simulate_case(finecfg,STANDARD_CASES[0])
        save_case(result,finecfg,destination)
        fine=result.arrays()
        fm=evaluate_cost(fine,finecfg,result.switch_events,low)
        cm=evaluate_cost(coarse,cfg,load_switch_events(path),low)
        item=dict(integration_steps_s=[cfg.simulation.integration_step_s,finecfg.simulation.integration_step_s],
                  arrays=differences(coarse,fine),
                  switch_events_identical=load_switch_events(path)==load_switch_events(destination),
                  delta_e_max_m=fm["maximum_gap_error_m"]-cm["maximum_gap_error_m"],
                  delta_extra_mean_a_loss_m=fm["extra_mean_a_loss_vs_low_m"]-cm["extra_mean_a_loss_vs_low_m"],
                  delta_switch_count=fm["activity"]["switch_count"]-cm["activity"]["switch_count"],
                  delta_high_equivalent_s=fm["activity"]["high_equivalent_total_s"]-cm["activity"]["high_equivalent_total_s"])
        report["refinement"][relative]=item
        print(relative,json.dumps({k:v for k,v in item.items() if k!="arrays"}),flush=True)
    write_json(out/"validation.json",report)


if __name__=="__main__":
    main()


