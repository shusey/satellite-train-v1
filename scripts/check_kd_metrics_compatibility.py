#!/usr/bin/env python3
"""Verify indexed residual postprocessing equals the pre-Kd cached baseline."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from sat_string.config import load_config
from sat_string.io import load_case_arrays,load_switch_events,prepare_output_directory,write_json
from sat_string.long_horizon import infer_episode_gap
from sat_string.residual import residual_analysis


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ["run","all-low","cached-residual","output"]:p.add_argument("--"+name,type=Path,required=True)
    args=p.parse_args();cfg=load_config(args.run/"normalized_config.yaml")
    arrays=load_case_arrays(args.run)
    with np.load(args.all_low/"results.npz") as f:low={key:f[key] for key in ["time_s","a_m"]}
    cached=json.loads((args.cached_residual/"residual_metrics.json").read_text(encoding="utf-8"))
    gap=infer_episode_gap(load_switch_events("outputs/evaluation_20261007/runs/central_14d"))
    report,states=residual_analysis(arrays,cfg,load_switch_events(args.run),gap_s=gap["gap_s"],low_arrays=low,
        relative_D_tolerance=cached["summary"]["relative_D_tolerance"])
    checks={key:report[key]==cached[key] for key in ["local_cycles","group_cycles","quiet_cycles"]}
    checks["summary"]=all(value==cached["summary"][key] for key,value in report["summary"].items())
    with np.load(args.cached_residual/"cycle_end_states.npz") as old:
        checks.update({"state_"+key:bool(np.array_equal(value,old[key])) for key,value in states.items()})
    out=prepare_output_directory(args.output);write_json(out/"compatibility.json",dict(passed=all(checks.values()),checks=checks,
        scope="postprocessing only; no physical simulation",source=str(args.run),cached=str(args.cached_residual)))
    print(dict(passed=all(checks.values()),checks=checks),flush=True)
    return 0 if all(checks.values()) else 1


if __name__=="__main__":raise SystemExit(main())
