#!/usr/bin/env python3
"""Separate true sampled recrossings, LOW entry above-on, and dwell-limited repeats."""
import argparse
from collections import Counter
import csv
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from sat_string.config import load_config
from sat_string.io import prepare_output_directory,write_json


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--analysis",type=Path,nargs="+",required=True);p.add_argument("--output",type=Path,required=True)
    args=p.parse_args();results=[]
    for directory in args.analysis:
        s=json.loads((directory/"summary.json").read_text(encoding="utf-8"))
        provenance=json.loads((directory/"provenance.json").read_text(encoding="utf-8"));cfg=load_config(Path(provenance["source"])/"normalized_config.yaml")
        kinds=Counter();upward=Counter();delayed=0;opposed=0
        with (directory/"threshold_events.csv").open(encoding="utf-8-sig") as f:
            for row in csv.DictReader(f):
                kinds[row["crossing_kind"]]+=1
                if row["crossing_kind"]=="upward_from_below":upward[row["crossing_contribution_class"]]+=1
                delayed+=float(row["command_minus_crossing_s"])>0
                opposed+=row["crossing_opposed_signs"]=="True"
        bound=2*(cfg.controller.dwell_time_s+cfg.controller.slew_time_s+cfg.controller.command_delay_s)
        at_bound=0;observed=0
        with (directory/"local_cycles.csv").open(encoding="utf-8-sig") as f:
            for row in csv.DictReader(f):
                if row["next_same_satellite_start_interval_s"]:
                    observed+=1;at_bound+=abs(float(row["next_same_satellite_start_interval_s"])-bound)<1e-9
        results.append(dict(kd_factor=s["kd_factor"],days=s["observation_end_s"]/86400,
            crossing_kinds=dict(kinds),upward_crossing_contribution_counts=dict(upward),
            command_later_than_final_crossing_sample=delayed,opposed_P_D_sign_crossing_states=opposed,
            ideal_minimum_LH_repeat_interval_s=bound,observed_uncensored_LH_intervals=observed,
            LH_intervals_at_dwell_slew_bound=at_bound,
            bound_interpretation="two slew + two dwell + two delays; control grid can lengthen it; current delays zero"))
    out=prepare_output_directory(args.output);write_json(out/"trigger_gates.json",results)
    print(json.dumps(results,indent=2),flush=True)


if __name__=="__main__":main()
