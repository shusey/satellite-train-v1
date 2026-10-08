#!/usr/bin/env python3
"""Add residual metrics to saved runs; reuse old aggregate metrics/grouping when supplied."""
import argparse
import csv
import json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from sat_string.config import load_config
from sat_string.evaluation import compute_evaluation_metrics
from sat_string.io import load_case_arrays,load_switch_events,load_metadata,prepare_output_directory,write_json
from sat_string.long_horizon import infer_episode_gap
from sat_string.residual import residual_analysis


def csv_rows(path,rows):
    if not rows:return
    flat=[]
    for row in rows:
        item=dict(row);prediction=item.pop("prediction",{})
        item.update({"forecast_"+k:v for k,v in prediction.items()})
        flat.append({k:json.dumps(v,ensure_ascii=False) if isinstance(v,(list,dict)) else v for k,v in item.items()})
    keys=list(dict.fromkeys(k for row in flat for k in row))
    with path.open("w",encoding="utf-8-sig",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=keys);writer.writeheader();writer.writerows(flat)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ["run","all-low","output"]:p.add_argument("--"+name,type=Path,required=True)
    p.add_argument("--cached-analysis",type=Path,help="reuse existing metrics.json and waves.json, without reanalysis")
    p.add_argument("--reference-14d",type=Path,default=Path("outputs/evaluation_20261007/runs/central_14d"))
    p.add_argument("--relative-D-tolerance",type=float,default=.1)
    args=p.parse_args()
    cfg=load_config(args.run/"normalized_config.yaml");lc=load_config(args.all_low/"normalized_config.yaml")
    left,right=cfg.to_dict(),lc.to_dict()
    for data in [left,right]:data["controller"].pop("dwell_time_s")
    if left!=right or load_metadata(args.all_low)["case"]["control_enabled"]:
        raise ValueError("All-LOW reference differs in a physically active input")
    if load_metadata(args.run)["case"]["disturbance_enabled"]!=load_metadata(args.all_low)["case"]["disturbance_enabled"]:
        raise ValueError("Disturbance flags differ")
    arrays=load_case_arrays(args.run);events=load_switch_events(args.run)
    with np.load(args.all_low/"results.npz") as f:low={"time_s":f["time_s"],"a_m":f["a_m"]}
    if not np.array_equal(arrays["time_s"],low["time_s"]):raise ValueError("Reference time grid differs")
    gap=infer_episode_gap(load_switch_events(args.reference_14d))
    groups=None
    if args.cached_analysis:
        metrics=json.loads((args.cached_analysis/"metrics.json").read_text(encoding="utf-8"))
        cached=json.loads((args.cached_analysis/"waves.json").read_text(encoding="utf-8"))
        if cached["observation_end_s"]!=arrays["time_s"][-1] or cfg.controller.dwell_time_s!=300:
            raise ValueError("Cached baseline analysis endpoint/configuration mismatch")
        groups=cached["waves"]
    else:metrics=compute_evaluation_metrics(arrays,cfg,events,low_baseline=low)
    report,states=residual_analysis(arrays,cfg,events,groups=groups,gap_s=gap["gap_s"],
        relative_D_tolerance=args.relative_D_tolerance,low_arrays=low)
    report["source_directory"]=str(args.run.resolve())
    report["state_array_file"]="cycle_end_states.npz"
    report["summary"].update(formation=metrics["formation"],switching=metrics["switching"],altitude=metrics["altitude"],settling=metrics["settling"])
    json.dumps(report,allow_nan=False)
    output=prepare_output_directory(args.output)
    write_json(output/"residual_metrics.json",report)
    write_json(output/"summary.json",report["summary"])
    csv_rows(output/"local_cycles.csv",report["local_cycles"])
    csv_rows(output/"group_residuals.csv",report["group_cycles"])
    csv_rows(output/"quiet_residuals.csv",report["quiet_cycles"])
    np.savez_compressed(output/"cycle_end_states.npz",**states)
    write_json(output/"provenance.json",dict(source=str(args.run.resolve()),all_low=str(args.all_low.resolve()),
        all_low_dwell_reuse_reason="Controller disabled; dwell never affects all-LOW area or dynamics",
        cached_analysis=str(args.cached_analysis.resolve()) if args.cached_analysis else None,
        configuration=cfg.to_dict(),grouping=gap,relative_D_tolerance=args.relative_D_tolerance))
    print(dict(dwell=cfg.controller.dwell_time_s,days=arrays["time_s"][-1]/86400,
        e_max=metrics["formation"]["max_abs_error_m"],rms=metrics["formation"]["rms_error_m"],
        switches=metrics["switching"]["n_transition_starts"],cycles=report["summary"]["local_completed_cycles"],
        global_endpoints=len(report["group_cycles"]),valid_group_predictions=report["summary"]["group_valid_finite_forecasts"],
        extra_loss=metrics["altitude"]["extra_loss_vs_all_low_m"],high_hold_s=report["summary"]["high_hold_times_s"]),flush=True)


if __name__=="__main__":main()
