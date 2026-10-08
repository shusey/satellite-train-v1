#!/usr/bin/env python3
"""Postprocess one long run: waves, P/D triggers, settling and cost."""
import argparse
import csv
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
import numpy as np
from sat_string.config import load_config
from sat_string.evaluation import EvaluationSettings,compute_evaluation_metrics,settling_metrics
from sat_string.io import load_case_arrays,load_metadata,load_switch_events,prepare_output_directory,write_json
from sat_string.long_horizon import analyze_long_horizon,infer_episode_gap


def write_csv(path,rows):
    if not rows:
        return
    fields=list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w",encoding="utf-8-sig",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=fields);writer.writeheader();writer.writerows(rows)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--run",type=Path,required=True)
    p.add_argument("--all-low",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--reference-14d",type=Path,default=Path("outputs/evaluation_20261007/runs/central_14d"))
    args=p.parse_args()
    cfg=load_config(args.run/"normalized_config.yaml")
    low_cfg=load_config(args.all_low/"normalized_config.yaml")
    if cfg!=low_cfg or load_metadata(args.all_low)["case"]["control_enabled"]:
        raise ValueError("All-LOW reference must match the entire input configuration")
    if load_metadata(args.run)["case"]["disturbance_enabled"]!=load_metadata(args.all_low)["case"]["disturbance_enabled"]:
        raise ValueError("Disturbance flags differ")
    arrays=load_case_arrays(args.run);events=load_switch_events(args.run)
    with np.load(args.all_low/"results.npz",allow_pickle=False) as archive:
        low={"time_s":archive["time_s"],"a_m":archive["a_m"]}
    gap=infer_episode_gap(load_switch_events(args.reference_14d))
    report,derived=analyze_long_horizon(arrays,cfg,events,low,gap["gap_s"])
    metrics=compute_evaluation_metrics(arrays,cfg,events,low_baseline=low)
    sensitivity=[]
    for eps in [1.,5.,10.,50.]:
        settings=EvaluationSettings(eps,cfg.recovery.velocity_tolerance_m_s,cfg.recovery.hold_time_s)
        sensitivity.append(settling_metrics(arrays["time_s"],arrays["gap_error_m"],arrays["gap_rate_m_s"],
            metrics["switching"],settings,evaluation_start_s=cfg.disturbance.start_time_s+cfg.disturbance.duration_s))
    output=prepare_output_directory(args.output)
    write_json(output/"metrics.json",metrics)
    write_json(output/"waves.json",report)
    write_json(output/"settling_sensitivity.json",{"settings_results":sensitivity})
    write_json(output/"provenance.json",dict(source_directory=str(args.run.resolve()),all_low_directory=str(args.all_low.resolve()),
               gap_reference=str(args.reference_14d.resolve()),episode_gap=gap,configuration=cfg.to_dict()))
    write_csv(output/"wave_summary.csv",[{k:v for k,v in row.items() if k not in ["outer_branch_fits"]} for row in report["waves"]])
    write_csv(output/"q_triggers.csv",report["triggers"])
    write_csv(output/"end_state_differences.csv",report["cycle_end_differences"])
    write_csv(output/"local_cycle_intervals.csv",report["local_cycle_intervals"])
    write_csv(output/"seven_day_blocks.csv",report["seven_day_blocks"])
    keep=(arrays["time_s"]>=7*86400)&(arrays["time_s"]<=10*86400)
    np.savez_compressed(output/"detail_7_to_10_days.npz",
        **{key:arrays[key][keep] for key in ["time_s","gap_error_m","gap_rate_m_s","q_m_s2","mode","area_m2","a_m"]},
        q_P_m_s2=derived["q_P_m_s2"][keep],q_D_m_s2=derived["q_D_m_s2"][keep])
    print(dict(days=arrays["time_s"][-1]/86400,e_max=metrics["formation"]["max_abs_error_m"],
        rms=metrics["formation"]["rms_error_m"],switches=metrics["switching"]["n_transition_starts"],
        last_switch_day=(metrics["switching"]["last_switch_time_s"] or 0)/86400,
        extra_loss_m=metrics["altitude"]["extra_loss_vs_all_low_m"],
        waves=[(w["wave_number"],w["start_time_s"]/86400,w["transition_starts"],w["E_n_active_max_error_m"],w["complete"]) for w in report["waves"]]),flush=True)


if __name__=="__main__":
    main()
