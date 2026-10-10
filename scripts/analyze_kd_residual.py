#!/usr/bin/env python3
"""Reuse saved residual analysis; add physical velocity distribution and P/D gates."""
import argparse
import json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from sat_string.config import load_config
from sat_string.io import load_case_arrays,load_switch_events,load_metadata,prepare_output_directory,write_json
from sat_string.evaluation import compute_evaluation_metrics
from sat_string.long_horizon import infer_episode_gap,group_episodes
from sat_string.residual import residual_analysis
from sat_string.kd_analysis import match_all_low,enhance_residual,threshold_events,return_gates,first_response,fixed_state_forecasts
from sat_string.continuation import restore_attitude
from analyze_dwell_residual import csv_rows


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ["run","all-low","output"]:p.add_argument("--"+name,type=Path,required=True)
    p.add_argument("--cached-residual",type=Path)
    p.add_argument("--reference-14d",type=Path,default=Path("outputs/evaluation_20261007/runs/central_14d"))
    p.add_argument("--dominance-ratio",type=float,default=2.)
    p.add_argument("--relative-D-tolerance",type=float,default=.1)
    args=p.parse_args()
    cfg=load_config(args.run/"normalized_config.yaml");lc=load_config(args.all_low/"normalized_config.yaml")
    if not match_all_low(cfg,lc) or load_metadata(args.all_low)["case"]["control_enabled"]:
        raise ValueError("All-LOW reference must match all physically active inputs")
    if load_metadata(args.run)["case"]["disturbance_enabled"]!=load_metadata(args.all_low)["case"]["disturbance_enabled"]:
        raise ValueError("Disturbance flags differ")
    a=load_case_arrays(args.run);events=load_switch_events(args.run)
    with np.load(args.all_low/"results.npz") as f:low={"time_s":f["time_s"],"a_m":f["a_m"]}
    if not np.array_equal(a["time_s"],low["time_s"]):raise ValueError("Reference grid differs")
    gap=infer_episode_gap(load_switch_events(args.reference_14d))
    if args.cached_residual:
        provenance=json.loads((args.cached_residual/"provenance.json").read_text(encoding="utf-8"))
        if provenance["configuration"]!=cfg.to_dict() or Path(provenance["source"]).resolve()!=args.run.resolve():
            raise ValueError("Cached residual configuration/source differs")
        report=json.loads((args.cached_residual/"residual_metrics.json").read_text(encoding="utf-8"))
        if report["summary"]["observation_end_s"]!=a["time_s"][-1]:raise ValueError("Cached horizon differs")
        if report["summary"]["relative_D_tolerance"]!=args.relative_D_tolerance:
            raise ValueError("Cached forecast tolerance differs; omit cache to change analysis tolerance")
        with np.load(args.cached_residual/"cycle_end_states.npz") as f:states={key:f[key] for key in f.files}
    else:
        report,states=residual_analysis(a,cfg,events,gap_s=gap["gap_s"],low_arrays=low,relative_D_tolerance=args.relative_D_tolerance)
        metrics=compute_evaluation_metrics(a,cfg,events,low_baseline=low)
        report["summary"].update({key:metrics[key] for key in ["formation","switching","altitude","settling"]})
    report=enhance_residual(report,states,cfg)
    triggers=threshold_events(a,cfg,events,args.dominance_ratio);gates=return_gates(a,cfg,events)
    groups=group_episodes(events,float(a["time_s"][-1]),gap["gap_s"])
    group_starts={min(x["command_time"] for x in group if x["to_mode"]=="H") for group in groups}
    for row in triggers:row["temporal_group_onset"]=row["command_time_s"] in group_starts
    first=first_response(a,cfg,events,gap["gap_s"],report,states)
    if cfg.derived.controller_kd_s_inv==2/21600. and report["group_cycles"]:
        endpoint=report["group_cycles"][0];k=endpoint["end_state_row"];time=endpoint["state_sample_time_s"]
        attitude=restore_attitude(cfg,events,time)
        first["fixed_baseline_state_forecasts"]=fixed_state_forecasts(states["gap_error_m"][k],
            states["gap_rate_m_s"][k],time,attitude.next_allowed_time,cfg)
        first["fixed_state_forecast_definition"]="same first-group baseline end-state and clocks; no changed-state simulation"
    report["summary"].update(kd_factor=cfg.derived.controller_kd_s_inv/(2/21600.),
        dominance_ratio=args.dominance_ratio,
        command_contribution_counts={kind:sum(x["command_contribution_class"]==kind for x in triggers)
            for kind in ["P_dominant","D_dominant","mixed","zero"]},
        crossing_contribution_counts={kind:sum(x["crossing_contribution_class"]==kind for x in triggers)
            for kind in ["P_dominant","D_dominant","mixed","zero"]},
        quiet_onset_command_contribution_counts={kind:sum(x["command_contribution_class"]==kind and x["quiet_packet_onset"] for x in triggers)
            for kind in ["P_dominant","D_dominant","mixed","zero"]},
        later_in_packet_command_contribution_counts={kind:sum(x["command_contribution_class"]==kind and not x["quiet_packet_onset"] for x in triggers)
            for kind in ["P_dominant","D_dominant","mixed","zero"]},
        temporal_group_onset_contribution_counts={kind:sum(x["command_contribution_class"]==kind and x["temporal_group_onset"] for x in triggers)
            for kind in ["P_dominant","D_dominant","mixed","zero"]},
        first_off_gate_blocked_cycles=sum(x["off_not_met_at_gate"] is True for x in gates),
        gate_blocked_completed_cycles=sum(x["off_not_met_at_gate"] is True and x["completed_by_tmax"] for x in gates),
        return_gate_censored_cycles=sum(x["return_right_censored"] for x in gates))
    report["source_directory"]=str(args.run.resolve());report["state_array_file"]="cycle_end_states.npz"
    json.dumps(report,allow_nan=False);json.dumps(first,allow_nan=False)
    out=prepare_output_directory(args.output)
    write_json(out/"residual_metrics.json",report);write_json(out/"summary.json",report["summary"])
    for key,name in [("local_cycles","local_cycles.csv"),("group_cycles","group_residuals.csv"),("quiet_cycles","quiet_residuals.csv")]:
        csv_rows(out/name,report[key])
    csv_rows(out/"threshold_events.csv",triggers);csv_rows(out/"return_gates.csv",gates)
    csv_rows(out/"first_local_cycles.csv",first.get("first_local_cycles",[]))
    write_json(out/"first_response.json",first)
    np.savez_compressed(out/"cycle_end_states.npz",**states)
    write_json(out/"provenance.json",dict(source=str(args.run.resolve()),all_low=str(args.all_low.resolve()),
        configuration=cfg.to_dict(),cached_residual=str(args.cached_residual) if args.cached_residual else None,
        grouping=gap,dominance_ratio=args.dominance_ratio,all_low_reuse_reason="Kd cannot change orbit/area when controller disabled"))
    s=report["summary"];print(dict(factor=s["kd_factor"],days=s["observation_end_s"]/86400,
        e_max=s["formation"]["max_abs_error_m"],rms=s["formation"]["rms_error_m"],
        cycles=s["local_completed_cycles"],transitions=s["switching"]["n_transition_starts"],
        extra_loss=s["altitude"]["extra_loss_vs_all_low_m"],blocked_gates=s["first_off_gate_blocked_cycles"],
        HIGH_hold_s=s["high_hold_times_s"],contribution=s["command_contribution_counts"]),flush=True)


if __name__=="__main__":main()
