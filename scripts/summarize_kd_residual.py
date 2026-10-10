#!/usr/bin/env python3
"""Compare saved Kd reports; no simulation or duplicate endpoint analysis."""
import argparse
import csv
import json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from sat_string.io import prepare_output_directory,write_json
from sat_string.io import load_switch_events
from sat_string.config import load_config
from sat_string.kd_analysis import first_day_metrics,velocity_distribution
from sat_string.long_horizon import decompose_q
from sat_string.residual import predict_crossing
from sat_string.continuation import restore_attitude
from summarize_dwell_residual import summarize,stats


def summarize_kd(report):
    row=summarize(report);s=report["summary"]
    row.update(kd_factor=s["kd_factor"],Kd_s_inv=s["Kd_s_inv"],Kp_s2_inv=s["Kp_s2_inv"],
        first_off_gate_blocked_cycles=s["first_off_gate_blocked_cycles"],
        gate_blocked_completed_cycles=s["gate_blocked_completed_cycles"],
        return_gate_censored_cycles=s["return_gate_censored_cycles"])
    for scope,key in [("local","local_cycles"),("quiet","quiet_cycles"),("group","group_cycles")]:
        for metric in ["V_res_all_m_s","V_res_rms_all_m_s","V_res_median_abs_all_m_s",
                "worsening_gap_fraction","worsening_weight_fraction","error_squared_mean_rate_m2_s"]:
            row[scope+"_"+metric]=stats([x[metric] for x in report[key]])
    for key in ["command_contribution_counts","crossing_contribution_counts",
            "quiet_onset_command_contribution_counts","later_in_packet_command_contribution_counts","temporal_group_onset_contribution_counts"]:
        row[key]=s[key]
    row["first_group"]=report["group_cycles"][0] if report["group_cycles"] else None
    row["last_quiet_endpoint"]=report["quiet_cycles"][-1] if report["quiet_cycles"] else None
    valid=[x for x in report["quiet_cycles"] if not x["right_censored"] and x["prediction_valid_on_observed_window"] and x["prediction"]["predicted_wait_s"] is not None]
    row["quiet_valid_forecast_status_counts"]={status:sum(x["prediction"]["status"]==status for x in valid)
        for status in ["future_upward_crossing","already_at_or_above_on"]}
    future=[x for x in valid if x["prediction"]["status"]=="future_upward_crossing"]
    row["quiet_future_upward_relative_error"]=stats([abs(x["relative_prediction_error"]) for x in future])
    row["quiet_future_upward_wait_day"]=stats([x["actual_next_wait_s"]/86400 for x in future])
    return row


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--analysis",type=Path,nargs="+",required=True)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--gain-axis",choices=["kd","kp"],default="kd")
    args=p.parse_args();rows=[]
    for directory in args.analysis:
        report=json.loads((directory/"residual_metrics.json").read_text(encoding="utf-8"))
        row=summarize_kd(report)
        if args.gain_axis=="kp":
            row["kp_factor"]=report["summary"]["kp_factor"]
            for key in ["off_before_gate_observed_cycles","off_before_gate_met_cycles"]:
                row[key]=report["summary"][key]
            with (directory/"threshold_events.csv").open(encoding="utf-8-sig") as f:triggers=list(csv.DictReader(f))
            row["activation_statistics"]={}
            for role,predicate in [("all",lambda x:True),("quiet_onset",lambda x:x["quiet_packet_onset"]=="True"),
                    ("temporal_group_onset",lambda x:x["temporal_group_onset"]=="True")]:
                selected=[x for x in triggers if predicate(x)]
                row["activation_statistics"][role]={key:stats([float(x[key]) for x in selected])
                    for key in ["command_max_abs_error_m","command_rms_error_m","command_incident_max_abs_error_m",
                        "command_B_error_m","command_B_rate_m_s","command_q_P_rate_m_s3",
                        "command_on_margin_m_s2","crossing_max_abs_error_m","crossing_incident_max_abs_error_m"]}
            threshold_on=load_config(Path(report["source_directory"])/"normalized_config.yaml").controller.threshold_on_m_s2
            row["quiet_P_growth_m_s3"]=stats([x["prediction"]["q_P_rate_m_s3"] for x in report["quiet_cycles"]])
            row["quiet_on_margin_m_s2"]=stats([threshold_on-
                x["predicted_satellite_q_m_s2"] for x in report["quiet_cycles"] if x["predicted_satellite_q_m_s2"] is not None])
        row["first_response"]=json.loads((directory/"first_response.json").read_text(encoding="utf-8"))
        source=Path(report["source_directory"])
        with np.load(source/"results.npz") as f:
            arrays={key:f[key] for key in ["time_s","gap_error_m"]};final_rate=f["gap_rate_m_s"][-1]
        cfg=load_config(source/"normalized_config.yaml");events=load_switch_events(source)
        row["first_24h"]=first_day_metrics(arrays,cfg,events)
        time=float(arrays["time_s"][-1]);error=arrays["gap_error_m"][-1];attitude=restore_attitude(cfg,events,time)
        kp,kd=cfg.derived.controller_kp_s2_inv,cfg.derived.controller_kd_s_inv
        pterm,dterm=decompose_q(error[None,:],final_rate[None,:],kp,kd);candidates=[]
        for j in range(cfg.formation.satellite_count):
            if attitude.mode[j]!=0 or attitude.pending_target[j]!=-1:continue
            pred=predict_crossing(pterm[0,j],dterm[0,j],time,attitude.next_allowed_time[j],kp,kd,
                cfg.controller.threshold_on_m_s2,cfg.controller.control_period_s)
            if pred["predicted_command_s"] is not None:candidates.append(dict(satellite_id=j+1,
                q_P_m_s2=float(pterm[0,j]),q_D_m_s2=float(dterm[0,j]),**pred))
        applicable=bool(np.all(attitude.mode==0) and np.all(attitude.pending_target==-1))
        row["observation_end_state"]=dict(time_s=time,all_low_no_pending=applicable,
            velocity=velocity_distribution(error,final_rate),
            earliest_conditional_forecast=min(candidates,key=lambda x:(x["predicted_command_s"],x["satellite_id"])) if candidates else None,
            interpretation="fresh constant-D extrapolation, not validated beyond tmax; invalid no-actuation assumption if any non-LOW/pending")
        row["analysis_directory"]=str(directory.resolve());rows.append(row)
    factor_key=args.gain_axis+"_factor"
    rows.sort(key=lambda r:(r["days"],r[factor_key]))
    out=prepare_output_directory(args.output);write_json(out/"comparison.json",rows)
    with (out/"comparison.csv").open("w",encoding="utf-8-sig",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(rows[0]));writer.writeheader()
        writer.writerows({k:json.dumps(v,ensure_ascii=False) if isinstance(v,(dict,list)) else v for k,v in r.items()} for r in rows)
    lines=[f"|日数|{args.gain_axis.capitalize()}倍率|最大誤差 [m]|全体RMS [m]|終了時誤差 [m]|遷移数|完了cycle数|最後の切替 [日]|追加高度損失 [m]|",
        "|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for r in rows:
        lines.append(f"|{r['days']:.0f}|{r[factor_key]:.1f}|{r['e_max_m']:.3f}|{r['rms_m']:.3f}|{r['final_max_error_m']:.3f}|{r['transitions']}|{r['completed_cycles']}|{r['last_switch_day']:.5f}|{r['extra_altitude_loss_m']:.3f}|")
    (out/"table.md").write_text("\n".join(lines)+"\n",encoding="utf-8")
    print("\n".join(lines),flush=True)


if __name__=="__main__":main()
