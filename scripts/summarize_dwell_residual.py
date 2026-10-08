#!/usr/bin/env python3
"""Small comparison tables from residual reports; no reanalysis or simulation."""
import argparse
import csv
import json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from sat_string.io import prepare_output_directory,write_json


def stats(values):
    values=[v for v in values if v is not None]
    return dict(n=len(values),min=float(np.min(values)) if values else None,
                median=float(np.median(values)) if values else None,max=float(np.max(values)) if values else None)


def summarize(report):
    s=report["summary"];local=report["local_cycles"];quiet=report["quiet_cycles"]
    finite=[x for x in quiet if not x["right_censored"] and x["prediction"]["predicted_wait_s"] is not None]
    valid=[x for x in finite if x["prediction_valid_on_observed_window"]]
    first=next((x for x in local if x["satellite_id"]==12 and x["local_cycle_number"]==1),None)
    row=dict(dwell_s=s["dwell_time_s"],days=s["observation_end_s"]/86400,
        e_max_m=s["formation"]["max_abs_error_m"],rms_m=s["formation"]["rms_error_m"],
        final_max_error_m=s["formation"]["final_max_abs_error_m"],
        transitions=s["switching"]["n_transition_starts"],completed_cycles=s["local_completed_cycles"],
        cycles_per_satellite=s["switching"]["lhl_cycles_per_satellite"],
        incomplete_cycles=s["incomplete_cycles"],last_switch_day=s["switching"]["last_switch_time_s"]/86400,
        extra_altitude_loss_m=s["altitude"]["extra_loss_vs_all_low_m"],
        mean_initial_altitude_loss_m=s["altitude"]["mean_altitude_loss_m"],
        high_stable_total_s=s["switching"]["high_stable_total_s"],
        high_equivalent_total_s=s["switching"]["high_equivalent_total_s"],
        dwell_limited_returns=s["dwell_limited_returns"],threshold_limited_returns=s["threshold_limited_returns"],
        high_hold_times_s=s["high_hold_times_s"],
        local_V_incident_m_s=stats([x["V_res_seeds_m_s"] for x in local]),
        local_QD_own_abs_m_s2=stats([abs(x["own_q_D_m_s2"]) for x in local]),
        quiet_V_all_m_s=stats([x["V_res_all_m_s"] for x in quiet]),
        quiet_QD_all_abs_m_s2=stats([x["Q_D_res_all_m_s2"] for x in quiet]),
        local_next_start_interval_day=stats([x["next_same_satellite_start_interval_s"]/86400 if x["next_same_satellite_start_interval_s"] is not None else None for x in local]),
        quiet_actual_wait_day=stats([x["actual_next_wait_s"]/86400 if x["actual_next_wait_s"] is not None else None for x in quiet]),
        quiet_finite_uncensored_forecasts=len(finite),quiet_valid_forecasts=len(valid),
        quiet_prediction_target_matches=sum(x["prediction_target_matches_actual"] for x in valid),
        quiet_abs_relative_prediction_error=stats([abs(x["relative_prediction_error"]) for x in valid]),
        local_valid_forecasts=s["local_valid_finite_forecasts"],
        off_to_LOW_command_delay_s=stats([x["off_cross_to_low_command_s"] for x in local]),
        quiet_kinematic_identity_abs_error_m_s2=stats([abs(x["delta_q_P_m_s2"]-x["q_P_change_from_integrated_D_m_s2"]) for x in quiet]),
        settling=s["settling"],
        first_sat12={key:first[key] for key in ["high_command_s","high_complete_s","low_command_s","low_complete_s",
            "V_res_seeds_m_s","V_res_all_m_s","own_q_P_m_s2","own_q_D_m_s2","own_q_m_s2",
            "high_equivalent_time_s","q_at_HIGH_complete_m_s2","q_at_LOW_command_m_s2","nonlow_satellite_ids"]} if first else None)
    return row


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--analysis",type=Path,nargs="+",required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args();rows=[]
    for path in args.analysis:
        row=summarize(json.loads((path/"residual_metrics.json").read_text(encoding="utf-8")))
        row["analysis_directory"]=str(path.resolve());rows.append(row)
    rows.sort(key=lambda r:(r["days"],r["dwell_s"]))
    output=prepare_output_directory(args.output);write_json(output/"comparison.json",rows)
    with (output/"comparison.csv").open("w",encoding="utf-8-sig",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader()
        writer.writerows({k:json.dumps(v,ensure_ascii=False) if isinstance(v,(list,dict)) else v for k,v in row.items()} for row in rows)
    lines=["|観測日数|dwell [s]|最大誤差 [m]|全体RMS [m]|終了時誤差 [m]|遷移数|完了cycle数|最後の切替 [日]|追加高度損失 [m]|",
           "|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for r in rows:
        lines.append(f"|{r['days']:.0f}|{r['dwell_s']:.0f}|{r['e_max_m']:.3f}|{r['rms_m']:.3f}|{r['final_max_error_m']:.3f}|{r['transitions']}|{r['completed_cycles']}|{r['last_switch_day']:.5f}|{r['extra_altitude_loss_m']:.3f}|")
    (output/"table.md").write_text("\n".join(lines)+"\n",encoding="utf-8")
    print("\n".join(lines),flush=True)


if __name__=="__main__":main()
