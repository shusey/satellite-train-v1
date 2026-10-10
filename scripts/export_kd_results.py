#!/usr/bin/env python3
"""Curate small tables/configs/events; full endpoint arrays remain in local runs."""
import argparse
import json
from pathlib import Path
import shutil
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from sat_string.io import prepare_output_directory,write_json
from summarize_dwell_residual import stats
from analyze_dwell_residual import csv_rows


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--comparison",type=Path,required=True);p.add_argument("--output",type=Path,required=True)
    p.add_argument("--gain-axis",choices=["kd","kp"],default="kd")
    p.add_argument("--factors",type=float,nargs="+",help="Export only new factors; keep the complete comparison table")
    args=p.parse_args();rows=json.loads((args.comparison/"comparison.json").read_text(encoding="utf-8"))
    out=prepare_output_directory(args.output)
    for name in ["comparison.json","comparison.csv","table.md"]:shutil.copy2(args.comparison/name,out/name)
    for row in rows:
        if args.factors is not None and row[args.gain_axis+'_factor'] not in args.factors:continue
        analysis=Path(row["analysis_directory"]);r=json.loads((analysis/"residual_metrics.json").read_text(encoding="utf-8"))
        provenance=json.loads((analysis/"provenance.json").read_text(encoding="utf-8"));run=Path(provenance["source"])
        case=prepare_output_directory(out/f"{args.gain_axis}_{row[args.gain_axis+'_factor']:g}_{row['days']:.0f}d")
        for name in ["summary.json","first_response.json","provenance.json"]:shutil.copy2(analysis/name,case/name)
        shutil.copy2(run/"normalized_config.yaml",case/"normalized_config.yaml")
        for name in [args.gain_axis+"_initialization.json","continuation.json"]:
            if (run/name).exists():shutil.copy2(run/name,case/name)
        # Sixty-day events include the complete thirty-day prefix; avoid duplicating it.
        if row["days"]==60:shutil.copy2(run/"switch_events.csv",case/"switch_events.csv")
        for name in ["group_residuals.csv","quiet_residuals.csv","first_local_cycles.csv"]+(["threshold_events.csv","off_before_dwell.csv"] if args.gain_axis=="kp" else []):
            if (analysis/name).exists():shutil.copy2(analysis/name,case/name)
        local=r["local_cycles"];summary=[];samples=[]
        for satellite in sorted({x["satellite_id"] for x in local}):
            own=[x for x in local if x["satellite_id"]==satellite]
            data=dict(satellite_id=satellite,completed_cycles=len(own))
            for metric in ["V_res_all_m_s","V_res_rms_all_m_s","V_res_median_abs_all_m_s",
                    "V_res_seeds_m_s","own_q_D_m_s2","worsening_gap_fraction","worsening_weight_fraction",
                    "actual_next_wait_s","next_same_satellite_start_interval_s"]:
                data[metric]=stats([x[metric] for x in own])
            data.update(HIGH_stable_total_s=sum(x["high_stable_time_s"] for x in own),
                HIGH_equivalent_total_s=sum(x["high_equivalent_time_s"] for x in own),
                threshold_limited_returns=sum(x["return_driver"]!="dwell_limited" for x in own))
            summary.append(data)
            samples.extend(x for i,x in enumerate(own) if i<3 or i>=len(own)-3)
        csv_rows(case/"per_satellite_residual_statistics.csv",summary)
        csv_rows(case/"local_cycle_samples.csv",samples)
        gates=[]
        import csv
        with (analysis/"return_gates.csv").open(encoding="utf-8-sig") as f:gates=list(csv.DictReader(f))
        csv_rows(case/"blocked_return_gates.csv",[x for x in gates if x["off_not_met_at_gate"]=="True"])
        write_json(case/"local_files.json",dict(full_cycle_csv=str(analysis/"local_cycles.csv"),
            full_threshold_csv=str(analysis/"threshold_events.csv"),full_return_gate_csv=str(analysis/"return_gates.csv"),
            signed_endpoint_npz=str(analysis/"cycle_end_states.npz"),physics_npz=str(run/"results.npz"),
            sample_policy="first three and last three completed local cycles per satellite; not the full distribution"))
    count=sum(args.factors is None or r[args.gain_axis+'_factor'] in args.factors for r in rows)
    print(f"Curated {count} cases: {out}",flush=True)


if __name__=="__main__":main()
