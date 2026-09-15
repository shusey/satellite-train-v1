#!/usr/bin/env python3
"""Baseline first, then a seven-point existing-grid threshold/cost experiment."""
from __future__ import annotations
import argparse
import csv
import hashlib
import importlib.metadata
import json
from pathlib import Path
import platform
import subprocess
import sys
from time import perf_counter
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))

import numpy as np
from sat_string.config import config_from_dict, load_config, dump_normalized_config
from sat_string.control_cost import evaluate_cost, linear_modes
from sat_string.control_cost_plotting import plot_cost_report
from sat_string.io import save_case, load_case_arrays, load_switch_events, write_json
from sat_string.metrics import compute_case_metrics, save_metrics
from sat_string.simulator import STANDARD_CASES, simulate_case
from sat_string.switching_wave import threshold_scan_values


def write_csv(path, rows):
    if not rows:
        return
    with Path(path).open("w",encoding="utf-8",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def flat_row(label, config, m):
    p=m["propagation"]
    return dict(case=label,h_on_m_s2=config.controller.threshold_on_m_s2,
                h_off_m_s2=config.controller.threshold_off_m_s2,
                e_max_m=m["maximum_gap_error_m"],rms_m=m["rms_gap_error_m"],
                high_stable_total_h=m["activity"]["high_stable_total_s"]/3600,
                high_equivalent_total_h=m["activity"]["high_equivalent_total_s"]/3600,
                switch_count=m["activity"]["switch_count"],mean_a_loss_m=m["mean_a_loss_m"],
                extra_mean_a_loss_vs_low_m=m["extra_mean_a_loss_vs_low_m"],
                extra_max_a_loss_vs_low_m=m["extra_max_a_loss_vs_low_m"],
                energy_loss_J=m["total_energy_loss_J"],
                extra_energy_loss_J=m["extra_energy_loss_vs_low_J"],
                reached_count=p["reached_satellite_count"],range_km=p["range_m"]/1000,
                speed_left_m_s=p["left"]["speed_m_s"],speed_right_m_s=p["right"]["speed_m_s"],
                first_unreached_left=p["left"]["first_unreached_satellite_id"],
                first_unreached_right=p["right"]["first_unreached_satellite_id"],
                distal_equivalent_h=m["distal_high_equivalent_s"]/3600,
                distal_switch_count=m["distal_switch_count"],
                distal_extra_loss_sum_m=m["distal_extra_a_loss_sum_m"])


def case_summary(path, config, low):
    arrays=load_case_arrays(path)
    events=load_switch_events(path)
    m=evaluate_cost(arrays,config,events,low)
    write_json(path/"control_cost.json",m)
    rows=[]
    for i in range(config.formation.satellite_count):
        rows.append(dict(
            satellite_id=i+1,
            first_switch_start_s=m["propagation"]["first_switch_start_s_by_satellite"][i],
            switch_count=m["activity"]["switch_count_by_satellite"][i],
            high_stable_s=m["activity"]["high_stable_s_by_satellite"][i],
            high_equivalent_s=m["activity"]["high_equivalent_s_by_satellite"][i],
            delta_a_m=m["delta_a_m_by_satellite"][i],
            delta_altitude_m=m["delta_altitude_m_by_satellite"][i],
            delta_specific_energy_J_kg=m["delta_specific_energy_J_kg_by_satellite"][i],
            delta_energy_J=m["delta_energy_J_by_satellite"][i],
            extra_a_loss_vs_low_m=m["extra_a_loss_vs_low_m_by_satellite"][i],
            extra_energy_loss_vs_low_J=m["extra_energy_loss_vs_low_J_by_satellite"][i],
        ))
    write_csv(path/"satellite_metrics.csv",rows)
    write_csv(path/"gap_metrics.csv",[
        dict(gap_id=i+1,e_max_m=v,rms_m=m["rms_gap_error_m_by_gap"][i])
        for i,v in enumerate(m["maximum_gap_error_m_by_gap"])])
    return m


def run(config, case, directory):
    directory.mkdir(parents=True,exist_ok=False)
    start=perf_counter()
    messages=[]
    result=simulate_case(config,case,log_callback=messages.append)
    save_case(result,config,directory)
    save_metrics(compute_case_metrics(result.arrays(),config,result.switch_events),directory)
    (directory/"run.log").write_text("\n".join(messages)+f"\nElapsed: {perf_counter()-start:.6f} s\n",encoding="utf-8")
    print(f"Completed {directory.name} in {perf_counter()-start:.1f} s",flush=True)


def postprocess(output,config):
    low=load_case_arrays(output/"baseline"/"all_low_drag_baseline")
    baseline_rows=[]
    cases={}
    for case in STANDARD_CASES:
        directory=output/"baseline"/case.name
        m=case_summary(directory,config,low)
        baseline_rows.append(flat_row(case.name,config,m))
        cases[case.name]=load_case_arrays(directory)
    write_csv(output/"baseline_comparison.csv",baseline_rows)
    rows=[]
    manifest_path=output/"sweep_manifest.json"
    if manifest_path.exists():
        manifest=json.loads(manifest_path.read_text(encoding="utf-8-sig"))
        for item in manifest["cases"]:
            directory=output/item["directory"]
            c=load_config(directory/"normalized_config.yaml")
            m=case_summary(directory,c,low)
            rows.append(flat_row(directory.name,c,m))
        # Discrete sampled nondominance; no claim of global optimality.
        for r in rows:
            r["nondominated_e_max_vs_loss"]=not any(
                other["e_max_m"]<=r["e_max_m"] and
                other["extra_mean_a_loss_vs_low_m"]<=r["extra_mean_a_loss_vs_low_m"] and
                (other["e_max_m"]<r["e_max_m"] or other["extra_mean_a_loss_vs_low_m"]<r["extra_mean_a_loss_vs_low_m"])
                for other in rows)
        write_csv(output/"threshold_sweep.csv",rows)
    modes=linear_modes(config)
    write_json(output/"linear_modes.json",modes)
    write_csv(output/"linear_modes.csv",modes["modes"])
    paths=plot_cost_report(cases,rows,config,output/"plots")
    write_json(output/"summary.json",dict(baseline=baseline_rows,sweep=rows,figures=paths))
    print(json.dumps(baseline_rows,indent=2),flush=True)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--config",type=Path,default=Path("configs/baseline.yaml"))
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--stage",choices=["baseline","sweep","postprocess","all"],default="all")
    args=parser.parse_args()
    config=load_config(args.config)
    output=args.output
    if args.stage in ("baseline","all"):
        output.mkdir(parents=True,exist_ok=False)
        dump_normalized_config(config,output/"normalized_config.yaml")
        root=Path(__file__).resolve().parents[1]
        core_names=["simulator.py","controller.py","config.py","dynamics.py","attitude.py","atmosphere.py","disturbance.py","integrator.py"]
        write_json(output/"provenance.json",dict(
            source_commit=subprocess.check_output(["git","rev-parse","HEAD"],cwd=root,text=True).strip(),
            python=platform.python_version(),
            packages={p:importlib.metadata.version(p) for p in ["numpy","matplotlib","PyYAML","imageio-ffmpeg"]},
            core_sha256={name:hashlib.sha256((root/"src"/"sat_string"/name).read_bytes()).hexdigest() for name in core_names},
            case_source=str(args.config),
            sweep_rule="seven values selected from existing grid; off=on-(baseline on-off)",
        ))
        # No sweep runs until all three baseline simulations have succeeded.
        for case in STANDARD_CASES:
            run(config,case,output/"baseline"/case.name)
        postprocess(output,config)
    else:
        if not (output/"baseline"/"all_low_drag_baseline"/"results.npz").exists():
            raise FileNotFoundError("Run baseline stage first")
        saved=load_config(output/"normalized_config.yaml")
        if saved.to_dict()!=config.to_dict():
            raise ValueError("Requested config differs from saved baseline")
    if args.stage in ("sweep","all"):
        if (output/"sweep_manifest.json").exists():
            raise FileExistsError("Sweep already exists; use postprocess to regenerate summaries")
        values=np.array([1.0,1.3,1.38,1.5,1.6,1.8,2.0])*1e-7
        if not all(np.any(np.isclose(v,threshold_scan_values(),rtol=0,atol=1e-20)) for v in values):
            raise ValueError("Sweep values must belong to existing repository grid")
        width=config.controller.threshold_on_m_s2-config.controller.threshold_off_m_s2
        manifest={"cases":[]}
        for h_on in values:
            if np.isclose(h_on,config.controller.threshold_on_m_s2,rtol=0,atol=1e-20):
                directory=output/"baseline"/"controlled_disturbed"
            else:
                data=config.to_dict()
                data["controller"].update(threshold_on_m_s2=float(h_on),threshold_off_m_s2=float(h_on-width))
                c=config_from_dict(data)
                directory=output/"sweep"/f"h_on_{h_on:.3e}"
                run(c,STANDARD_CASES[0],directory)
            manifest["cases"].append(dict(directory=directory.relative_to(output).as_posix()))
        write_json(output/"sweep_manifest.json",manifest)
    if args.stage in ("sweep","postprocess","all"):
        postprocess(output,config)


if __name__=="__main__":
    main()
