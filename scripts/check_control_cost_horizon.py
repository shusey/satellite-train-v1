#!/usr/bin/env python3
"""Check whether 24-hour cost savings survive the repository's 72-hour horizon."""
import argparse
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from sat_string.config import load_config,config_from_dict
from sat_string.io import load_case_arrays,write_json
from sat_string.simulator import STANDARD_CASES
from run_control_cost_analysis import run,case_summary,flat_row,write_csv


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--output",type=Path,required=True)
    p.add_argument("--horizon-config",type=Path,default=Path("runs/standard_suite/normalized_config.yaml"))
    args=p.parse_args()
    source=load_config(args.output/"normalized_config.yaml")
    end=load_config(args.horizon_config).simulation.duration_s
    directory=args.output/"horizon_check"
    directory.mkdir(parents=True,exist_ok=False)
    spec=[("A_low",source,STANDARD_CASES[2]),
          ("B_current",source,STANDARD_CASES[0])]
    for name in ["h_on_1.800e-07","h_on_2.000e-07"]:
        spec.append((name,load_config(args.output/"sweep"/name/"normalized_config.yaml"),STANDARD_CASES[0]))
    rows=[]
    for label,cfg,case in spec:
        data=cfg.to_dict()
        data["simulation"]["duration_s"]=end
        c=config_from_dict(data)
        run(c,case,directory/label)
        low=load_case_arrays(directory/"A_low")
        m=case_summary(directory/label,c,low)
        rows.append(flat_row(label,c,m))
    write_csv(directory/"comparison.csv",rows)
    write_json(directory/"summary.json",dict(
        horizon_duration_s=end,horizon_source=str(args.horizon_config),
        only_changed_physical_experiment_field="simulation.duration_s",
        rows=rows))
    from sat_string.control_cost_plotting import plot_horizon_cost
    plot_horizon_cost({label:load_case_arrays(directory/label) for label,_,_ in spec}, source, args.output/"plots")
    for row in rows:
        print(row,flush=True)


if __name__=="__main__":
    main()

