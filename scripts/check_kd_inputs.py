#!/usr/bin/env python3
"""Check each saved Kd input against the frozen baseline; no simulation."""
import argparse
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from sat_string.config import load_config
from sat_string.io import load_metadata,prepare_output_directory,write_json


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--runs",type=Path,required=True);p.add_argument("--output",type=Path,required=True)
    p.add_argument("--gain-axis",choices=["kd","kp"],default="kd")
    p.add_argument("--factors",type=float,nargs="+",help="Expected new factors; defaults retain existing experiments")
    args=p.parse_args();baseline=load_config("tests/fixtures/regression_24h/config.yaml")
    factors=args.factors if args.factors is not None else ([.3,2.,2.5,3.] if args.gain_axis=="kp" else [.3,3.])
    checks={}
    for path in sorted(args.runs.glob("*/normalized_config.yaml")):
        cfg=load_config(path);field="kp_s2_inv" if args.gain_axis=="kp" else "kd_s_inv"
        value=getattr(cfg.derived,"controller_"+field);factor=value/getattr(baseline.derived,"controller_"+field)
        expected=baseline.to_dict();expected["controller"][field]=value
        expected["derived"]["controller_"+field]=value
        expected["simulation"]["duration_s"]=cfg.simulation.duration_s
        case=load_metadata(path.parent)["case"]
        checks[path.parent.name]=dict([(args.gain_axis+"_factor",factor),("days",cfg.simulation.duration_s/86400)],
            only_gain_and_duration_changed=cfg.to_dict()==expected,
            control_and_disturbance_enabled=case["control_enabled"] and case["disturbance_enabled"],
            requested_factor=factor in factors)
        checks[path.parent.name]["only_"+args.gain_axis+"_and_duration_changed"]=cfg.to_dict()==expected
    passed=bool(checks) and all(all(r[k] for k in ["only_gain_and_duration_changed","control_and_disturbance_enabled","requested_factor"]) for r in checks.values())
    out=prepare_output_directory(args.output);write_json(out/"inputs.json",dict(passed=passed,runs=checks))
    print(dict(passed=passed,runs=len(checks)),flush=True)
    return 0 if passed else 1


if __name__=="__main__":raise SystemExit(main())
