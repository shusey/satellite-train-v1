#!/usr/bin/env python3
"""Verify old file hashes, frozen inputs and exact reused simulation prefixes."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from sat_string.config import load_config
from sat_string.io import load_switch_events,prepare_output_directory,write_json


def digest(path):
    h=hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda:stream.read(1024*1024),b""):h.update(block)
    return h.hexdigest()


def prefix_check(short,long,cut=None):
    with np.load(short/"results.npz") as a,np.load(long/"results.npz") as b:
        t=a["time_s"];end=float(t[-1] if cut is None else cut);count=int(np.searchsorted(t,end+1e-9))
        checks={key:bool(np.array_equal(a[key][:count],b[key][:count])) for key in a.files}
    original=[e for e in load_switch_events(short) if e["command_time"]<=end]
    extended=[e for e in load_switch_events(long) if e["command_time"]<=end]
    return dict(source=str(short),destination=str(long),through_s=end,arrays=checks,
                events_unchanged=original==extended,passed=all(checks.values()) and original==extended)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root",type=Path,default=Path("outputs/dwell_residual_20261008"))
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    manifest_path=args.root/"preserved_inputs_sha256.json"
    manifest=json.loads(manifest_path.read_text(encoding="utf-8"))
    changed=[path for path,value in manifest.items() if not Path(path).is_file() or digest(Path(path))!=value]
    frozen=load_config("tests/fixtures/regression_24h/config.yaml")
    config_checks={};prefixes=[]
    reference=Path("outputs/long_horizon_q_20261008/regression_24h/suite/controlled_disturbed")
    for dwell in [150,600]:
        for days in [30,60]:
            run=args.root/"runs"/f"dwell_{dwell}_{days}d"
            expected=frozen.to_dict();expected["controller"]["dwell_time_s"]=float(dwell);expected["simulation"]["duration_s"]=float(days)*86400
            config_checks[run.name]=load_config(run/"normalized_config.yaml").to_dict()==expected
        short=args.root/"runs"/f"dwell_{dwell}_30d";long=args.root/"runs"/f"dwell_{dwell}_60d"
        prefixes.append(prefix_check(reference,short,20420.))
        prefixes.append(prefix_check(short,long))
    report=dict(preserved_old_files=len(manifest),changed_or_missing_old_files=changed,
        input_manifest_sha256=digest(manifest_path),configuration_checks=config_checks,prefix_checks=prefixes,
        passed=not changed and all(config_checks.values()) and all(x["passed"] for x in prefixes))
    output=prepare_output_directory(args.output);write_json(output/"preservation.json",report)
    print(dict(passed=report["passed"],old_files=len(manifest),prefixes=len(prefixes)),flush=True)
    return 0 if report["passed"] else 1


if __name__=="__main__":raise SystemExit(main())
