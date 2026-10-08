#!/usr/bin/env python3
"""Check old-file hashes and exact saved prefixes, without simulation."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from sat_string.io import prepare_output_directory,write_json,load_switch_events


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--manifest",type=Path,required=True)
    p.add_argument("--runs",type=Path,required=True)
    p.add_argument("--output",type=Path,required=True)
    args=p.parse_args()
    manifest=json.loads(args.manifest.read_text(encoding="utf-8"))
    changed=[]
    for name,expected in manifest.items():
        path=Path(name)
        if not path.is_file():changed.append(name);continue
        with path.open("rb") as stream:
            if hashlib.file_digest(stream,"sha256").hexdigest()!=expected:changed.append(name)
    prefixes=[]
    for provenance in sorted(args.runs.glob("*/continuation.json")):
        data=json.loads(provenance.read_text(encoding="utf-8"))
        source=Path(data["source_directory"]);destination=provenance.parent
        checks={}
        with np.load(source/"results.npz") as before,np.load(destination/"results.npz") as after:
            length=len(before["time_s"])
            for key in before.files:
                checks[key]=bool(np.array_equal(before[key],after[key][:length]))
        events=load_switch_events(source);new=load_switch_events(destination)
        event_prefix=events==new[:len(events)]
        prefixes.append(dict(source_directory=str(source),destination_directory=str(destination),
                             arrays=checks,events_match=event_prefix,
                             passed=all(checks.values()) and event_prefix))
    report=dict(old_files_checked=len(manifest),changed_or_missing_files=changed,prefixes=prefixes,
                passed=not changed and bool(prefixes) and all(p["passed"] for p in prefixes))
    output=prepare_output_directory(args.output);write_json(output/"preservation.json",report)
    print(dict(passed=report["passed"],old_files=len(manifest),continuations=len(prefixes)))
    return 0 if report["passed"] else 1


if __name__=="__main__":raise SystemExit(main())
