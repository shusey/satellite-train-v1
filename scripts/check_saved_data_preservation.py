#!/usr/bin/env python3
"""Verify a pre-work SHA256 manifest for an offline study with no continuations."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from sat_string.io import prepare_output_directory, write_json


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();manifest=json.loads(args.manifest.read_text(encoding='utf-8'))
    if not manifest:raise ValueError('A nonempty pre-work manifest is required')
    changed=[]
    for name,expected in manifest.items():
        path=Path(name)
        if not path.is_file():changed.append(name);continue
        with path.open('rb') as f:
            if hashlib.file_digest(f,'sha256').hexdigest()!=expected:changed.append(name)
    out=prepare_output_directory(args.output)
    report=dict(passed=not changed,old_files_checked=len(manifest),changed_or_missing_files=changed)
    write_json(out/'preservation.json',report);print(report,flush=True)
    return 0 if report['passed'] else 1


if __name__=='__main__':raise SystemExit(main())
