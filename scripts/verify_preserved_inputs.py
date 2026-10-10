#!/usr/bin/env python3
"""Read-only SHA256 check of files recorded before research work."""
import argparse
import hashlib
import json
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'src'))
from sat_string.io import write_json


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest',type=Path,required=True);p.add_argument('--output',type=Path,required=True)
    args=p.parse_args()
    if args.output.exists():raise FileExistsError(args.output)
    manifest=json.loads(args.manifest.read_text(encoding='utf-8'));root=Path.cwd().resolve();changes=[]
    for relative,expected in manifest.items():
        path=(root/relative).resolve()
        if not path.is_relative_to(root):raise ValueError('Manifest path outside repository')
        if not path.is_file():changes.append(dict(path=relative,status='missing'));continue
        with path.open('rb') as f:actual=hashlib.file_digest(f,'sha256').hexdigest()
        if actual!=expected:changes.append(dict(path=relative,status='changed',expected=expected,actual=actual))
    result=dict(passed=not changes,files_checked=len(manifest),changes=changes)
    write_json(args.output,result);print(result,flush=True)
    if changes:raise RuntimeError('Existing research files changed')


if __name__=='__main__':main()
