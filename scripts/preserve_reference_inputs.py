"""Record/check old tracked hashes and old local-run file sizes/mtimes."""
import argparse,hashlib,json,subprocess
from run_reference_dispersion import ROOT,DATA
from sat_string.io import write_json


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    p=argparse.ArgumentParser();p.add_argument('--check',action='store_true');args=p.parse_args()
    manifest=ROOT/'outputs/reference_dispersion_20261012/preserved_inputs.json'
    if not args.check:
        if manifest.exists():raise FileExistsError(manifest)
        names=subprocess.check_output(['git','ls-tree','-r','--name-only','HEAD'],cwd=ROOT,text=True).splitlines()
        hashes={name:sha(ROOT/name) for name in names}
        stats={str(f.relative_to(ROOT)):[f.stat().st_size,f.stat().st_mtime_ns]
            for directory in ('outputs','runs') if (ROOT/directory).exists()
            for f in (ROOT/directory).rglob('*') if f.is_file() and 'reference_dispersion_20261012' not in f.parts}
        write_json(manifest,dict(base_commit=subprocess.check_output(['git','rev-parse','HEAD'],cwd=ROOT,text=True).strip(),
            tracked_sha256=hashes,old_run_size_mtime_ns=stats))
        print('Recorded',len(hashes),'tracked hashes and',len(stats),'old local files')
    else:
        old=json.loads(manifest.read_text());changed=[]
        for name,digest in old['tracked_sha256'].items():
            if not (ROOT/name).exists() or sha(ROOT/name)!=digest:changed.append(name)
        for name,stat in old['old_run_size_mtime_ns'].items():
            f=ROOT/name
            if not f.exists() or [f.stat().st_size,f.stat().st_mtime_ns]!=stat:changed.append(name)
        report=dict(passed=not changed,changed=changed,base_commit=old['base_commit'],
            tracked_hashes_checked=len(old['tracked_sha256']),old_run_files_stat_checked=len(old['old_run_size_mtime_ns']),
            limitation='Tracked contents checked by SHA256; large old local runs checked by size+mtime, not duplicate full hashing')
        write_json(DATA/'preservation.json',report);print(report)
        if changed:raise SystemExit(1)


if __name__=='__main__':main()
