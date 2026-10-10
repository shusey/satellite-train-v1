"""Download only explicitly selected public files, preserving an immutable cache."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import hashlib
import json
from pathlib import Path
import urllib.request


def download(entry, root, timeout=45):
    path=root/entry['filename']
    if path.exists():
        data=path.read_bytes()
        if entry.get('sha256') and hashlib.sha256(data).hexdigest()!=entry['sha256']:
            raise ValueError('Cached file checksum differs from manifest')
        return dict(**entry, status='cached', bytes=len(data), sha256=hashlib.sha256(data).hexdigest())
    try:
        req=urllib.request.Request(entry['url'],headers={'User-Agent':'satellite-train-research/1.0'})
        with urllib.request.urlopen(req, timeout=timeout) as response:
            data=response.read()
            content_type=response.headers.get('Content-Type')
        if entry.get('sha256') and hashlib.sha256(data).hexdigest()!=entry['sha256']:
            raise ValueError('Downloaded file checksum differs from manifest')
        path.parent.mkdir(parents=True,exist_ok=True)
        with path.open('xb') as stream: stream.write(data)
        return dict(**entry,status='downloaded',bytes=len(data),sha256=hashlib.sha256(data).hexdigest(),content_type=content_type)
    except Exception as exc:
        return dict(**entry,status='failed',error=str(exc))


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--manifest',type=Path,required=True)
    p.add_argument('--cache',type=Path,required=True)
    p.add_argument('--receipt',type=Path,required=True)
    p.add_argument('--timeout',type=int,default=45)
    args=p.parse_args()
    if args.receipt.exists():raise FileExistsError(args.receipt)
    entries=json.loads(args.manifest.read_text(encoding='utf-8-sig'))
    for entry in entries:
        target=(args.cache/entry['filename']).resolve()
        if not target.is_relative_to(args.cache.resolve()):raise ValueError('Path outside cache')
    with ThreadPoolExecutor(max_workers=4) as pool:
        result=list(pool.map(lambda e:download(e,args.cache,args.timeout),entries))
    args.receipt.parent.mkdir(parents=True,exist_ok=True)
    args.receipt.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding='utf-8')
    print(json.dumps(result,indent=2,ensure_ascii=False))


if __name__=='__main__':main()
