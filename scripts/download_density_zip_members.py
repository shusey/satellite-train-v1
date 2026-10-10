"""Read a public ZIP directory and download only selected monthly members.

HTTP byte ranges avoid downloading the multi-year PANGAEA archives. ZIP CRC is
verified, and extracted files are written exclusively (never overwritten).
"""
import argparse
import hashlib
import io
import json
from pathlib import Path
import struct
import urllib.request
import zipfile
import zlib


def byte_range(url, start=None, end=None, suffix=None):
    spec=f'bytes=-{suffix}' if suffix else f'bytes={start}-{end}'
    request=urllib.request.Request(url,headers={'Range':spec,'User-Agent':'satellite-train-research/1.0'})
    with urllib.request.urlopen(request,timeout=60) as response:
        if response.status!=206:raise RuntimeError('Server must honor byte ranges; full archive refused')
        content_range=response.headers['Content-Range']
        data=response.read()
    return data,content_range


def directory(url):
    tail,content_range=byte_range(url,suffix=65536)
    at=tail.rfind(b'PK\x05\x06')
    if at<0:raise ValueError('Missing ZIP end record')
    end=list(struct.unpack('<4s4H2LH',tail[at:at+22]))
    size,offset=end[5:7]
    central,_=byte_range(url,offset,offset+size-1)
    end[6]=0
    archive=zipfile.ZipFile(io.BytesIO(central+struct.pack('<4s4H2LH',*end)))
    return archive.infolist(),content_range


def member(url, info):
    header,_=byte_range(url,info.header_offset,info.header_offset+29)
    fields=struct.unpack('<4s5H3L2H',header)
    if fields[0]!=b'PK\x03\x04':raise ValueError('Invalid local ZIP header')
    offset=info.header_offset+30+fields[-2]+fields[-1]
    compressed,_=byte_range(url,offset,offset+info.compress_size-1)
    data=zlib.decompress(compressed,-15) if info.compress_type==8 else compressed
    if len(data)!=info.file_size or zlib.crc32(data)!=info.CRC:raise ValueError('ZIP member integrity failure')
    return data


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--url',required=True);p.add_argument('--cache',type=Path,required=True)
    p.add_argument('--member',action='append',default=[]);p.add_argument('--receipt',type=Path,required=True)
    args=p.parse_args()
    if args.receipt.exists():raise FileExistsError(args.receipt)
    infos,content_range=directory(args.url)
    result=dict(url=args.url,archive_content_range=content_range,members=[dict(name=i.filename,bytes=i.file_size,compressed_bytes=i.compress_size) for i in infos],selected=[])
    for name in args.member:
        info=next(i for i in infos if i.filename==name)
        path=args.cache/Path(name).name
        if path.exists():data=path.read_bytes();status='cached'
        else:
            data=member(args.url,info);args.cache.mkdir(parents=True,exist_ok=True)
            with path.open('xb') as f:f.write(data)
            status='downloaded'
        if len(data)!=info.file_size or zlib.crc32(data)!=info.CRC:raise ValueError('Cached ZIP member differs')
        result['selected'].append(dict(name=name,path=str(path),bytes=len(data),status=status,sha256=hashlib.sha256(data).hexdigest()))
    args.receipt.parent.mkdir(parents=True,exist_ok=True)
    args.receipt.write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(result['selected'] if args.member else result['members'],indent=2))


if __name__=='__main__':main()
