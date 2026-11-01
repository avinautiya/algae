"""Bounded metadata/range audit; never downloads a full large tar implicitly."""
import argparse
import hashlib
import json
from pathlib import Path
import urllib.request
import re

RECORD='https://zenodo.org/api/records/13073856'


def request(url,limit,headers=None):
    req=urllib.request.Request(url,headers={'User-Agent':'glacier-research-data-audit/1.0',**(headers or {})})
    with urllib.request.urlopen(req,timeout=25) as response:
        if headers and 'Range' in headers and response.status!=206:
            raise RuntimeError('range response not206; refusing full archive')
        if headers and 'Range' in headers:
            match=re.fullmatch(r'bytes=(\d+)-(\d+)',headers['Range'])
            received=re.fullmatch(r'bytes (\d+)-(\d+)/(\d+|\*)',response.headers.get('Content-Range',''))
            if not match or not received or match.groups()!=received.groups()[:2]:
                raise RuntimeError('Content-Range does not match requested byte interval')
        data=response.read(limit+1)
        if len(data)>limit:raise RuntimeError('response exceeds declared byte cap')
        return data,dict(response.headers)


def audit(out,max_headers=100):
    if isinstance(max_headers,bool) or not isinstance(max_headers,int) or not 1<=max_headers<=1000:
        raise ValueError('max_headers must be an integer from 1 to 1000')
    out.mkdir(parents=True,exist_ok=True)
    payload,_=request(RECORD,300000)
    meta=json.loads(payload);(out/'zenodo_13073856_metadata.json').write_bytes(payload)
    files=meta['files'];target=next(f for f in files if f['key']=='supplementary-data-to-article.tar')
    url=target['links']['self'];offset=0;entries=[];status='header_budget_reached'
    for _ in range(max_headers):
        header,headers=request(url,512,{'Range':f'bytes={offset}-{offset+511}'})
        if len(header)!=512:raise RuntimeError('incomplete tar header')
        if not any(header):status='end_of_archive';break
        name=header[:100].split(b'\0',1)[0].decode('utf-8',errors='replace')
        prefix=header[345:500].split(b'\0',1)[0].decode('utf-8',errors='replace')
        if prefix:name=prefix+'/'+name
        size=int(header[124:136].strip(b'\0 ') or b'0',8)
        # Tar header checksum validation, independent of the remote transport.
        declared=int(header[148:156].strip(b'\0 ') or b'0',8)
        computed=sum(header[:148])+8*32+sum(header[156:])
        if declared!=computed:raise RuntimeError('tar header checksum mismatch')
        entries.append(dict(name=name,size=size,header_offset=offset,data_offset=offset+512,type=chr(header[156]) if header[156] else '0'))
        offset+=512+((size+511)//512)*512
    record=dict(source_record=RECORD,source_doi='10.5281/zenodo.13073856',
                metadata_sha256=hashlib.sha256(payload).hexdigest(),archive=target,
                status=status,entries=entries,
                warning='range header inventory is not a full archive checksum or extracted-data validation')
    (out/'profile_archive_inventory.json').write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(dict(status=status,headers=len(entries),entries=entries),indent=2),flush=True)


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--output',type=Path,required=True);p.add_argument('--max-headers',type=int,default=100)
    a=p.parse_args();audit(a.output,a.max_headers)
