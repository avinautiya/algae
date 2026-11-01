"""Fetch selected published tar members with strict byte ranges and provenance.

The source archive MD5 is recorded, not claimed verified for a partial download.
No GIS or data file from the archive is executed.
"""
import argparse
import hashlib
import json
from pathlib import Path
from audit_zenodo_profiles import request

MEMBERS = ('README.txt', 'cliff-profiles/cliff_profile_point_clouds.csv',
           'cliff-profiles/cliff_profile_lines.gpkg', 'accuracy/survey-error-statistics.csv')


def fetch(root):
    inventory = json.loads((root/'profile_archive_inventory.json').read_text())
    url=inventory['archive']['links']['self']
    records=[]
    for name in MEMBERS:
        entry=next(e for e in inventory['entries'] if e['name']==name)
        if entry['type']!='0' or not 0<entry['size']<=30_000_000:
            raise ValueError('unexpected member type or byte budget')
        lo=entry['data_offset'];hi=lo+entry['size']-1
        payload,headers=request(url,entry['size'],{'Range':f'bytes={lo}-{hi}'})
        if len(payload)!=entry['size']:raise RuntimeError('truncated member')
        target=root/'published_profiles'/Path(name).name
        target.parent.mkdir(parents=True,exist_ok=True)
        target.write_bytes(payload)
        records.append(dict(member=name,bytes=len(payload),local_file=str(target.name),
                            sha256=hashlib.sha256(payload).hexdigest(),
                            content_range=headers.get('Content-Range') or headers.get('content-range')))
        print(f"Fetched {name}: {len(payload)} bytes",flush=True)
    (root/'profile_download_manifest.json').write_text(json.dumps(dict(
        source_doi=inventory['source_doi'],source_archive=inventory['archive'],
        archive_checksum_verified=False,
        verification='HTTP byte-range, tar-header checksum and member length; locally computed SHA256. Full archive MD5 NOT verified.',
        members=records),indent=2)+'\n')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--root',type=Path,required=True)
    fetch(p.parse_args().root)
