"""Recover selected public author-mirror CSVs; never execute archive content."""
import argparse
import hashlib
import json
from pathlib import Path, PurePosixPath
import struct
from urllib.request import Request, urlopen
import zlib

URL = 'https://people.ee.ethz.ch/~werderm/pm-data/PlaineMorteLake_hydro_thermonynamics.zip'
SIZE = 337848233
ETAG = '"142327a9-5cc08324a4e70"'
PREFIX = 'PlaineMorteLake_hydro_thermonynamics'

def get_range(start, end):
    req = Request(URL, headers={'Range': f'bytes={start}-{end}'})
    with urlopen(req, timeout=60) as response:
        if response.status != 206:
            raise ValueError('Server did not provide a bounded range; stop rather than download all')
        if response.headers.get('Content-Range') != f'bytes {start}-{end}/{SIZE}':
            raise ValueError('Unexpected source archive range/version')
        if response.headers.get('ETag') != ETAG:
            raise ValueError('Archive version changed; audit before accepting it')
        data = response.read(end-start+2)
    if len(data) != end-start+1:
        raise ValueError('Response size mismatch')
    return data

def directory(tail, start):
    pos = tail.rfind(b'PK\x05\x06')
    if pos < 0:
        raise ValueError('Missing end directory')
    end = struct.unpack_from('<4s4H2LH', tail, pos)
    if end[1] or end[2] or end[3] != end[4]:
        raise ValueError('Unsupported split archive')
    off = end[6]-start
    entries = []
    while tail[off:off+4] == b'PK\x01\x02':
        v = struct.unpack_from('<4s6H3L5H2L', tail, off)
        n, extra, comment = v[10:13]
        name = tail[off+46:off+46+n].decode('utf8')
        entries.append(dict(name=name, method=v[4], crc32=v[7], compressed=v[8], uncompressed=v[9], offset=v[-1]))
        off += 46+n+extra+comment
    if len(entries) != end[4]:
        raise ValueError('Directory entry count mismatch')
    return entries

def extract(entries, blocks, output):
    records = []
    for e in entries:
        name = e['name']
        selected = (('/data/Channel_elevation/' in name or '/products/' in name) and name.endswith('.csv')) or name == PREFIX+'/README.txt'
        if not selected:
            continue
        path = PurePosixPath(name)
        if path.is_absolute() or '..' in path.parts or path.parts[0] != PREFIX:
            raise ValueError('Unsafe member path')
        for start, block in blocks:
            offset = e['offset']-start
            if offset < 0 or offset+30 > len(block):
                continue
            v = struct.unpack_from('<4s5H3L2H', block, offset)
            if v[0] != b'PK\x03\x04' or v[2] & 1:
                raise ValueError('Invalid or encrypted local member')
            a = offset+30+v[-2]+v[-1]
            local_name = block[offset+30:offset+30+v[-2]].decode('utf8')
            if local_name != name:
                raise ValueError('Local/directory names disagree')
            if a+e['compressed'] > len(block):
                continue
            packed = block[a:a+e['compressed']]
            if e['method'] not in (0, 8):
                raise ValueError('Unsupported compression')
            data = zlib.decompress(packed, -15) if e['method'] == 8 else packed
            if len(data) != e['uncompressed'] or zlib.crc32(data) & 0xffffffff != e['crc32']:
                raise ValueError('Member integrity mismatch')
            rel = Path(*path.parts[1:])
            dest = output/rel
            dest.parent.mkdir(parents=True, exist_ok=True)
            dest.write_bytes(data)
            records.append(dict(e, path=rel.as_posix(), sha256=hashlib.sha256(data).hexdigest(), verified_crc32=True))
            break
        else:
            raise ValueError('Selected member outside acquired ranges: '+name)
    return records

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    if args.output.exists() and any(args.output.iterdir()):
        raise ValueError('Use an empty output directory to preserve earlier source records')
    tail_start = SIZE-131072
    tail = get_range(tail_start, SIZE-1)
    entries = directory(tail, tail_start)
    blocks = [(0, get_range(0, 8191)), (331460392, get_range(331460392, 337773946))]
    args.output.mkdir(parents=True, exist_ok=True)
    records = extract(entries, blocks, args.output)
    manifest = dict(mirror=URL, archive_total_bytes=SIZE, archive_etag=ETAG, full_archive_downloaded=False,
                    central_directory_entries=len(entries), members=records)
    (args.output/'SOURCE_MANIFEST.json').write_text(json.dumps(manifest, indent=2)+'\n')
    print(f'Verified {len(records)} selected members; no archive code executed')

if __name__ == '__main__':
    main()
