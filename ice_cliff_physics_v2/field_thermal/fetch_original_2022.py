"""Reproduce bounded original-data retrieval from recorded immutable URLs/hashes."""
import argparse,hashlib,json,urllib.request
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def main():
 p=argparse.ArgumentParser();p.add_argument('--out',type=Path,required=True);a=p.parse_args();m=json.loads((ROOT/'original_2022/SOURCE_MANIFEST.json').read_text());a.out.mkdir(parents=True,exist_ok=True)
 for item in m['members']:
  path=a.out/item['path']
  if path.exists():
   if hashlib.sha256(path.read_bytes()).hexdigest()!=item['sha256']:raise ValueError('existing output differs; refused overwrite')
   continue
  with urllib.request.urlopen(item['url'],timeout=30) as res:data=res.read(2000001)
  if len(data)!=item['bytes'] or hashlib.sha256(data).hexdigest()!=item['sha256']:raise ValueError('source version/hash mismatch')
  path.write_bytes(data)
 print('All pinned original CSV/config files verified')
if __name__=='__main__':main()
