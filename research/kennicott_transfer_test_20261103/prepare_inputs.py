"""Restore verified feature CSVs from compact, hash-bound transport blobs."""
import gzip,hashlib,json
from pathlib import Path
ROOT=Path(__file__).resolve().parent
if __name__=='__main__':
 manifest=json.loads((ROOT/'FEATURE_TRANSPORT.json').read_text())
 for entry in manifest['files']:
  target=ROOT/entry['path'];compressed=ROOT/entry['transport_path']
  if target.exists():
   if hashlib.sha256(target.read_bytes()).hexdigest()!=entry['sha256']:raise ValueError('Existing feature differs: '+entry['path'])
   continue
  raw=compressed.read_bytes()
  if hashlib.sha256(raw).hexdigest()!=entry['transport_sha256']:raise ValueError('Transport checksum mismatch')
  data=gzip.decompress(raw)
  if hashlib.sha256(data).hexdigest()!=entry['sha256']:raise ValueError('Feature checksum mismatch')
  target.parent.mkdir(parents=True,exist_ok=True)
  with target.open('xb')as f:f.write(data)
 print('All feature inputs verified')
