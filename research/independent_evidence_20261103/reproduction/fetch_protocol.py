"""Fetch exact public protocol3362374 and Table3; existing directories refuse."""
import argparse,hashlib,json,urllib.request
from pathlib import Path
NAMES=['WGDCG_MeltModel_Intercomparison_Experiment_Description_FINAL_RUNS.pdf','README_glacierdata_updated.txt','MMIC_model_run_dates.xlsx']
def get(url,maxbytes):
 with urllib.request.urlopen(url,timeout=60) as r:b=r.read(maxbytes+1)
 if len(b)>maxbytes:raise ValueError('public resource exceeds declared fetch cap')
 return b
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
 if a.out.exists():raise ValueError('new output directory required')
 raw=get('https://zenodo.org/api/records/3362374',2_000_000);d=json.loads(raw)
 files={f['key']:f for f in d['files']};prepared={}
 for name in NAMES:
  f=files[name]
  if f['size']>1_000_000:raise ValueError('unexpected protocol size')
  b=get(f['links']['self'],1_000_000);algorithm,expected=f['checksum'].split(':')
  if algorithm!='md5' or len(b)!=f['size'] or hashlib.md5(b).hexdigest()!=expected:raise ValueError('protocol checksum mismatch')
  prepared[name]=b
 prepared['published_table3.xlsx']=get('https://tc.copernicus.org/articles/20/1895/2026/tc-20-1895-2026-t03.xlsx',1_000_000)
 prepared['PROTOCOL_RECORD.json']=raw
 a.out.mkdir(parents=True)
 for name,b in prepared.items():(a.out/name).write_bytes(b)
 print(json.dumps(dict(files=list(prepared),license=d['metadata'].get('license'),sha256={k:hashlib.sha256(b).hexdigest() for k,b in prepared.items()}),indent=2))
if __name__=='__main__':main()
