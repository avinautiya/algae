"""Read verified Changri Nup component data without inventing melt labels.

No training, target-dependent snow features or model comparison is performed.
"""
import argparse,csv,datetime as dt,hashlib,io,json,math,re,zipfile
from pathlib import Path
MD5='4797deeda5c93a1895f08efb42183086'
ROOT='CN-5470_20152017/'
FORBIDDEN={'SNOW','T_z','RH_z','FF_z','DIST','HEIGHT','ABL','ACC'}
ALLOWED=['T','RH','FF','DIR','P','SWIN','SWOUT','LWIN','LWOUT','PP']
NEPAL=dt.timezone(dt.timedelta(hours=5,minutes=45))
def scalar(s):
 s=s.strip()
 if s.casefold() in {'nan','','na','n/a'}:return None
 if re.fullmatch(r'[+-]?\d+,\d+[eE][+-]?\d+',s):s=s.replace(',','.')
 v=float(s)
 if not math.isfinite(v):raise ValueError('unsupported nonfinite source value')
 return v
def time_utc(row):
 # Actual source declares Nepal time. No naive-machine-time conversion.
 hour=int(row.get('HOUR','0'));minute=int(row.get('MINUTE','0'))
 if hour==24 and minute!=0:raise ValueError('24-hour timestamp requires minute0')
 t=dt.datetime(*(int(row[k]) for k in ['YEAR','MONTH','DAY']),0 if hour==24 else hour,minute,tzinfo=NEPAL)
 if hour==24:t+=dt.timedelta(days=1)
 return t.astimezone(dt.timezone.utc).isoformat()
def parse(raw,header_index=0):
 data=list(csv.reader(io.StringIO(raw.decode('utf8',errors='replace')),delimiter=';'))
 if len(data)<=header_index:raise ValueError('missing header')
 header=[s.strip() for s in data[header_index]]
 if len(header)!=len(set(header)):raise ValueError('duplicate source columns')
 result=[];last=None
 for cells in data[header_index+1:]:
  if not cells or not any(s.strip() for s in cells):continue
  if len(cells)!=len(header):raise ValueError('column count mismatch')
  row=dict(zip(header,cells));row={('YEAR' if k.casefold()=='year' else 'MONTH' if k.casefold()=='month' else 'DAY' if k.casefold()=='day' else k):v for k,v in row.items()}
  stamp=time_utc(row)
  if last is not None and stamp<=last:raise ValueError('source timestamps not strictly increasing')
  last=stamp;result.append((stamp,row))
 return result
def load(path):
 data=Path(path).read_bytes()
 if hashlib.md5(data).hexdigest()!=MD5:raise ValueError('source archive checksum mismatch')
 with zipfile.ZipFile(io.BytesIO(data)) as z:
  fn=ROOT+'CN-5470_20152017_60.csv';vn=ROOT+'CN-5470_20152017_VAL_60.csv';dn=ROOT+'CN-5470_20152017_VAL_dailysurfaceheight.csv'
  forcing=parse(z.read(fn));validation=parse(z.read(vn),3);daily=parse(z.read(dn),1)
  observations={stamp:row for stamp,row in validation};hourly=[]
  for stamp,r in forcing:
   v=observations.get(stamp);hourly.append(dict(timestamp_utc=stamp,
    features={k:scalar(r[k]) for k in ALLOWED},
    forbidden_target_related_features_excluded=sorted(FORBIDDEN),
    interface_flux_W_m2=None if v is None else dict(site1=scalar(v['Flux -10 cm']),site2=scalar(v['Flux -8cm'])),
    probe_temperature_C=None if v is None else {k:scalar(v[k]) for k in ['T1 -10cm','T2 -5cm','T3 -1cm','T4 -8cm','T5 -1cm']}))
  daily_rows=[dict(timestamp_utc=stamp,processed_lowering_m=scalar(r['ABL (m)']),processed_accumulation_m=scalar(r['ACC (m)']),surface_height_m=scalar(r['HEIGHT (m)']),interpretation='thresholded SR50 surface displacement, NOT independently measured pure ice melt') for stamp,r in daily]
  hashes={n:hashlib.sha256(z.read(n)).hexdigest() for n in [fn,vn,dn,ROOT+'CN-5470_20152017_AWS_Metadata.docx']}
 return dict(source_md5=MD5,member_sha256=hashes,hourly=hourly,daily=daily_rows,
  primary_melt_scoring_authorized=False,
  open_gates=['resolve daily timestamp/averaging convention','separate snow accumulation/debris settlement from true ice melt','resolve heat-flux sign and sensor migration before flux scoring','freeze train/evaluation groups and confirm no prior target-based model tuning'],
  label='source-aligned component data audit, no predictive performance established')
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('archive',type=Path);p.add_argument('out',type=Path);a=p.parse_args()
 if a.out.exists():raise ValueError('output exists; never overwrite observations')
 result=load(a.archive);a.out.parent.mkdir(parents=True,exist_ok=True)
 with a.out.open('x') as f:json.dump(result,f,indent=2,allow_nan=False)
 print(json.dumps(dict(hourly_rows=len(result['hourly']),daily_rows=len(result['daily']),primary_melt_scoring_authorized=False)))
if __name__=='__main__':main()
