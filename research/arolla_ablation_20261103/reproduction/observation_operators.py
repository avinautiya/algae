"""Inventory distinct measurement operators without selecting by model fit."""
import csv,datetime as dt,hashlib,io,json,zipfile
from pathlib import Path

def inventory(path):
 b=Path(path).read_bytes()
 if hashlib.md5(b).hexdigest()!='4797deeda5c93a1895f08efb42183086':raise ValueError('archive checksum mismatch')
 with zipfile.ZipFile(io.BytesIO(b)) as z:
  raw=z.read('CN-5470_20152017/CN-5470_20152017_VAL_dailysurfaceheight.csv')
  daily=list(csv.reader(io.StringIO(raw.decode()),delimiter=';'))
  rows=[r for r in daily[2:] if dt.date(2016,4,16)<=dt.date(*map(int,r[:3]))<=dt.date(2016,7,24)]
  hourly_raw=z.read('CN-5470_20152017/CN-5470_20152017_VAL_60.csv')
  hourly=list(csv.reader(io.StringIO(hourly_raw.decode('latin1')),delimiter=';'))
  first=next(r for r in hourly[4:] if r[:5]==['2016','4','16','1','0'])
  last=next(r for r in hourly[4:] if r[:5]==['2016','7','24','23','0'])
 return dict(source_sha256=hashlib.sha256(b).hexdigest(),calendar_rows=len(rows),
  thresholded_positive_daily_lowering_sum_m=sum(float(r[4]) for r in rows),
  thresholded_positive_daily_height_gain_sum_m=sum(float(r[5]) for r in rows),
  daily_HEIGHT_endpoint_lowering_m=-(float(rows[-1][6])-float(rows[0][6])),
  daily_DIST_endpoint_increase_m=float(rows[-1][3])-float(rows[0][3]),
  hourly_DIST_labelled_endpoint_increase_m=float(last[-1])-float(first[-1]),
  raw_hourly_start_label_nepal='2016-04-16 01:00',raw_hourly_end_label_nepal='2016-07-24 23:00',
  operators_are_not_interchangeable=True,selected_operator_by_score=None,
  pure_ice_melt_authorized=False,
  unresolved='Daily means, thresholded positive changes and hourly endpoint differences have different supports and snow/settlement sensitivity. A source-defined analysis operator is needed; agreement with a model/table cannot select it.')
if __name__=='__main__':
 import argparse
 p=argparse.ArgumentParser();p.add_argument('archive',type=Path);p.add_argument('out',type=Path);a=p.parse_args()
 if a.out.exists():raise ValueError('output exists')
 r=inventory(a.archive);a.out.write_text(json.dumps(r,indent=2)+'\n');print(json.dumps(r,indent=2))
