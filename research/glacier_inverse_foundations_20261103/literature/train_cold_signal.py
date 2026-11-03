"""Existing Tower1 TRAIN-temperature eligibility only; no camera or melt access."""
import argparse,csv,datetime as dt,hashlib,json
from pathlib import Path

def audit(path):
 path=Path(path);a=dt.datetime(2022,12,2,3);b=dt.datetime(2022,12,16,3)
 with path.open() as stream:rows=[r for r in csv.DictReader(stream) if a<=dt.datetime.fromisoformat(r['SOURCE_LOCAL_LABEL'])<b]
 result={}
 for col in ['DT1','DT4','DT5']:
  v=[float(r[col]) for r in rows]
  result[col]=dict(rows=len(v),minimum_C=min(v),maximum_C=max(v),negative_rows=sum(x<0 for x in v),rows_below_minus04C=sum(x<-.4 for x in v),rows_within_plusminus04C=sum(abs(x)<=.4 for x in v))
 return dict(source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),start=a.isoformat(),end_exclusive=b.isoformat(),sensor_accuracy_C=.4,diagnostics=result,ice_side_temperature_independently_measured=False,field_reference_excitation_established=False,camera_files_read=False)
if __name__=='__main__':
 p=argparse.ArgumentParser();p.add_argument('process');a=p.parse_args();print(json.dumps(audit(a.process),indent=2))
