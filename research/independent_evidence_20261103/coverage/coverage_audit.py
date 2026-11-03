"""Calendar-based support audit; never fits or scores held-out sensor values."""
import argparse,datetime as dt,hashlib,json,math
from pathlib import Path
NEPAL=dt.timezone(dt.timedelta(hours=5,minutes=45))
TRAIN0=dt.datetime(2015,12,1,tzinfo=NEPAL)
TEST0=dt.datetime(2015,12,29,tzinfo=NEPAL)
END=dt.datetime(2016,1,12,tzinfo=NEPAL)
def finite(v):
 try:return type(v) in (float,int) and math.isfinite(v)
 except OverflowError:return False
def audit(rows):
 out={site:{period:[] for period in ['training','evaluation']} for site in ['site1','site2']};n={site:{period:0 for period in ['training','evaluation']} for site in out};seen=set();last=None
 for r in rows:
  t=dt.datetime.fromisoformat(r['timestamp_utc'])
  if t.tzinfo is None or t.utcoffset() is None:raise ValueError('explicit timestamp offset required')
  if t in seen or (last is not None and t<=last):raise ValueError('duplicate or unordered source records')
  seen.add(t);last=t
  if not TRAIN0<=t<END:continue
  period='training' if t<TEST0 else 'evaluation'
  temp=r.get('probe_temperature_C');flux=r.get('interface_flux_W_m2')
  for site,keys in [('site1',['T1 -10cm','T2 -5cm','T3 -1cm']),('site2',['T4 -8cm','T5 -1cm'])]:
   n[site][period]+=1
   if isinstance(temp,dict) and isinstance(flux,dict) and all(finite(temp.get(k)) for k in keys) and finite(flux.get(site)):
    out[site][period].append(t.astimezone(dt.timezone.utc).isoformat())
 records={s:{p:dict(source_records=n[s][p],complete_records=len(ids),support_sha256=hashlib.sha256(json.dumps(ids,separators=(',',':')).encode()).hexdigest()) for p,ids in periods.items()} for s,periods in out.items()}
 return dict(calendar_timezone='UTC+05:45',training_start=TRAIN0.isoformat(),evaluation_start=TEST0.isoformat(),evaluation_end_exclusive=END.isoformat(),support=records,heldout_values_scored=False,melt_validation_established=False,scope='nonperformance availability audit; sign, intervals, migration and model applicability require separate review')
def main():
 p=argparse.ArgumentParser();p.add_argument('aligned',type=Path);p.add_argument('out',type=Path);a=p.parse_args()
 if a.out.exists():raise ValueError('output exists')
 raw=a.aligned.read_bytes();data=json.loads(raw);result=audit(data['hourly']);result['aligned_source_sha256']=hashlib.sha256(raw).hexdigest();a.out.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps(result,indent=2))
if __name__=='__main__':main()
