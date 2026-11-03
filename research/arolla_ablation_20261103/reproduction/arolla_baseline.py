"""Verified published Arolla predictions; no observed melt is loaded or scored."""
import csv,datetime as dt,hashlib,json,math,re
from pathlib import Path
ROOT=Path(__file__).resolve().parent
LOCAL=dt.timezone(dt.timedelta(hours=2))
EXPECTED_START=dt.datetime(2010,7,4,12,tzinfo=LOCAL)
EXPECTED_END=dt.datetime(2010,9,12,17,tzinfo=LOCAL)
def timestamp(text):
 m=re.fullmatch(r'(\d{2})-([A-Za-z]{3})-(\d{4}) (\d{2}):(\d{2}):(\d{2})',text.strip())
 if not m:raise ValueError('unexpected published timestamp format')
 day,month,year,h,minute,second=m.groups();months=['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
 return dt.datetime(int(year),months.index(month)+1,int(day),int(h),int(minute),int(second),tzinfo=LOCAL)
def load(root=None):
 root=Path(root or ROOT/'arolla_outputs');manifest=json.loads((ROOT/'AROLLA_OUTPUT_MANIFEST.json').read_text())
 for e in manifest['members']:
  p=root/e['member']
  if hashlib.sha256(p.read_bytes()).hexdigest()!=e['sha256']:raise ValueError('published source hash mismatch')
 labels=[timestamp(t) for t in (root/'Dates_Rows/Dates_hourly_ARO.csv').read_text().splitlines()]
 expected=[EXPECTED_START+dt.timedelta(hours=i) for i in range(1686)]
 if labels!=expected or labels[-1]!=EXPECTED_END:raise ValueError('source-prescribed hourly axis mismatch')
 matrix=[[float(x) for x in r] for r in csv.reader((root/'DEBCF/ARO/DEBCF_ARO_EXP1_EXP2.csv').read_text().splitlines())]
 if len(matrix)!=len(labels) or any(len(r)!=13 or not all(math.isfinite(x) and x>=0 for x in r) for r in matrix):raise ValueError('invalid published melt matrix')
 return [dict(label_local=t.isoformat(),label_utc=t.astimezone(dt.timezone.utc).isoformat(),model_melt_mwe_per_step=r[0]) for t,r in zip(labels,matrix)]
def map_window(rows,start,end,*,label_convention):
 """Only construct row mapping; never aggregate observed or model test totals.

 Both interpretations must be retained if the authoritative label convention
 is unknown. Picking one because it improves score is not allowed.
 """
 if label_convention not in {'interval_ending','interval_starting'}:raise ValueError('explicit step-label interpretation required')
 if not isinstance(start,dt.datetime) or not isinstance(end,dt.datetime) or start.utcoffset() is None or end.utcoffset() is None or end<=start:raise ValueError('ordered offset-aware endpoints required')
 if start.minute or start.second or start.microsecond or end.minute or end.second or end.microsecond:raise ValueError('fractional-hour visits require a separate observation operator')
 included=[];intervals=[]
 for i,r in enumerate(rows):
  label=dt.datetime.fromisoformat(r['label_local']);a=label-dt.timedelta(hours=1) if label_convention=='interval_ending' else label;b=a+dt.timedelta(hours=1)
  if a>=start and b<=end:included.append(i);intervals.append((a,b))
 if not intervals or intervals[0][0]!=start or intervals[-1][1]!=end or any(a[1]!=b[0] for a,b in zip(intervals,intervals[1:])):raise ValueError('window not fully supported')
 return dict(start_local=start.astimezone(LOCAL).isoformat(),end_local=end.astimezone(LOCAL).isoformat(),label_convention=label_convention,row_indices_zero_based=included,hours=len(included),prediction_values_aggregated=False,observations_used=False)
def describe():
 rows=load();start=dt.datetime(2010,8,28,tzinfo=LOCAL);end=dt.datetime(2010,9,11,tzinfo=LOCAL)
 return dict(model='DEBCF',actual_site_debris_thickness_cm=6,actual_thickness_column_one_based=1,unit='m_water_equivalent_per_step',run_rows=len(rows),run_start=rows[0]['label_local'],run_end=rows[-1]['label_local'],utc_offset_hours=2,source_forcing_indices_one_based=[1030,2715],example_window_mappings=[map_window(rows,start,end,label_convention=c) for c in ['interval_ending','interval_starting']],metadata_gates=['actual stake visit times uncertain; enumerate independently declared visit-time bounds','source step-label convention unresolved; retain both interpretations until justified','source k provenance blank in site metadata; do not call every supplied property independently measured'],observed_targets_loaded=False,test_scores_computed=False)
if __name__=='__main__':
 p=ROOT/'AROLLA_BASELINE_MAPPING.json'
 if p.exists():raise ValueError('output exists')
 result=describe();p.write_text(json.dumps(result,indent=2)+'\n');print(json.dumps({k:v for k,v in result.items() if k!='example_window_mappings'},indent=2))
