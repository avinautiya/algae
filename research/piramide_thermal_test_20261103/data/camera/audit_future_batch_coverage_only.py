"""Source eligibility only. Never emits future camera outcome values."""
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
from matio import load_from_mat
import sys
sys.path.insert(0,str(Path(__file__).parents[1]))
from read_verified_process import read_process
P=Path(__file__).parent;source=P/'Ablation_Tower1.mat'
manifest=json.loads((P/'SEALED_SOURCE_MANIFEST.json').read_text());expected=next(x['sha256'] for x in manifest['members'] if x['file']==source.name)
assert hashlib.sha256(source.read_bytes()).hexdigest()==expected
obj=load_from_mat(source,raw_data=True)['tracks'].properties['any']
names=[str(x.item()) for x in obj['varNames'][0,0].ravel()];assert names==['distance','ablation']
timestamps=np.asarray(obj['rowTimes'][0,0].properties['data']).ravel();times=pd.DatetimeIndex(timestamps.astype('int64').astype('datetime64[ms]'))
# Only finiteness is consumed; native numeric values are never returned or scored.
finite=np.isfinite(np.asarray(obj['data'][0,0].ravel()[0]).ravel())
nominals=['2023-01-13T03:00','2023-01-27T03:00','2023-02-10T03:00','2023-02-24T03:00']
ends=[]
for nominal in nominals:
 t=pd.Timestamp(nominal);ids=np.flatnonzero(np.abs((times-t).total_seconds())<=43200)
 valid=[i for i in ids if finite[i]]
 chosen=min(valid,key=lambda i:(abs((times[i]-t).total_seconds()),times[i])) if valid else None
 exact=times.get_indexer([t])[0]
 ends.append({'nominal_end':t.isoformat(),'nominal_is_native':bool(exact>=0),'nominal_isfinite':bool(finite[exact]) if exact>=0 else None,
    'selected_end':times[chosen].isoformat() if chosen is not None else None,
    'nearby_time_finiteness':[{'source_local_label':times[i].isoformat(),'isfinite':bool(finite[i])} for i in ids]})
thermal=read_process(P.parent/'decoded','debtemp','Tower1');air=read_process(P.parent/'decoded','trh','Tower1')

def coverage(frame,columns,start,end):
 s,e=pd.Timestamp(start),pd.Timestamp(end)
 frame=frame[columns];valid=frame.notna().all(axis=1)&np.isfinite(frame).all(axis=1)
 labels=frame.index[valid];inside=labels[(labels>=s)&(labels<=e)]
 left=labels[labels<=s];right=labels[labels>=e]
 bounded=bool(len(left) and len(right))
 bridge=list(inside)
 if len(left):bridge.append(left[-1])
 if len(right):bridge.append(right[0])
 ordered=pd.DatetimeIndex(sorted(set(bridge)))
 if len(ordered)>1:
  gaps=np.diff(ordered.values).astype('timedelta64[s]').astype(float)/3600
  relevant=[g for g,a,b in zip(gaps,ordered[:-1],ordered[1:]) if b>s and a<e]
  maxgap=max(relevant,default=0.)
 else:maxgap=None
 start_gap=max(0.,(labels[0]-s).total_seconds()/3600) if len(labels) else None
 end_gap=max(0.,(e-labels[-1]).total_seconds()/3600) if len(labels) else None
 if not bounded:maxgap=max([v for v in [maxgap,start_gap,end_gap] if v is not None],default=None)
 return {'columns':columns,'native_finite_inside_rows':len(inside),'source_bounds_cover_interval':bounded,'maximum_gap_hours':maxgap,
    'coverage_pass_maxgap18h_and_no_extrapolation':bool(bounded and maxgap is not None and maxgap<=18),
    'source_first_label':labels[0].isoformat(),'source_last_label':labels[-1].isoformat()}

intervals=[];start='2022-12-30T00:00:00'
for x in ends:
 end=x['selected_end']
 if end is None:
  intervals.append({'start':start,'nominal_end':x['nominal_end'],'end':None,'source_eligible':False,'reason':'No finite camera endpoint withinfixed12hsearch'});start=x['nominal_end'];continue
 a=coverage(thermal,['DT1','DT2','DT3','DT4','DT5'],start,end);b=coverage(air,['TA_3'],start,end)
 intervals.append({'start':start,'nominal_end':x['nominal_end'],'end':end,'thermal':a,'air_temperature_2m':b,
   'source_eligible':a['coverage_pass_maxgap18h_and_no_extrapolation'] and b['coverage_pass_maxgap18h_and_no_extrapolation']})
 start=end
out={'scope':'FUTURE camera finiteness and process-feature coverage ONLY; no future outcome values/statistics/errors',
 'source_camera_sha256':expected,'source_naive_local_clock_only':True,
 'selection_rule':'For each fixed nominal14dayanchor, nearestfinite native camera distance timestamp within±12h; tieearlier; no interpolation; no performance selection',
 'known_batch_start':'2022-12-30T00:00:00','future_anchors':ends,'proposed_intervals':intervals,
 'new_target_values_exposed':False,'future_melt_predictions_scored':False}
(P/'FUTURE_BATCH_COVERAGE_ONLY.json').write_text(json.dumps(out,indent=2)+'\n')
print('Future eligibility audit saved; timestamps/coverage only.')
for r in intervals:print(r['start'],r['end'],r['source_eligible'],r.get('air_temperature_2m'))
