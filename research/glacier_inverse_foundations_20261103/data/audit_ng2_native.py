"""Thermal source diagnostics only; no ice-ablation targets exist in this file."""
from pathlib import Path
import json,hashlib
import numpy as np,pandas as pd
P=Path(__file__).parent;f=P/'NG2_Td_Dec14-Apr16.txt';r=json.loads((P/'NG2_RECORD.json').read_text());m=r['files'][0]
assert len(f.read_bytes())==m['size'];assert 'md5:'+hashlib.md5(f.read_bytes()).hexdigest()==m['checksum']
a=np.loadtxt(f);assert a.shape[1]==17
idx=pd.DatetimeIndex(pd.to_datetime(pd.DataFrame(a[:,:6].astype(int),columns=['year','month','day','hour','minute','second'])))
assert idx.is_unique and idx.is_monotonic_increasing
z=np.array([.01,.2,.4,.6,.8,1.,1.2,1.4,1.6,1.8,2.]);T=a[:,6:];cols=[f'T_{x:g}m_C' for x in z]
frame=pd.DataFrame(T,index=idx,columns=cols);frame.to_csv(P/'NG2_NATIVE_THERMAL.csv',index_label='SOURCE_LOCAL_LABEL')

def runs(mask):
 result=[];i=0
 while i<len(mask):
  if not mask[i]:i+=1;continue
  j=i+1
  while j<len(mask) and mask[j] and (idx[j]-idx[j-1]).total_seconds()<=21601.0:j+=1
  result.append({'start':idx[i].isoformat(),'end':idx[j-1].isoformat(),'samples':j-i,'duration_days':(idx[j-1]-idx[i]).total_seconds()/86400})
  i=j
 return sorted(result,key=lambda x:x['duration_days'],reverse=True)
subzero=np.all(T[:,z>=1.] < -1.,axis=1)
report={'source_doi':'10.5281/zenodo.3935685','source_sha256':hashlib.sha256(f.read_bytes()).hexdigest(),'md5_verified':True,
 'rows':len(idx),'first_native_label':idx[0].isoformat(),'last_native_label':idx[-1].isoformat(),'timezone':'Not specified by source record; no UTC claim',
 'depths_m':list(z),'source_thickness_m':2.,'lowermost_probe_source_designation':'Installed at debris-ice interface; source metadata, not an independent borehole ice profile',
 'sensor_accuracy_C':.25,'accuracy_is_not_sigma':True,'step_seconds_counts':{str(k):int(v) for k,v in zip(*np.unique(np.diff(idx.values).astype('timedelta64[s]').astype(int),return_counts=True))},
 'nonfinite_temperature_values':int((~np.isfinite(T)).sum()),'minimum_temperature_by_depth_C':dict(zip(cols,T.min(axis=0).tolist())),
 'maximum_temperature_by_depth_C':dict(zip(cols,T.max(axis=0).tolist())),
 'interface_below_minus1_count':int((T[:,-1]<-1).sum()),'interface_below_minus_accuracy_count':int((T[:,-1]<-.25).sum()),
 'deep_domain_1to2m_all_probes_below_minus1_longest_runs':runs(subzero)[:5],
 'melt_or_flux_targets_present':False,'independent_bulk_C_present':False,
 'fresh_physical_test':'Lowerdebris1-2m pluscoldice response; measuredT1m topboundary, train/evalinterior+interface temperatures; avoidclaimingmeltaccuracywithoutindependentmelt.',
 'limitations':['6h sampling cannotresolvehourlydiurnalflux','Noiceprobe belowinterface inrecord','NoindependentabsoluteC/heatflux','Sourcepaper flags NG2 near-surfacebias; lower-domainboundaryavoidsusingthatprobe']}
(P/'NG2_THERMAL_ELIGIBILITY.json').write_text(json.dumps(report,indent=2)+'\n')
print('NG2 native source verified',len(idx),'rows',idx[0],idx[-1]);print('Longest fullycold1-2m intervals:',report['deep_domain_1to2m_all_probes_below_minus1_longest_runs'][:3])
