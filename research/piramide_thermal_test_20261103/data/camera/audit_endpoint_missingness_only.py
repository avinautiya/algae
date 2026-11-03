"""Coverage-only audit: never emits finite camera distance or ablation values."""
from pathlib import Path
import hashlib,json
import numpy as np,pandas as pd
from matio import load_from_mat
p=Path(__file__).parent;source=p/'Ablation_Tower1.mat'
m=json.loads((p/'SEALED_SOURCE_MANIFEST.json').read_text());expected=next(x['sha256'] for x in m['members'] if x['file']==source.name)
assert hashlib.sha256(source.read_bytes()).hexdigest()==expected
raw=load_from_mat(source,raw_data=True)['tracks'].properties['any']
names=[str(x.item()) for x in raw['varNames'][0,0].ravel()];assert names==['distance','ablation']
epoch=np.asarray(raw['rowTimes'][0,0].properties['data']).ravel();times=pd.DatetimeIndex(epoch.astype('int64').astype('datetime64[ms]'))
cells=raw['data'][0,0].ravel();native_distance=np.asarray(cells[0]).ravel()
converted=load_from_mat(source,raw_data=False)['tracks']
assert np.array_equal(np.isfinite(native_distance),np.isfinite(converted['distance'].to_numpy()))
old=pd.Timestamp('2022-12-30T03:00');start=pd.Timestamp('2022-12-16T03:00')
indices=np.flatnonzero(np.abs((times-old).total_seconds())<=12*3600)
rows=[{'source_local_label':times[i].isoformat(),'distance_isfinite':bool(np.isfinite(native_distance[i])),'distance_isnan':bool(np.isnan(native_distance[i]))} for i in indices]
valid=[i for i in indices if np.isfinite(native_distance[i])]
assert valid,'No finite endpoint within predeclared12h search'
selected=min(valid,key=lambda i:(abs((times[i]-old).total_seconds()),times[i]))
checks=[]
for label in [start,old,times[selected]]:
 i=times.get_indexer([label])[0];assert i>=0
 checks.append({'source_local_label':label.isoformat(),'distance_isfinite':bool(np.isfinite(native_distance[i])),'distance_isnan':bool(np.isnan(native_distance[i]))})
out={'source_sha256':expected,'original_start':start.isoformat(),'original_end':old.isoformat(),'selected_end':times[selected].isoformat(),'selection_rule':'Nearest finite native camera distance endpoint within +/-12h of originalDec30T03; minimize absolute clock difference; tieearlier; no interpolation; training unchanged','endpoint_checks':checks,'nearby_timestamps_finiteness_only':rows,'raw_MCOS_and_converted_finiteness_equal':True,'missing_value_native_nan_confirmed':checks[1]['distance_isnan'],'new_target_values_exposed':False,'prediction_error_computed':False,'scope':'Coverage-only revision; publish all amended predictions before releasing selected observed value'}
(p/'MISSINGNESS_ONLY_ENDPOINT_AUDIT.json').write_text(json.dumps(out,indent=2)+'\n')
print('Finiteness-only audit saved. Selected native end label:',out['selected_end'])
