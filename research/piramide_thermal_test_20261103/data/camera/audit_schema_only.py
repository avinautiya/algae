"""Inspect ONLY native schema and timestamp arrays. Never dereference data cells."""
import json,hashlib
from pathlib import Path
import numpy as np
from matio import load_from_mat
p=Path(__file__).parent;rows=[]
for j in [1,2,3]:
 path=p/f'Ablation_Tower{j}.mat';d=load_from_mat(path,raw_data=True)
 o=d['tracks'];a=o.properties['any'];assert o.classname=='timetable'
 def prop(k):return a[k][0,0]
 names=[str(x.item()) for x in prop('varNames').ravel()]
 units=[str(x.item()) if np.asarray(x).size else '' for x in prop('varUnits').ravel()]
 rt=prop('rowTimes');v=np.asarray(rt.properties['data']).ravel();assert np.isfinite(v).all() and not np.iscomplexobj(v)
 assert (v==np.rint(v)).all();times=v.astype('int64').astype('datetime64[ms]')
 tz=rt.properties.get('tz',rt.properties.get('tmz'));assert np.asarray(tz).size==0
 # Timestamp-only exports are allowed. No camera column is exported here.
 timesfile=p/f'Tower{j}_NATIVE_TIMES_ONLY.json';timesfile.write_text(json.dumps([str(t) for t in times],indent=2)+'\n')
 rows.append({'tower':j,'source_sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'variable_names':names,
  'variable_units':units,'numRows':int(prop('numRows').item()),'numVars':int(prop('numVars').item()),
  'first_source_local_label':str(times[0]),'last_source_local_label':str(times[-1]),'native_timezone':'',
  'timestamps_unique':bool(len(np.unique(times))==len(times)), 'timestamps_increasing':bool(np.all(np.diff(v)>0)),
  'timestamps_file':timesfile.name,'target_data_cells_dereferenced':False})
(p/'SCHEMA_AND_CLOCK_ONLY.json').write_text(json.dumps(rows,indent=2)+'\n')
print(json.dumps(rows,indent=2))
