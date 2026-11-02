from pathlib import Path
import json,sys,traceback
ROOT=Path('/private/tmp/cliff-upgrade-20261101/package');sys.path.insert(0,str(ROOT))
import model
from face_geometry.face_profile import face_partition
OUT=Path('/private/tmp/cliff-upgrade-20261101/reviews')
d=json.loads((OUT/'thermal_review_cases_frozen.json').read_text());cases={x['id']:x['configuration'] for x in d['cases']};log=[]
original_remap=model.remap_columns
original_geom=model.advance_with_partition

def remap(p,total,wet,inc):
 try:return original_remap(p,total,wet,inc)
 except Exception as e:
  log.append(dict(component='column_remap',old_total=p.total_area_m2,new_total=total,difference=total-p.total_area_m2,
    wet_old=p.wet_area_m2,dry_old=p.dry_area_m2,new_wet=wet,reason=str(e),trace=traceback.format_exc()))
  raise

def geom(p,energy,floor,stage,cover):
 try:return original_geom(p,energy,floor,stage,cover)
 except Exception as e:
  part=face_partition(p,stage,cover)
  log.append(dict(component='geometry',energies=energy,class_areas=part['totals'],water_stage=stage,cover=cover,
    floor=p.floor_z_m,reason=str(e),trace=traceback.format_exc()))
  raise
model.remap_columns=remap;model.advance_with_partition=geom
for id in ['T11','T12','T16','T19']:
 try:model.forecast(cases[id])
 except Exception:pass
 if log:log[-1]['id']=id
(OUT/'thermal_failure_details.json').write_text(json.dumps(log,indent=2)+'\n')
print(json.dumps(log,indent=2))
