"""Explicit-drainage regression fixtures for the four previously failing cases."""
from pathlib import Path
import copy,hashlib,json,sys,traceback,math
ROOT=Path('/private/tmp/cliff-upgrade-20261101/package');OUT=Path('/private/tmp/cliff-upgrade-20261101/reviews');sys.path.insert(0,str(ROOT))
from model import forecast
old=json.loads((OUT/'thermal_review_cases_pre_fix_frozen.json').read_text());allcases={r['id']:r['configuration'] for r in old['cases']}
configs=[]
for name in ['T11','T12','T16','T19']:
 c=copy.deepcopy(allcases[name]);c['upgrades']['meltwater_policy']=dict(mode='downstream_free_drainage',max_generated_discharge_m3_s=.001,provenance='SYNTHETIC separated downstream meltwater drain, no mixing with prescribed warm stream or backwater; drainage capacity independently imposed for implementation regression')
 configs.append(dict(id=name,configuration=c))
source=['model.py','face_thermal.py','thermal/column_partition.py'];before={s:hashlib.sha256((ROOT/s).read_bytes()).hexdigest() for s in source}
(OUT/'thermal_regression_contract_cases_frozen.json').write_text(json.dumps(dict(classification='implementation regressions under an explicit separated-drainage assumption; not new environmental evidence',original_case_file_sha256=hashlib.sha256((OUT/'thermal_review_cases_pre_fix_frozen.json').read_bytes()).hexdigest(),source_sha256=before,cases=configs),indent=2)+'\n')
results=[]
for case in configs:
 try:
  r=forecast(case['configuration']);mismatch=0.;area_mismatch=0.
  for row in r['series']:
   p=row['debris_columns_post_geometry'];s=row['port_states']['face_debris']
   if p is not None:
    profiles=[x for x in [p['wet_profile_C'],p['dry_profile_C']] if x is not None];n=len(profiles[0])
    energy=p['capacity_J_m2_K']/n*((p['wet_area_m2']*math.fsum(p['wet_profile_C']) if p['wet_profile_C'] is not None else 0.)+(p['dry_area_m2']*math.fsum(p['dry_profile_C']) if p['dry_profile_C'] is not None else 0.))
    err=abs(energy-s['debris_energy_j'])/max(1.,abs(energy));mismatch=max(mismatch,err)
    area_mismatch=max(area_mismatch,abs(p['area_reference_m2']-s['area_m2']))
    if err>1e-12 or p['area_reference_m2']!=s['area_m2']:raise RuntimeError('post-geometry column/scalar enthalpy or canonical-area mismatch')
  scale=max(1.,abs(r['totals']['atmospheric_J'])+r['totals']['water_in_J'])
  results.append(dict(id=case['id'],status='completed',checks=r['checks'],max_normalized_post_geometry_enthalpy_mismatch=mismatch,max_canonical_area_mismatch_m2=area_mismatch,normalized_global_energy_residual=abs(r['checks']['global_energy_residual_J'])/scale,
   steps=len(r['series']),max_actual_height_bands=max(row['profile_height_band_count'] or 0 for row in r['series'])))
 except Exception as e:results.append(dict(id=case['id'],status='failure',reason=repr(e),trace=traceback.format_exc()))
after={s:hashlib.sha256((ROOT/s).read_bytes()).hexdigest() for s in source}
out=dict(results=results,source_sha256_before=before,source_sha256_after=after,source_changed_during_review=before!=after)
(OUT/'thermal_regression_contract_results.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(out,indent=2))
