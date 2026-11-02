"""Independent small-output arithmetic and adversarial ROI schema review.

Does not mutate owner files or launch any flow calculation. Synthetic malformed
bundles are written only into a temporary directory and removed afterward.
"""
import copy,hashlib,importlib.util,json,tempfile
from pathlib import Path
from types import SimpleNamespace
import numpy as np
BASE=Path(__file__).resolve().parents[1]
ROOT=Path(__file__).resolve().parent

def main():
 path=BASE/'roi/matched_roi.py';spec=importlib.util.spec_from_file_location('reviewed_matched_roi',path);mod=importlib.util.module_from_spec(spec)
 import sys
 sys.modules[spec.name]=mod;spec.loader.exec_module(mod)
 result=dict(reviewed_roi_source_sha256=hashlib.sha256(path.read_bytes()).hexdigest(),reviewed_phase_source_sha256=hashlib.sha256((BASE/'phase/swe2d_v5.py').read_bytes()).hexdigest(),cases=[],counterexamples={})
 summary=json.loads((BASE/'roi/results/final_physical_smoke/SUMMARY.json').read_text())
 result['dependency_sha_checks']={}
 source_lookup={q.name:q for q in [BASE/'roi/matched_roi.py',BASE/'roi/run_v5_verification.py',BASE/'phase/swe2d_v5.py',BASE/'phase/swe2d_v4_pinned.py']}
 for rec in summary['records']:
  name=rec['id'];mf=BASE/'roi/results/final_physical_smoke'/(name+'.json');manifest=json.loads(mf.read_text())
  result['dependency_sha_checks'][name]={n:hashlib.sha256(source_lookup[n].read_bytes()).hexdigest()==h for n,h in manifest['run_metadata']['source_sha256'].items()}
  try:mod.load_bundle(mf);load_status='accepted'
  except ValueError as exc:load_status=str(exc)
  with np.load(mf.with_suffix('.npz'),allow_pickle=False) as z:
   a,b,mask=z['initial_bed_m'],z['current_bed_m'],z['mask'];dx=rec['observation']['dx_m'];roi=float(((a-b)[mask]).sum()*dx**2);total=float((a-b)[1:-1,1:-1].sum()*dx**2)
  led=rec['ledgers'];row=dict(id=name,bundle_load_status=load_status,roi_recomputed_m3=roi,roi_record_difference_m3=roi-rec['observation']['geometric_melt_m3'],roi_heat_difference_m3=roi-rec['observation']['heat_implied_melt_m3'],whole_geometry_minus_phase_m3=total-led['ice_volume_lost_m3'],generated_water_minus_density_ratio_ice_m3=led['generated_water_m3']-.917*led['ice_volume_lost_m3'],latent_minus_rhoiLvolume_J=led['morphological_latent_J']-917*334000*led['ice_volume_lost_m3'],source_sha_matches_current=manifest['run_metadata']['solver_source_sha256']==result['reviewed_phase_source_sha256'])
  result['cases'].append(row)
 s=SimpleNamespace(nx=2,ny=2,dx=.1,zb=np.zeros((4,4)),ice=np.ones((4,4),dtype=bool))
 ref=mod.capture(s,roi_id='review_fixture',bounds_m=(0.,.2,0.,.2),frame='synthetic_xy',provenance={'kind':'synthetic no field evidence'})
 s.zb[1:-1,1:-1]=-.01;r=mod.observe(s,ref,physical_time_s=1.)
 md=dict(solver_version='synthetic',solver_source_sha256=hashlib.sha256(b'synthetic').hexdigest(),configuration={'mode':'synthetic'},evidence_class='conditional_numerical_run')
 with tempfile.TemporaryDirectory() as tmp:
  p=Path(tmp)/'ok.json';mod.save_bundle(p,s,ref,r,run_metadata=md);original=json.loads(p.read_text())
  for name,edit in [('nan_volume',lambda d:d['record'].update(geometric_melt_m3=float('nan'))),('wrong_bounds',lambda d:d['record'].update(bounds_m=[0.,999.,0.,999.])),('forged_field_evidence_label',lambda d:d['run_metadata'].update(evidence_class='independent_field_validation')),('nan_heat',lambda d:d['record'].update(heat_implied_melt_m3=float('nan'),geometry_heat_residual_m3=float('nan')))]:
   d=copy.deepcopy(original);edit(d);p.write_text(json.dumps(d))
   try:mod.load_bundle(p);accepted=True
   except (ValueError,TypeError,KeyError):accepted=False
   result['counterexamples'][name+'_accepted']=accepted
  p.write_text(json.dumps(original))
  a=copy.deepcopy(r);a['area_m2']=float('nan')
  try:mod.compare(a,r,relative_reporting_min_volume_m3=1e-6);accepted=True
  except (ValueError,TypeError,KeyError):accepted=False
  result['counterexamples']['nan_comparison_area_accepted']=accepted
  bad=Path(tmp)/'badconfig.json';badmd=copy.deepcopy(md);badmd['configuration']['x']=float('nan')
  try:mod.save_bundle(bad,s,ref,r,run_metadata=badmd);refused=False
  except (ValueError,TypeError,KeyError):refused=True
  result['counterexamples']['nonfinite_configuration_refused']=refused
  result['counterexamples']['refused_save_left_orphan_npz']=bad.with_suffix('.npz').exists()
 (ROOT/'INDEPENDENT_INTEGRATION_REVIEW.json').write_text(json.dumps(result,indent=2,allow_nan=False)+'\n');print(json.dumps(result,indent=2,allow_nan=False))
if __name__=='__main__':main()
