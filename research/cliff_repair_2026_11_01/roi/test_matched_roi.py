import hashlib,json,tempfile,unittest
from pathlib import Path
from types import SimpleNamespace
import numpy as np
from matched_roi import capture,observe,compare,save_bundle,load_bundle

def state(dx=.25):
 n=round(1/dx);return SimpleNamespace(dx=dx,nx=n,ny=n,zb=np.zeros((n+2,n+2)),ice=np.ones((n+2,n+2),dtype=bool))
def ref(s):return capture(s,roi_id='fixed_test',bounds_m=(0.,1.,0.,1.),frame='idealised_test_m',provenance=dict(source='synthetic test, not a glacier'))
class ROITests(unittest.TestCase):
 def fixture_bundle(self,tmp):
  s=state();r=ref(s);s.zb[1:-1,1:-1]=-.01;a=observe(s,r,physical_time_s=1.,heat_implied_volume_m3=.01)
  meta=dict(solver_version='synthetic_v5_fixture',solver_source_sha256=hashlib.sha256(b'test').hexdigest(),configuration={'test':True},evidence_class='conditional_numerical_run');p=Path(tmp)/'case.json';save_bundle(p,s,r,a,run_metadata=meta);return p
 def test_matched_volume_multiple_grids(self):
  records=[]
  for dx in [.25,.125]:
   s=state(dx);r=ref(s);s.zb[1:-1,1:-1]=-.01;records.append(observe(s,r,physical_time_s=1.,heat_implied_volume_m3=.01))
  c=compare(*records,relative_reporting_min_volume_m3=1e-5);self.assertAlmostEqual(c['signed_absolute_error_m3'],0.);self.assertAlmostEqual(c['relative_error'],0.)
 def test_ghosts_not_counted(self):
  s=state();r=ref(s);s.zb[0]=-99;s.zb[-1]=-99;s.zb[:,0]=-99;s.zb[:,-1]=-99;self.assertEqual(observe(s,r,physical_time_s=1.)['geometric_melt_m3'],0.)
 def test_mismatch_is_not_scored(self):
  s=state();r=ref(s);a=observe(s,r,physical_time_s=1.);b=dict(a,physical_time_s=2.)
  with self.assertRaises(ValueError):compare(a,b,relative_reporting_min_volume_m3=1e-5)
 def test_unaligned_and_nonice_refuse(self):
  s=state()
  with self.assertRaises(ValueError):capture(s,roi_id='a',bounds_m=(0.,.9,0.,1.),frame='m',provenance={'x':1})
  s.ice[1,1]=False
  with self.assertRaises(ValueError):ref(s)
 def test_ledger_mismatch_refuses(self):
  s=state();r=ref(s);s.zb[1:-1,1:-1]=-.01
  with self.assertRaises(ValueError):observe(s,r,physical_time_s=1.,heat_implied_volume_m3=.02)
 def test_zero_reference_absolute_only(self):
  s=state();r=ref(s);a=observe(s,r,physical_time_s=1.);s.zb[1:-1,1:-1]=-.001;b=observe(s,r,physical_time_s=1.)
  c=compare(b,a,relative_reporting_min_volume_m3=1e-5);self.assertIsNone(c['relative_error']);self.assertAlmostEqual(c['signed_absolute_error_m3'],.001)
 def test_bundle_integrity_and_no_overwrite(self):
  s=state();r=ref(s);s.zb[1:-1,1:-1]=-.01;a=observe(s,r,physical_time_s=1.,heat_implied_volume_m3=.01)
  meta=dict(solver_version='synthetic_v5_fixture',solver_source_sha256=hashlib.sha256(b'test').hexdigest(),configuration={'test':True},evidence_class='conditional_numerical_run')
  with tempfile.TemporaryDirectory() as tmp:
   p=Path(tmp)/'case.json';save_bundle(p,s,r,a,run_metadata=meta);m=load_bundle(p);self.assertAlmostEqual(m['record']['geometric_melt_m3'],.01)
   with self.assertRaises(ValueError):save_bundle(p,s,r,a,run_metadata=meta)
   data=p.with_suffix('.npz');data.write_bytes(data.read_bytes()+b'tamper')
   with self.assertRaises(ValueError):load_bundle(p)
 def test_float64_and_reference_not_mutated(self):
  s=state();r=ref(s);self.assertEqual(r.bed_m.dtype,np.float64);s.zb[1,1]=-1;self.assertEqual(r.bed_m[1,1],0.);self.assertFalse(r.bed_m.flags.writeable)
 def test_nan_comparison_refuses(self):
  s=state();a=observe(s,ref(s),physical_time_s=1.)
  for k in ['area_m2','geometric_melt_m3','physical_time_s']:
   b=dict(a,**{k:float('nan')})
   with self.assertRaises(ValueError):compare(a,b,relative_reporting_min_volume_m3=1e-5)
 def test_invalid_config_leaves_no_files(self):
  with tempfile.TemporaryDirectory() as tmp:
   s=state();r=ref(s);a=observe(s,r,physical_time_s=1.);p=Path(tmp)/'case.json';meta=dict(solver_version='v5',solver_source_sha256=hashlib.sha256(b'test').hexdigest(),configuration={'bad':float('nan')},evidence_class='conditional_numerical_run')
   with self.assertRaises(ValueError):save_bundle(p,s,r,a,run_metadata=meta)
   self.assertEqual(list(Path(tmp).iterdir()),[])
 def test_manifest_tampering_refuses(self):
  for change in ['nan_volume','bounds','heat','evidence','residual']:
   with tempfile.TemporaryDirectory() as tmp:
    p=self.fixture_bundle(tmp);m=json.loads(p.read_text())
    if change=='nan_volume':m['record']['geometric_melt_m3']=float('nan')
    elif change=='bounds':m['record']['bounds_m']=[0.,.5,0.,2.]
    elif change=='heat':m['record']['heat_implied_melt_m3']=999.
    elif change=='residual':m['record']['geometry_heat_residual_m3']=.1
    else:m['run_metadata']['evidence_class']='field_validation'
    p.write_text(json.dumps(m))
    with self.assertRaises(ValueError):load_bundle(p)
 def test_duplicate_manifest_key_refuses(self):
  with tempfile.TemporaryDirectory() as tmp:
   p=self.fixture_bundle(tmp);p.write_text('{"schema_version":1,'+p.read_text()[1:])
   with self.assertRaises(ValueError):load_bundle(p)
 def test_huge_metadata_number_refuses(self):
  s=state();a=observe(s,ref(s),physical_time_s=1.);b=dict(a,area_m2=10**1000)
  with self.assertRaises(ValueError):compare(a,b,relative_reporting_min_volume_m3=1e-5)
if __name__=='__main__':unittest.main()
