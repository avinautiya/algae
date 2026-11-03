import unittest
from eligibility_v2 import check,COMMON,INDEPENDENT
def fixture():return dict(**{k:True for k in COMMON},metric='final_cumulative_percentage_error')
class Tests(unittest.TestCase):
 def test_baseline_reproduction_not_independence(self):
  r=check(fixture(),claim='published_baseline_reproduction');self.assertTrue(r['metadata_gate_passed']);self.assertFalse(r['accuracy_improvement_established']);self.assertFalse(check(fixture(),claim='independent_pure_ice_melt_test')['metadata_gate_passed'])
 def test_daily_mae_does_not_reproduce_paper(self):
  f=fixture();f['metric']='daily_mae';self.assertFalse(check(f,claim='published_baseline_reproduction')['metadata_gate_passed'])
 def test_optimized_properties_not_independent(self):
  f=fixture();f.update({k:True for k in INDEPENDENT});f['pure_ice_melt_observation_verified']=True;f['parameter_calibration_excludes_test_targets']=False;self.assertFalse(check(f,claim='independent_pure_ice_melt_test')['metadata_gate_passed'])
 def test_false_string_never_passes(self):
  f=fixture();f['date_axis_verified']='false';self.assertFalse(check(f,claim='published_baseline_reproduction')['metadata_gate_passed'])
 def test_snow_lowering_not_pure_ice(self):
  f=fixture();f.update({k:True for k in INDEPENDENT});f['pure_ice_melt_observation_verified']=False;self.assertFalse(check(f,claim='independent_pure_ice_melt_test')['metadata_gate_passed'])
 def test_missing_model_mapping_refuses(self):
  f=fixture();del f['model_version_mapping_verified'];self.assertFalse(check(f,claim='published_baseline_reproduction')['metadata_gate_passed'])
if __name__=='__main__':unittest.main()
