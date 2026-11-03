import unittest
from resolve_target import error_percent,density_sensitivity
class Tests(unittest.TestCase):
 def data(self):return dict(model_sum_m_water_equivalent=.7186453836834981,thresholded_surface_lowering_m=1.0233947,published_Table3_percent=-23.3,published_score_reproduced=False)
 def test_protocol_result(self):self.assertAlmostEqual(error_percent(.7186453836834981,1.0233947,900),-21.97586417445369)
 def test_apparent_match_cannot_authorize_density(self):
  a=density_sensitivity(self.data());self.assertEqual(round(a['sensitivity_only'][1]['error_percent'],1),-23.3);self.assertFalse(a['published_score_reproduced_under_authorized_conditional_rule']);self.assertFalse(a['independent_melt_validation'])
 def test_implied_target_not_observation(self):self.assertTrue(density_sensitivity(self.data())['implied_target_is_not_an_observation'])
 def test_bad_quantities(self):
  for x in [0.,-1.,float('nan'),True]:
   with self.assertRaises(ValueError):error_percent(x,1.,900.)
if __name__=='__main__':unittest.main()
