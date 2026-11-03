import unittest
from training_capacity_audit import compute
class Checks(unittest.TestCase):
 def test_training_only_bounds(self):
  r=compute();self.assertEqual(r['rows'],672);self.assertEqual(r['all_three_valid_rows'],672)
  c=r['finite_sensible_capacity_scenes'];v=[x['C_J_m3K'] for x in c.values()]
  self.assertTrue(v[0]<v[1]<v[2]);self.assertAlmostEqual(v[2]-v[0],.44*(4180000-1206),6)
 def test_not_field_capacity_or_phase_validation(self):
  r=compute();self.assertFalse(r['melt_validation_ready']);self.assertIn('not a volume-weighted',r['proxy_definition'])
  self.assertEqual(r['maximum_pore_phase_latent_J_m3_assuming_water_volume_phi'],146960000)
if __name__=='__main__':unittest.main()
