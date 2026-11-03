import math,unittest
from neutral_profile import *
class NeutralProfileTests(unittest.TestCase):
 def test_exact_log_profile(self):
  vals=[3+2*math.log(z)for z in HEIGHTS];b,weights,r=log_slope(vals);self.assertAlmostEqual(b,2);self.assertAlmostEqual(r,1);self.assertAlmostEqual(sum(weights),0)
 def test_constant_T_RH_gives_zero_heat_and_vapor(self):
  r=neutral_flux([1+math.log(z)for z in HEIGHTS],[10]*4,[30]*4,68000)
  self.assertAlmostEqual(r['H_up_W_m2'],0,places=9);self.assertAlmostEqual(r['LE_up_W_m2'],0,places=9)
 def test_warm_moist_lower_air_flux_up(self):
  r=neutral_flux([1+math.log(z)for z in HEIGHTS],[15-math.log(z)for z in HEIGHTS],[40]*4,68000)
  self.assertGreater(r['H_up_W_m2'],0);self.assertGreater(r['LE_up_W_m2'],0)
 def test_reject_wind_inversion_as_supported_neutral(self):
  r=neutral_flux([4,3,2,1],[10]*4,[30]*4,68000);self.assertIsNone(r['H_up_W_m2']);self.assertTrue(r['status'].startswith('unsupported'))
 def test_pressure_and_RH_validation(self):
  with self.assertRaises(ValueError):specific_humidity(10,120,68000)
  with self.assertRaises(ValueError):specific_humidity(10,30,float('nan'))
if __name__=='__main__':unittest.main()
