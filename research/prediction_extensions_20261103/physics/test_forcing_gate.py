import unittest,math
from forcing_gate import screen,water_ice_flux,RHOI,LF

def config():
 return dict(required_gap_repair_m=.25-(.07/((2/3)**1.5*math.sqrt(9.81)*.5))**(2/3),horizon_s=21*86400.,pulse_discharge_m3_s=.14,bypass_area_m2=1.,toe_h_upper_W_m2_K=200.,toe_temperature_upper_K=.5,bypass_h_lower_W_m2_K=200.,bypass_temperature_lower_K=.5,additional_toe_ice_flux_upper_W_m2=0.,additional_bypass_ice_flux_lower_W_m2=0.,evidence=dict(bound_sources=['synthetic verification only'],equal_rating_controls_verified=True,bounds_hold_over_window=True))
class Tests(unittest.TestCase):
 def test_reproduces_no_extra_heat_bound(self):
  r=screen(**config());self.assertTrue(r['ruled_out_within_window']);self.assertAlmostEqual(r['maximum_repair_m'],(100-water_ice_flux(.07,1,200,.5))*21*86400/(RHOI*LF))
 def test_extra_heat_not_proof_of_success(self):
  c=config();c['additional_toe_ice_flux_upper_W_m2']=20.;r=screen(**c);self.assertFalse(r['ruled_out_within_window']);self.assertIn('does not predict',r['interpretation'])
 def test_zero_flow_and_heat(self):self.assertEqual(water_ice_flux(0,1,200,.5),0.);self.assertEqual(water_ice_flux(.1,1,0,.5),0.)
 def test_flux_never_exceeds_film_or_supply(self):
  for q in [1e-8,1e-4,.1,1]:
   f=water_ice_flux(q,1,200,.5);self.assertLessEqual(f,100+1e-12);self.assertLessEqual(f,1000*4217*q*.5+1e-12)
 def test_bounds_contract_required(self):
  c=config();c['evidence']={}
  with self.assertRaises(ValueError):screen(**c)
 def test_refuses_invalid_units_or_numbers(self):
  for k,v in [('horizon_s',0),('bypass_area_m2',-1),('toe_temperature_upper_K',float('nan')),('pulse_discharge_m3_s',True)]:
   c=config();c[k]=v
   with self.assertRaises(ValueError):screen(**c)
if __name__=='__main__':unittest.main()
