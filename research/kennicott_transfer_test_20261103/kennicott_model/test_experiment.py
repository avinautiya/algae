import unittest
from unittest.mock import patch
import experiment as e
class CalibrationTests(unittest.TestCase):
 def fake(self,plan,stake,method,C,**kwargs):return {'physical_melt_mmwe':C/10000}
 def test_one_physical_parameter_bisection_and_no_rescaling(self):
  with patch.object(e,'forward',self.fake):r=e.calibrate_capacity({},e.PRIMARY,150.)
  self.assertAlmostEqual(r['capacity_J_m3K'],1.5e6,delta=20)
  self.assertFalse(r['physical_energy_ledger_rescaled']);self.assertFalse(r['source_mass_density_reconverted'])
 def test_refuses_outside_bracket(self):
  with patch.object(e,'forward',self.fake):
   with self.assertRaises(e.UnsupportedModel):e.calibrate_capacity({},e.PRIMARY,500.)
 def test_endpoint_is_allowed(self):
  with patch.object(e,'forward',self.fake):r=e.calibrate_capacity({},e.PRIMARY,80.)
  self.assertEqual(r['capacity_J_m3K'],.8e6)
 def test_detects_nonmonotone_inside_bracket(self):
  def bad(plan,stake,method,C,**kwargs):return {'physical_melt_mmwe':C/10000 if C in[.8e6,3.5e6]else 1000}
  with patch.object(e,'forward',bad):
   with self.assertRaises(e.UnsupportedModel):e.calibrate_capacity({},e.PRIMARY,150.)
 def test_pointmass_does_not_get_density_converted(self):
  self.assertFalse(e.SPEC['mass_target_density_reconversion'])
class StructuralForwardTests(unittest.TestCase):
 def test_P12_lower_probe_is_not_imposed_inside_stake(self):
  plan={'sites':{'P12':{'inverse_profile':'profile','forcing_plane_m':.09,'stake_thickness_m':.29}},'profiles':{'profile':{'alpha_layered_m2_s':[1e-6,2e-6],'alpha_homogeneous_m2_s':1e-6,'depths_m':[.165,.24,.315]}}}
  boundary={'metadata':{'blocked_RMSE_C':1.},'boundary_C':[5.,6.],'proxy_boundary_C':[4.,5.],'is_proxy':[False,False],'times_utc':['2023-01-01T00:00:00Z','2023-01-01T01:00:00Z']}
  captured={}
  def fake(**kw):
   captured.update(kw);n=len(kw['surface_temperature_C']);return {'physical_melt_mmwe':0.,'trace':[{'melt_mmwe':0.}for _ in range(n)]}
  with patch.object(e,'boundary_for',lambda *args:boundary),patch.object(e,'simulate',fake):r=e.forward(plan,'P12',e.PRIMARY,1.5e6)
  self.assertAlmostEqual(captured['debris_depth_m'],.20);self.assertEqual(len(captured['surface_temperature_C']),12)
  self.assertTrue(all(5<x<=6 for x in captured['surface_temperature_C']))
  self.assertTrue(r['no_below_ice_thermistor_forcing']);self.assertEqual(r['stake_ice_depth_m'],.29)
 def test_geometry_corrected_dose_baselines(self):
  plan={'sites':{'P12':{'forcing_plane_m':.09,'stake_thickness_m':.29}}};b={'boundary_C':[2,4],'air_C':[1,3]}
  with patch.object(e,'boundary_for',lambda *args:b):d=e.index_doses(plan,'P12')
  self.assertEqual(d['buried_temp_exposure'],10800.);self.assertAlmostEqual(d['conventional_gradient'],54000.)
  self.assertAlmostEqual(d['depth_corrected_air_pdd'],36000.);self.assertEqual(d['constant_rate'],3600.)
if __name__=='__main__':unittest.main()

class ExactLayerResistanceTests(unittest.TestCase):
 def test_cell_crossing_preserves_total_series_resistance(self):
  for n in[8,16,32]:
   alpha=e.layer_cell_alpha(.09,.29,.24,[.4e-6,.9e-6],n);dz=.2/n
   self.assertAlmostEqual(sum(dz/a for a in alpha),.15/.4e-6+.05/.9e-6,places=7)
