import unittest
from storage_flux import *
class StorageFluxTests(unittest.TestCase):
 def test_p1_storage_independent_polygon_integral(self):
  h=.12;g=.02;T4=4.;T5=1.;C=2e6;top=(T4+T5)/2
  self.assertAlmostEqual(slab_sensible_J_m2(T4,T5,h,g,C),C*((h/2)*(top+T5)/2+g*T5/2))
 def test_signed_gradient_not_clipped(self):self.assertEqual(gradient_integral([0,10],[0,0],[1,1],.1),-100)
 def test_exact_ledger_and_train_recovery(self):
  y=10.;J=[100,200,300];H=[200,300,-100];r=predict(y,J,H);k=r['calibrated_conductivity_W_mK']
  self.assertAlmostEqual(r['interval_net_melt_equivalent_mmwe'][0],y)
  for j,d,e in zip(J,H,r['interval_net_energy_J_m2']):self.assertAlmostEqual(k*j-d,e)
 def test_cooling_storage_releases_energy(self):self.assertAlmostEqual(predict(10,[1,1],[0,-L])['interval_net_melt_equivalent_mmwe'][1],11)
 def test_negative_energy_flag_without_clipping(self):
  r=predict(10,[1,-1],[0,0]);self.assertTrue(r['negative_interval_phase_conversion_unsupported']);self.assertLess(r['interval_net_melt_equivalent_mmwe'][1],0)
 def test_bad_inputs(self):
  with self.assertRaises(ValueError):slab_sensible_J_m2(1,2,.1,-.1,2e6)
  with self.assertRaises(ValueError):gradient_integral([0,0],[1,2],[0,0],.1)
if __name__=='__main__':unittest.main()
