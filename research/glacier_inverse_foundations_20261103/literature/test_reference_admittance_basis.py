import math,unittest
from reference_admittance_basis import ice_impedance,length_response,recover_two_real_laplace,warm_gauge
class ReferenceBasis(unittest.TestCase):
 def test_warm_gauge_is_exact(self):
  r=warm_gauge(.96,1.5e6,.03,4.)
  self.assertEqual(r['alpha'],r['scaled_alpha']);self.assertEqual(r['contact_length'],r['scaled_contact_length'])
 def test_known_ice_reference_breaks_warm_scale_symmetry(self):
  s=1/86400.;base=length_response(s,.96,.03,2.1,917*2100);scaled=length_response(s,3.84,.03/4,2.1,917*2100)
  self.assertGreater(abs(base-scaled),.01)
 def test_noise_free_two_frequency_recovery(self):
  s1,s2=1/86400.,1/21600.;k,C,R=.96,1.5e6,.03;ki,Ci=2.1,917*2100
  r=recover_two_real_laplace(s1,s2,length_response(s1,k,R,ki,Ci).real,length_response(s2,k,R,ki,Ci).real,k/C,ki,Ci)
  self.assertAlmostEqual(r['conductivity'],k);self.assertAlmostEqual(r['capacity'],C);self.assertAlmostEqual(r['contact_R'],R)
 def test_finite_domain_limits(self):
  s=1/86400.;ki,Ci=2.1,917*2100
  self.assertLess(abs(ice_impedance(s,ki,Ci,20)-ice_impedance(s,ki,Ci)),1e-12)
  self.assertAlmostEqual(ice_impedance(1e-12,ki,Ci,2).real,2/ki,places=5)
 def test_scaling_ice_reference_too_restores_ambiguity(self):
  s=1/86400.;first=length_response(s,.96,.03,2.1,917*2100);second=length_response(s,3.84,.03/4,8.4,4*917*2100)
  self.assertAlmostEqual(abs(first-second),0,places=12)
if __name__=='__main__':unittest.main()
