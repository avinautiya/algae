import unittest
import numpy as np
from adaptive_flux import phase_deliver
class PhaseLedgerTests(unittest.TestCase):
 def test_positive_heat_latent_limit(self):
  q=np.array([1.,3.,7.]);m,c=phase_deliver(q);np.testing.assert_array_equal(m,q);self.assertEqual(c,0.)
 def test_cold_deficit_delays_melt(self):
  m,c=phase_deliver([-10.,3.,5.,4.]);np.testing.assert_array_equal(m,[0.,0.,0.,2.]);self.assertEqual(c,0.)
 def test_arbitrary_signed_sequence_energy_identity(self):
  for q in [[1.,-3.,1.],[-1.,-2.,5.,-7.],[0.],[]]:
   m,c=phase_deliver(q);self.assertAlmostEqual(sum(q),sum(m)+c);self.assertTrue((m>=0).all());self.assertLessEqual(c,0.)
 def test_prefix_causality(self):
  a=phase_deliver([-10.,3.,5.,4.,8.])[0];b=phase_deliver([-10.,3.,5.,100.,-100.])[0];np.testing.assert_array_equal(a[:3],b[:3])
 def test_invalid_input(self):
  for q in [[float('nan')],[float('inf')],[[1.,2.]]]:
   with self.assertRaises(ValueError):phase_deliver(q)
if __name__=='__main__':unittest.main()
