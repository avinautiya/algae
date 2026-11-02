import unittest
import numpy as np
from run import node,integral,read,stamp
class Checks(unittest.TestCase):
 def test_utc(self):self.assertEqual(stamp('1970-01-01 00:00:00'),0)
 def test_node_equilibrium(self):
  t=np.arange(1000)*300.;np.testing.assert_allclose(node(t,np.ones(1000)*8,.2,1e-6,4),4,atol=1e-12)
 def test_node_relaxation(self):
  t=np.arange(100)*300.;a=1e-6;L=.2;y=node(t,np.ones(100)*8,L,a,0)
  np.testing.assert_allclose(y,4*(1-np.exp(-4*a*t/L**2)),atol=1e-12)
 def test_integral_phase_split(self):
  t=np.array([0,1,2.]);p,n,a,b=integral(t,np.array([2.,2,2]),.25,1.25)
  self.assertAlmostEqual(p,2);self.assertEqual(n,0)
 def test_no_synthetic_probe(self):
  t,z,T=read('2020b');self.assertEqual(len(z),2);np.testing.assert_allclose(z,[.05,.15])
if __name__=='__main__':unittest.main()
