import math,unittest
from radiation_reservoir import step,radiation_proxy,SIGMA
class Checks(unittest.TestCase):
 def test_steady_warm_state(self):
  e,g=step(100*3600,100,7200,3600);self.assertAlmostEqual(e,360000);self.assertAlmostEqual(g,720000)
 def test_cooling_then_rewarming_repays_deficit(self):
  e,g=step(0,-100,3600,3600);self.assertEqual(e,-360000);self.assertEqual(g,0)
  e,g=step(e,100,3600,3600);self.assertEqual(e,0);self.assertEqual(g,0)
  e,g=step(e,100,3600,3600);self.assertGreater(e,0);self.assertGreater(g,0)
 def test_positive_to_negative_event(self):
  tau=10.;e0=20.;q=-2.;hit=tau*math.log(2)
  e,g=step(e0,q,20,tau);self.assertAlmostEqual(e,q*(20-hit));self.assertAlmostEqual(g,e0+q*hit)
 def test_composition_constant_forcing(self):
  for tau in [0.,1.,1000.,86400.]:
   for e0 in [-1e4,0,1e4]:
    for q in [-20,0,20]:
     e,g=step(e0,q,3600,tau);e1,g1=step(e0,q,1000,tau);e2,g2=step(e1,q,2600,tau)
     self.assertAlmostEqual(e,e2,7);self.assertAlmostEqual(g,g1+g2,7)
 def test_radiative_equilibrium(self):
  T=280.;self.assertAlmostEqual(radiation_proxy(0,SIGMA*T**4,T-273.15),0,12)
 def test_zero_tau_cold_store(self):
  e,g=step(-100,10,20,0);self.assertEqual(e,0);self.assertEqual(g,100)
 def test_invalid(self):
  with self.assertRaises(ValueError):step(0,0,0,1)
  with self.assertRaises(ValueError):radiation_proxy(-1,300,2)
if __name__=='__main__':unittest.main()
