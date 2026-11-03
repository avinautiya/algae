import cmath,math,unittest
from contact_screen import transfer,same_sensor_resistance
class ContactBasisTests(unittest.TestCase):
 def test_zero_frequency_series_resistance(self):
  self.assertEqual(transfer(0,2,1e6,.5,.1),complex(1/(.25+.1)))
  self.assertAlmostEqual(transfer(1e-12,2,1e6,.5,.1).real,1/(.25+.1),places=8)
 def test_zero_contact_recovers_standard_diffusion_transfer(self):
  w=2*math.pi/86400;k=2.;C=1e6;d=.5;lam=cmath.sqrt(1j*w*C/k)
  self.assertAlmostEqual(abs(transfer(w,k,C,d,0)-k*lam/cmath.sinh(lam*d)),0.,places=14)
 def test_clearance_and_contact_steady_nonidentifiability(self):
  self.assertAlmostEqual(same_sensor_resistance(.01,.02,1),same_sensor_resistance(.02,.01,1))
 def test_contact_is_not_a_frequency_independent_gain(self):
  ratios=[transfer(2*math.pi/(h*3600),2,1e6,.5,.03)/transfer(2*math.pi/(h*3600),2,1e6,.5,0) for h in [6,24]]
  self.assertGreater(abs(ratios[0]-ratios[1]),.01)
 def test_invalid_parameters_refused(self):
  for args in [(0,2,1e6,.5,-1),(0,float('nan'),1e6,.5,0),(0,2,0,.5,0)]:
   with self.assertRaises(ValueError):transfer(*args)
if __name__=='__main__':unittest.main()
