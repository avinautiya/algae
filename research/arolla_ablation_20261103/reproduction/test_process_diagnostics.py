import datetime as dt, math, unittest
from process_diagnostics import SIGMA,equivalent_flux,radiometric_temperature,process_window
class Tests(unittest.TestCase):
 def test_latent_dimensional_identity(self):
  self.assertAlmostEqual(equivalent_flux(.1,30*86400)*30*86400,100*334000)
  self.assertLess(equivalent_flux(-.1,86400),0)
  for v in [0,-1,float('nan')]:
   with self.assertRaises(ValueError):equivalent_flux(.1,v)
 def test_radiometry_reflected_downwelling(self):
  kelvin=283.15;eps=.93;incoming=280
  outgoing=eps*SIGMA*kelvin**4+(1-eps)*incoming
  self.assertAlmostEqual(radiometric_temperature(outgoing,incoming,eps),10.,places=10)
 def test_incomplete_support_refuses(self):
  with self.assertRaises(ValueError):process_window([],dt.datetime(2015,1,1),dt.datetime(2015,1,2))
 def test_invalid_sentinel_refuses(self):
  t=dt.datetime(2015,1,1)
  row=dict(time=t,T=1,FF=1,RH=50,SWIN=10,SWOUT=2,LWIN=-9999,LWOUT=300,PP=0,SNOW=0)
  with self.assertRaises(ValueError):process_window([row],t,t+dt.timedelta(hours=1))
if __name__=='__main__':unittest.main()
