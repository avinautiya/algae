import unittest
from .spectroscopy import OMEGA,transfer
from .profile_fit import fit_geometry,rank_geometry
from .test_spectroscopy import fv_harmonic
class Tests(unittest.TestCase):
 def test_nondefault_geometry_fit_null_rank(self):
  L,x=.23,.185;alpha,b=9e-7,.027;obs=[]
  for w in [0,OMEGA,2*OMEGA]:
   z=fv_harmonic(alpha,b,w,272,L=L,x=x);obs.append({'omega_rad_s':w,'ratio_real':z.real,'ratio_imag':z.imag,'ratio_error_scale':.01})
  fit=fit_geometry(obs,L=L,x=x);null=fit_geometry(obs,L=L,x=x,null=True);wrong=fit_geometry(obs,L=.17,x=.1)
  self.assertLess(abs(fit['alpha_m2_s']/alpha-1),.003);self.assertLess(abs(fit['b_kR_m']-b),.0003);self.assertGreater(null['weighted_objective'],fit['weighted_objective']);self.assertGreater(wrong['weighted_objective'],fit['weighted_objective']);self.assertGreater(rank_geometry(fit['alpha_m2_s'],fit['b_kR_m'],obs,L=L,x=x)['singular_values'][1],1.)
 def test_probe_below_ice_refused(self):
  with self.assertRaises(ValueError):fit_geometry([],L=.20,x=.225)
if __name__=='__main__':unittest.main()
