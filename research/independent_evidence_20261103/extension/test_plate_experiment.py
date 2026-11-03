import unittest
import numpy as np
from plate_experiment import coefficient,prepare,models_for_site,epoch,START,CUT,validate_hourly,lagged
from dirichlet_slab import PreparedSlab
class Tests(unittest.TestCase):
 def test_hourly_lag_contract(self):
  validate_hourly([0,3600,7200])
  for a in [[0,1800,3600],[0,3600,10800],[0,0,3600]]:
   with self.assertRaises(ValueError):validate_hourly(a)
 def test_ridge_features_never_use_interior(self):
  a=np.column_stack([np.arange(100.),np.zeros(100),np.arange(100.)+2]);b=a.copy();b[:,1]=9999
  np.testing.assert_equal(lagged(a),lagged(b))

 def test_exact_flux_scale(self):self.assertAlmostEqual(coefficient([1,2,3],[2,4,6])['raw_signal_coefficient_W_m_K'],2.)
 def test_signed_orientation_is_not_clipped(self):
  d=coefficient([1,2],[-3,-6]);self.assertEqual(d['training_inferred_orientation'],-1);self.assertEqual(d['effective_conductivity_abs_W_m_K'],3.)
 def test_no_excitation_or_nonfinite_refuses(self):
  for x,y in [([0,0],[1,1]),([1,float('nan')],[1,1]),([1],[1,2])]:
   with self.assertRaises(ValueError):coefficient(x,y)
 def test_analytic_steady_bottom_gradient(self):
  t=np.arange(0.,10001.,1000.);s=PreparedSlab(t,np.ones(len(t)),np.zeros(len(t)),[.05],[.01,.05,.1],[1.,5/9,0.],modes=128,max_gap_s=1000.)
  np.testing.assert_allclose(s.bottom_downward_gradient(1e-6).gradient_K_m,np.full(len(t),1/.09),atol=1e-12)
 def test_heldout_flux_never_changes_fits(self):
  times=np.arange(epoch(START)-2*86400,epoch(CUT)+3*86400,3600.);top=3+np.sin(np.arange(len(times))/24);bottom=np.zeros(len(times));mid=top*5/9
  s=PreparedSlab(times-times[0],top,bottom,[.05],[.01,.05,.1],[top[0],mid[0],0.],modes=64,max_gap_s=3600.)
  q=.5*s.bottom_downward_gradient(1e-6).gradient_K_m+np.cos(np.arange(len(times))/11)
  p=dict(times=times,temp=np.column_stack([top,mid,bottom]),slab=s,q=q.copy(),train=(times>=epoch(START))&(times<epoch(CUT)),gradient=top/.09)
  a,fa=models_for_site(p,1e-6);p['q'][times>=epoch(CUT)]+=999.;b,fb=models_for_site(p,1e-6)
  self.assertEqual(fa,fb)
  for name in a:np.testing.assert_equal(a[name],b[name])
if __name__=='__main__':unittest.main()
