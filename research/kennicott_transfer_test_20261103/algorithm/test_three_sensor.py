import math,unittest
from three_sensor import *

def oracle(n=241):
 # Independent consistent-mass, Crank-Nicolson three-node thermal column.
 z=[.05,.15,.35];hl=z[1]-z[0];hr=z[2]-z[1];alpha=[.8e-6,1.7e-6];g1=alpha[0]/hl;g2=alpha[1]/hr
 m1,m2,m3=hl/6,(hl+hr)/3,hr/6;step=3600.;t=[j*step for j in range(n)]
 def left(x):return 8+4*math.sin(x*2*math.pi/86400)
 def right(x):return .5+.2*math.sin(x*2*math.pi/(86400*1.7)+.4)
 T=[[left(0),4,right(0)]]
 for j in range(n-1):
  a,b=t[j],t[j+1];prev=T[-1];u,v=left(b),right(b)
  interior=((m2/step-(g1+g2)/2)*prev[1]-m1*(u-prev[0])/step-m3*(v-prev[2])/step+g1*(u+prev[0])/2+g2*(v+prev[2])/2)/(m2/step+(g1+g2)/2)
  T.append([u,interior,v])
 return z,t,T,alpha

class ThreeSensorTests(unittest.TestCase):
 def test_exact_independent_CN_recovery(self):
  z,t,T,truth=oracle();r=fit(z,t,T)
  self.assertTrue(r['supported_point_estimate'])
  for a,b in zip(r['alpha_m2_s'],truth):self.assertAlmostEqual(a/b,1,places=10)
 def test_steady_rank_degeneracy_refuses_point_flux(self):
  z=[0,.1,.3];t=[j*3600 for j in range(73)];T=[[10,7,1]for _ in t];r=fit(z,t,T)
  self.assertFalse(r['practically_identifiable_under_declared_threshold'])
  with self.assertRaises(UnsupportedModel):reconstruct(z,t,T,ice_depth_m=.4,basal_temperature_C=0.)
 def test_flat_zero_rank(self):
  r=fit([0,.1,.3],[j*3600 for j in range(25)],[[0,0,0]for _ in range(25)])
  self.assertEqual(r['diagnostics']['singular_ratio'],0);self.assertFalse(r['supported_point_estimate'])
 def test_signed_control_volume_ledger_and_aggregation(self):
  z,t,T,_=oracle();r=reconstruct(z,t,T,ice_depth_m=.42,basal_temperature_C=0.)
  trace=r['normalized_heat_trace']
  for row in trace:self.assertAlmostEqual(row['normalized_top_heat_K_m']-row['normalized_debris_storage_change_K_m'],row['normalized_interface_heat_K_m'])
  self.assertAlmostEqual(aggregate(trace,t[0],t[-1]),r['normalized_total_interface_heat_K_m'])
  with self.assertRaises(UnsupportedModel):aggregate(trace,t[0]+1,t[-1])
 def test_gap_storage_matches_independent_polygon_trapezoids(self):
  middle,deep,h,g,tb=3,1,.2,.07,-.3
  independent=(h/2)*(((middle+deep)/2)+deep)/2+g*(deep+tb)/2
  self.assertAlmostEqual(slab_moment_K_m(middle,deep,h,g,tb),independent)
 def test_capacity_scale_quotient(self):
  z,t,T,_=oracle();r=reconstruct(z,t,T,ice_depth_m=.42,basal_temperature_C=0.);trace=r['normalized_heat_trace'];h0=aggregate(trace,t[0],t[120]);h1=aggregate(trace,t[120],t[-1])
  for c in [1e6,2e6,3e6]:self.assertAlmostEqual(scale_quotient(100,c*h0,c*h1),scale_quotient(100,h0,h1))
 def test_phase_not_inferred_from_temperature_fit(self):
  z,t,T,_=oracle();r=reconstruct(z,t,T,ice_depth_m=.42,basal_temperature_C=0.)
  with self.assertRaises(UnsupportedModel):absolute_melt_mmwe(r,1.35e6,temperate_ice_confirmed=False,all_other_interface_terms_zero_confirmed=True)
  cold=dict(r,basal_temperature_C=-1.)
  with self.assertRaises(UnsupportedModel):absolute_melt_mmwe(cold,1.35e6,temperate_ice_confirmed=True,all_other_interface_terms_zero_confirmed=True)
 def test_past_fit_cannot_see_future_features(self):
  z,t,T,_=oracle();asof=t[144];r=fit_past(z,t,T,asof_seconds=asof,lookback_seconds=72*3600,minimum_history_seconds=72*3600)
  changed=[row[:]for row in T]
  for j in range(145,len(T)):changed[j]=[999,-999,999]
  s=fit_past(z,t,changed,asof_seconds=asof,lookback_seconds=72*3600,minimum_history_seconds=72*3600)
  self.assertEqual(r['alpha_m2_s'],s['alpha_m2_s']);self.assertEqual(r['fit_last_source_seconds'],asof)
 def test_source_reversal_gap_and_unknown_gap_refusals(self):
  z,t,T,_=oracle()
  with self.assertRaises(ValueError):fit(z,[0,3600,0],T[:3])
  with self.assertRaises(UnsupportedModel):fit(z,[0,3600,10800],T[:3])
  with self.assertRaises(UnsupportedModel):reconstruct(z,t,T,ice_depth_m=.3,basal_temperature_C=0.)
  with self.assertRaises(ValueError):fit(z,t,T,sensor_absolute_bound_C=float('nan'))
 def test_nnls_face_solution(self):
  alpha,diag=nnls2([[1,0],[0,1]],[2,-3]);self.assertEqual(alpha,(2.,0.))
if __name__=='__main__':unittest.main()

class AdditionalCausalityTests(unittest.TestCase):
 def test_future_nan_does_not_change_past_fit(self):
  z,t,T,_=oracle();asof=t[144];a=fit_past(z,t,T,asof_seconds=asof,lookback_seconds=72*3600,minimum_history_seconds=72*3600)
  for j in range(145,len(T)):T[j]=[float('nan')]*3
  b=fit_past(z,t,T,asof_seconds=asof,lookback_seconds=72*3600,minimum_history_seconds=72*3600)
  self.assertEqual(a['alpha_m2_s'],b['alpha_m2_s'])
