import math,unittest
from thermal_path import Path,evaluate,example,CW
class Tests(unittest.TestCase):
 def test_global_energy(self):
  for q in [.0001,.01,1.]:
   for u in [0.,10.,10000.]:
    r=evaluate(Path(q,u,200.,1.,.5));self.assertLess(abs(r['residual_W']),1e-12*max(1,r['inlet_W']))
 def test_serial_subdivision(self):
  a=evaluate(Path(.01,100.,200.,1.,.5));b=evaluate(Path(.01,0.,300.,1.,.5));self.assertAlmostEqual(a['upstream_melt_W']+a['target_melt_W'],b['target_melt_W'],places=12)
 def test_zero_flow(self):self.assertEqual(evaluate(Path(0.,100.,200.,1.,.5))['target_melt_W'],0.)
 def test_ranking_reversal_is_local_not_energy_creation(self):
  r=example();self.assertGreater(r['no_upstream_paths']['gap_drift_m_day'],0.);self.assertLess(r['upstream_bypass_path']['gap_drift_m_day'],0.)
 def test_exact_critical_boundary(self):
  e=example();p=Path(**e['bypass_direct']);c=e['critical_upstream_bypass_conductance_W_K'];r=evaluate(Path(p.flow_m3_s,c,p.target_conductance_W_K,p.target_area_m2,p.inlet_superheat_K));self.assertAlmostEqual(r['incision_m_s'],e['no_upstream_paths']['toe']['incision_m_s'],places=18)
if __name__=='__main__':unittest.main()
