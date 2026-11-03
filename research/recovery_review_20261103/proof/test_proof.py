import unittest
from check_proof import cap,H,HC,GAP,SECONDS,root,Q0,melt,small_area_repair,require_fixed_geometry
class Tests(unittest.TestCase):
 def test_declared_no_go_all_cases(self):
  for m in [1,2,4,8]:self.assertLess(cap(Q0*m)*SECONDS,GAP-HC)
 def test_larger_discharge_tightens_upper_bound(self):
  self.assertTrue(all(cap(Q0*a)>cap(Q0*b) for a,b in [(1,2),(2,4),(4,8)]))
 def test_vary_h_not_fixed_film(self):
  h=root(lambda x:cap(Q0,x)*SECONDS,GAP-HC,200,100000)
  self.assertGreater(h/H,10);self.assertLess(h/H,20);self.assertAlmostEqual(cap(Q0,h)*SECONDS,GAP-HC,places=12)
 def test_film_bound(self):
  for q in [.00001,.07,1.]:
   for area in [.001,.5,100.]:self.assertLessEqual(melt(q,area),H*.5/(917*334000)*(1+1e-12))
 def test_slight_repair_does_not_imply_recovery(self):
  r=small_area_repair();self.assertGreater(r['repair_m'],0.);self.assertLess(r['repair_m'],GAP-HC);self.assertFalse(r['reaches_critical_gap']);self.assertAlmostEqual(r['repair_m'],.00010343479994640026,places=12)
 def test_fixed_geometry_gate(self):
  a=dict(divider_width_m=.25,total_pool_width_m=.75,lane_width_m=.25,pool_length_m=.5,reach_length_m=.5,tail_length_m=.25)
  self.assertTrue(require_fixed_geometry([a,a.copy()]))
  with self.assertRaises(ValueError):require_fixed_geometry([a,dict(a,divider_width_m=.125,total_pool_width_m=.625)])
if __name__=='__main__':unittest.main()
