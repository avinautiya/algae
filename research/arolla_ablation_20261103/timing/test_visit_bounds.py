import random,unittest
from visit_bounds import integrate,ratio,exact_bounds,residual_uncertainty_bound
class Tests(unittest.TestCase):
 def test_constant_calendar(self):
  s=[(0.,46.,1.)];r=exact_bounds(s,(0.,1.),(31.,32.),(45.,46.));self.assertAlmostEqual(r['minimum_ratio'],13/32);self.assertAlmostEqual(r['maximum_ratio'],.5)
 def test_shared_middle_extrema(self):
  rng=random.Random(3);s=[(i,i+1,rng.random()*4) for i in range(46)];r=exact_bounds(s,(0.,1.),(31.,32.),(45.,46.))
  for _ in range(100):
   x=ratio(s,rng.random(),31+rng.random(),45+rng.random());self.assertLessEqual(r['minimum_ratio']-1e-12,x);self.assertGreaterEqual(r['maximum_ratio']+1e-12,x)
 def test_gap_never_filled(self):
  with self.assertRaises(ValueError):integrate([(0,1,1),(2,3,1)],0,3)
 def test_zero_training_refuses(self):
  with self.assertRaises(ValueError):ratio([(0,2,0),(2,3,1)],0,1,3)
 def test_negative_driver_refuses(self):
  with self.assertRaises(ValueError):integrate([(0,1,-1)],0,1)
 def test_reading_vs_interval_error(self):
  self.assertEqual(residual_uncertainty_bound(.5,90,association='each_shared_reading'),270);self.assertEqual(residual_uncertainty_bound(.5,90,association='each_interval'),135)
if __name__=='__main__':unittest.main()
