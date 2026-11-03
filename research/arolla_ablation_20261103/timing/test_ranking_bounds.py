import itertools,unittest
from ranking_bounds import rank_bound
class Tests(unittest.TestCase):
 def test_nominal_win_not_robust_win(self):
  r=rank_bound(846,315,.4,.45,90,association='each_shared_reading');self.assertLess(r['nominal_error_a'],r['nominal_error_b']);self.assertEqual(r['status'],'ranking_not_identified_under_bounds')
 def test_exact_reading_error_radius(self):
  train,test,fa,fb,d=100.,80.,.5,.6,2.;r=rank_bound(train,test,fa,fb,d,association='each_shared_reading');vals=[]
  for e0,e1,e2 in itertools.product([-d,d],repeat=3):vals.append((fa+fb)*(train+e1-e0)-2*(test+e2-e1))
  self.assertAlmostEqual(min(vals),r['midpoint_residual_bounds'][0]);self.assertAlmostEqual(max(vals),r['midpoint_residual_bounds'][1])
 def test_zero_uncertainty_matches_absolute_errors(self):
  r=rank_bound(100,80,.5,.6,0,association='each_interval');self.assertEqual(r['status'],'b_has_lower_absolute_error_throughout_bounds')
 def test_identical_is_not_improvement(self):
  self.assertEqual(rank_bound(100,80,.5,.5,1,association='each_interval')['status'],'identical_predictions')
 def test_unknown_and_nonpositive_training_refuse(self):
  for assoc in ['unknown','each_shared_reading']:
   with self.assertRaises(ValueError):rank_bound(1,1,.5,.6,1,association=assoc)
if __name__=='__main__':unittest.main()
