import itertools,math,random,unittest
from audit import worst_error,dominance,fixed_forecast

def independent_corners(x,y,fa,fb,d,association):
 out=[]
 if association=='per_reading':
  for e0,e1,e2 in itertools.product((-d,d),repeat=3):
   truex=x+e1-e0;truey=y+e2-e1
   out.append((fa*truex-truey,fb*truex-truey))
 else:
  for etrain,etest in itertools.product((-d,d),repeat=2):
   out.append((fa*(x+etrain)-(y+etest),fb*(x+etrain)-(y+etest)))
 return out
class TestBounds(unittest.TestCase):
 def test_analytic_matches_independent_corners(self):
  rng=random.Random(431)
  for association in ['per_reading','per_interval']:
   for _ in range(500):
    d=rng.uniform(.01,100);x=rng.uniform(2.1*d,10*d);y=rng.uniform(2.1*d,10*d);fa,fb=rng.uniform(0,3),rng.uniform(0,3)
    corners=independent_corners(x,y,fa,fb,d,association)
    got=worst_error(x,y,fa,d,association)
    self.assertAlmostEqual(got['worst_absolute_error_mmwe'],max(abs(a) for a,b in corners),places=9)
    self.assertAlmostEqual(got['signed_error_bounds_mmwe'][0],min(a for a,b in corners),places=9)
    self.assertAlmostEqual(got['signed_error_bounds_mmwe'][1],max(a for a,b in corners),places=9)
    signs=[a*a-b*b for a,b in corners]
    status=dominance(x,y,fa,fb,d,association)['status']
    expected='a_strictly_better_everywhere' if max(signs)<0 else 'b_strictly_better_everywhere' if min(signs)>0 else 'ranking_changes_with_admissible_measurement_error'
    self.assertEqual(status,expected)
 def test_fixed_forecast_does_not_recalibrate_under_measurement_error(self):
  for association in ['per_reading','per_interval']:
   r=180 if association=='per_reading' else 90
   result=fixed_forecast(846,315,.4,90,association)
   errors=[abs(.4*846-(315+e)) for e in [-r,r]]
   self.assertAlmostEqual(result['worst_absolute_error_mmwe'],max(errors))
   self.assertLess(result['worst_absolute_error_mmwe'],worst_error(846,315,.4,90,association)['worst_absolute_error_mmwe'])
 def test_shared_middle_reading_sign(self):
  # This vertex simultaneously increases calibration and decreases heldout truth.
  c=independent_corners(846,315,.4,.5,90,'per_reading')
  expected=.4*(846+180)-(315-180)
  self.assertAlmostEqual(max(a for a,b in c),expected)
 def test_equal_models_and_boundary_tie(self):
  self.assertEqual(dominance(10,10,1,1,1,'per_interval')['status'],'identical_predictions')
  self.assertEqual(dominance(10,10,0,1,0,'per_interval')['status'],'b_strictly_better_everywhere')
 def test_unknown_or_nonphysical_uncertainty_refuses(self):
  for args in [(100,100,1,90,'per_reading'),(100,100,1,1,'unknown'),(100,100,float('nan'),1,'per_reading')]:
   with self.assertRaises(ValueError):worst_error(*args)
if __name__=='__main__':unittest.main()
