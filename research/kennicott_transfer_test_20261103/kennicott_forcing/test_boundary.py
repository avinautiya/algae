import unittest,datetime as dt
from .boundary import _at,ols,_seconds,_label,build_boundary
class Tests(unittest.TestCase):
 def test_linear_interpolation(self):self.assertEqual(_at([0,2],[1,5],1,2),3)
 def test_refuse_extrapolation(self):
  with self.assertRaises(ValueError):_at([0,2],[1,5],3,2)
 def test_refuse_long_gap(self):
  with self.assertRaises(ValueError):_at([0,3],[1,5],1,2)
 def test_regression(self):self.assertEqual(ols([0,1,2],[2,5,8]),(2,3))
 def test_identical_air_refused(self):
  with self.assertRaises(ValueError):ols([1,1],[1,2])
 def test_utc_round_trip(self):self.assertEqual(_label(_seconds('2023-06-08T20:00:00Z')),'2023-06-08T20:00:00Z')
 def test_nonutc_refused(self):
  with self.assertRaises(ValueError):_seconds('2023-06-08T20:00:00-08:00')
 def test_source_boundary_and_shift(self):
  a=build_boundary('KDT_2020b','2020-07-25T20:00:00Z','2020-09-06T20:00:00Z')
  b=build_boundary('KDT_2020b','2020-07-25T20:00:00Z','2020-09-06T20:00:00Z',proxy_error_multiplier=1)
  self.assertEqual(len(a['boundary_C']),1033);self.assertEqual(a['metadata']['proxy_hours'],156)
  self.assertEqual(a['metadata']['fit_rows'],337);self.assertEqual(a['metadata']['blocked_fit_rows'],168)
  for x,y,flag in zip(a['boundary_C'],b['boundary_C'],a['is_proxy']):
   self.assertAlmostEqual(y-x,2*a['metadata']['blocked_RMSE_C'] if flag else 0.,places=12)
  self.assertEqual(len(a['proxy_boundary_C']),len(a['air_C']))
 def test_unknown_source_refused(self):
  with self.assertRaises(ValueError):build_boundary('unknown','2023-06-08','2023-06-09')
if __name__=='__main__':unittest.main()
