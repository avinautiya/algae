import unittest
import numpy as np
from ratio_enclosure import enclosure
class RatioTests(unittest.TestCase):
 def test_independent_vertex_enumeration(self):
  # x,y in[1,4], x+y<=5: all three polygon vertices enumerated independently.
  vs=[(1,1),(1,4),(4,1)]
  expected=[100*(2*x+y)/(x+3*y)for x,y in vs]
  r=enclosure([[1,1]],[5],[(1,4),(1,4)],[1,3],[2,1],100)
  self.assertAlmostEqual(r['lower_mmwe'],min(expected));self.assertAlmostEqual(r['upper_mmwe'],max(expected))
 def test_proportional_goal_identified_without_parameter_uniqueness(self):
  r=enclosure([],[],[(1,4),(2,9)],[1,2],[3,6],100,1e-8)
  self.assertAlmostEqual(r['lower_mmwe'],300);self.assertAlmostEqual(r['upper_mmwe'],300);self.assertFalse(r['abstain'])
 def test_common_parameter_scaling_cancels(self):
  a=enclosure([],[],[(1,4),(2,9)],[1,2],[3,1],100)
  b=enclosure([],[],[(7,28),(14,63)],[1,2],[3,1],100)
  self.assertAlmostEqual(a['lower_mmwe'],b['lower_mmwe']);self.assertAlmostEqual(a['upper_mmwe'],b['upper_mmwe'])
 def test_infeasible_not_certain(self):
  self.assertEqual(enclosure([[1]],[0],[(1,2)],[1],[1],100)['status'],'infeasible_assumptions')
 def test_zero_denominator_and_negative_heat_refused(self):
  with self.assertRaises(ValueError):enclosure([],[],[(0,1)],[1],[1],100)
  with self.assertRaises(ValueError):enclosure([],[],[(1,2)],[1],[-1],100)
if __name__=='__main__':unittest.main()

class CoefficientBoxTests(unittest.TestCase):
 def test_independent_corner_enumeration(self):
  from ratio_enclosure import coefficient_box_enclosure
  from itertools import product
  cl,cu=[1,2],[2,3];fl,fu=[2,1],[3,4]
  r=coefficient_box_enclosure([],[],[(1,4),(2,5)],cl,cu,fl,fu,100)
  samples=[100*(ff[0]*x+ff[1]*y)/(cc[0]*x+cc[1]*y)for x,y in product([1,4],[2,5])for cc in product(*zip(cl,cu))for ff in product(*zip(fl,fu))]
  self.assertAlmostEqual(r['lower_mmwe'],min(samples));self.assertAlmostEqual(r['upper_mmwe'],max(samples))
 def test_negative_parameter_refused(self):
  from ratio_enclosure import coefficient_box_enclosure
  with self.assertRaises(ValueError):coefficient_box_enclosure([],[],[(-1,2)],[1],[2],[1],[2],100)
