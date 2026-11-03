import unittest
from decision_rule import evaluate
class Decision(unittest.TestCase):
 def m(self,es):return {'absolute_errors_mmwe':es,'MAE_mmwe':sum(es)/3}
 def test_requires_all_baselines(self):
  m={'physical_layered_local':self.m([1,1,1]),'homogeneous_empirical_scale':self.m([2,2,2]),'other':self.m([.9,.9,.9])};self.assertFalse(evaluate(m)['flag'])
 def test_two_stakes_and_mae(self):
  m={'physical_layered_local':self.m([1,1,3]),'homogeneous_empirical_scale':self.m([2,2,4])};self.assertTrue(evaluate(m)['flag'])
 def test_zero_tie_not_gain(self):
  m={'physical_layered_local':self.m([0,0,0]),'homogeneous_empirical_scale':self.m([0,0,0])};self.assertFalse(evaluate(m)['flag'])
 def test_missing_comparator_not_win(self):self.assertFalse(evaluate({'physical_layered_local':self.m([1,1,1])})['flag'])
if __name__=='__main__':unittest.main()
