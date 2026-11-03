import unittest
from transfer_certificate import certificate,passive_rc_doses
class Tests(unittest.TestCase):
 def test_exact_training_scalar_scale_cancels(self):
  a=certificate([{'id':'a','training_dose':2,'future_dose':3},{'id':'b','training_dose':20,'future_dose':30}],100,0)
  self.assertEqual(a['lower'],150);self.assertEqual(a['upper'],150);self.assertFalse(a['abstain'])
 def test_minimax_lower_bound_witness(self):
  a=certificate([{'id':'a','training_dose':1,'future_dose':1},{'id':'b','training_dose':1,'future_dose':3}],10,9)
  self.assertEqual(a['minimax_point'],20);self.assertEqual(a['conditional_minimax_absolute_error'],10);self.assertTrue(a['abstain'])
 def test_passive_memory_counterexample(self):
  rows=[]
  for tau in [.1,10.]:
   t,f=passive_rc_doses(tau);rows.append({'id':str(tau),'training_dose':t,'future_dose':f})
  a=certificate(rows,100,5)
  self.assertGreater(a['upper']-a['lower'],100)
 def test_bad_and_empty_refused(self):
  for rows in [[],[{'id':'a','training_dose':0,'future_dose':1}],[{'id':'a','training_dose':1,'future_dose':float('nan')}]]:
   with self.assertRaises(ValueError):certificate(rows,1,1)
if __name__=='__main__':unittest.main()
