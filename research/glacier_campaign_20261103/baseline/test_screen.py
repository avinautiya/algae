import unittest
from screen import run,rates,partition,H,T,RHOI,L,HC,GAP,Q0
class Checks(unittest.TestCase):
 def test_twelve_cases_no_21day_reset(self):
  r=run();self.assertEqual(len(r['cases']),12);self.assertFalse(any(x['can_repair_in_21days_under_ideal_assumptions'] for x in r['cases']))
 def test_heat_limit_and_monotonicity(self):
  es=[rates(q,1.) for q in [0,.001,.01,.1,1.]];self.assertEqual(es,sorted(es));self.assertLess(es[-1],H*T/(RHOI*L))
 def test_partition(self):
  for q in [Q0,2*Q0,8*Q0]:
   a,b=partition(q,GAP);self.assertAlmostEqual(a+b,q);self.assertGreaterEqual(b,q/2)
 def test_bound_on_entire_repair_interval(self):
  for at in [1.,.5,.25]:
   for q in [Q0,2*Q0,4*Q0,8*Q0]:
    upper=H*T/(RHOI*L)-rates(q/2,1.)
    for j in range(101):
     g=HC+(GAP-HC)*j/100;qt,qb=partition(q,g)
     self.assertLessEqual(rates(qt,at)-rates(qb,1.),upper+1e-20)
if __name__=='__main__':unittest.main()
