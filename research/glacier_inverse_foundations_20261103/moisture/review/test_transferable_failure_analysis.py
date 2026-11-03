import math,unittest
from transferable_failure_analysis import distribution,fosd,expfactor,increasing_convex_order,icx_moisture_error_radius
class FOSDTests(unittest.TestCase):
 def test_dominating_distribution_all_beta(self):
  a=distribution([0,.1,.2],[1,2,1]);b=distribution([.1,.2,.3],[1,2,1])
  self.assertTrue(fosd(b,a)['test_first_order_stochastically_dominates_train'])
  for i in range(401):
   beta=i/10;factor=expfactor(b,a,beta)
   self.assertAlmostEqual(factor,math.exp(beta*.1),places=10)
   self.assertGreaterEqual(factor,1)
   # PDD=12, observed10; any factor>=1 cannot improve error.
   self.assertGreaterEqual(abs(12*factor-10),2)
 def test_mean_increase_is_not_fosd(self):
  a=distribution([.1,.2],[1,1]);b=distribution([0,.4],[1,1])
  self.assertFalse(fosd(b,a)['test_first_order_stochastically_dominates_train'])
  self.assertGreater(sum(x*w for x,w in b),sum(x*w for x,w in a))
 def test_icx_can_pass_when_fosd_fails(self):
  a=distribution([.1,.2],[1,1]);b=distribution([0,.4],[1,1])
  self.assertFalse(fosd(b,a)['test_first_order_stochastically_dominates_train'])
  self.assertTrue(increasing_convex_order(b,a)['test_increasing_convex_dominates_train'])
  for i in range(401):self.assertGreaterEqual(expfactor(b,a,i/10),1.)
 def test_higher_mean_does_not_suffice_for_icx(self):
  a=distribution([0,1],[1,1]);b=distribution([.7],[1])
  self.assertFalse(increasing_convex_order(b,a)['test_increasing_convex_dominates_train'])
  self.assertLess(expfactor(b,a,40),1.)
 def test_error_radius_exact_translation(self):
  a=distribution([0,.1,.2],[1,2,1]);b=distribution([.1,.2,.3],[1,2,1])
  # Mean and every stop-loss ordering fail above .1/2 symmetric shift.
  r=icx_moisture_error_radius(b,a)['maximum_symmetric_pointwise_theta_error_radius']
  self.assertAlmostEqual(r,.05,places=10)
 def test_zero_weights_and_equal_distribution(self):
  a=distribution([0,2,10],[1,1,0]);b=distribution([2,0],[1,1]);self.assertTrue(fosd(a,b)['test_first_order_stochastically_dominates_train'])
  for beta in [0,.1,1,10]:self.assertAlmostEqual(expfactor(a,b,beta),1.)
if __name__=='__main__':unittest.main()
