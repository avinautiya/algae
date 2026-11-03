import unittest,math
import numpy as np
from moving_resistance import incremental_melt,RHO,LATENT
class MovingResistanceTests(unittest.TestCase):
 def test_known_train_and_uniform_heat_deceleration(self):
  r=incremental_melt([1000.]*5,80.,.5);m=r['period_melt_mmwe']
  self.assertAlmostEqual(m[0],80.,places=10)
  self.assertTrue(all(b<a for a,b in zip(m,m[1:])))
 def test_classical_integrated_stefan_energy(self):
  doses=np.array([1000.,2000.,0.,3000.]);r=incremental_melt(doses,80.,.5);H=np.array(r['depth_end_m'])
  lhs=.5*RHO*LATENT*(H[1:]**2-H[:-1]**2)
  np.testing.assert_allclose(lhs,r['inferred_k_W_mK']*doses,rtol=1e-12,atol=1e-8)
  self.assertEqual(r['period_melt_mmwe'][2],0.)
 def test_large_depth_recovers_fixed_resistance_ratio(self):
  r=incremental_melt([1000.,2000.,500.],80.,1e7)
  np.testing.assert_allclose(r['period_melt_mmwe'],[80.,160.,40.],rtol=1e-7)
 def test_invalid_inputs(self):
  for v in [([0.,1.],80.,.5),([1.,float('nan')],80.,.5),([1.,2.],80.,0.),([1.,2.],-1.,.5)]:
   with self.assertRaises(ValueError):incremental_melt(*v)
if __name__=='__main__':unittest.main()
