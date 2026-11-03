import unittest
import numpy as np
from scipy.linalg import expm
from adjoint_bound import right_absorption_weights,initial_heat_bound
class AdjointTests(unittest.TestCase):
 def test_uniform_chain_linear_coordinate(self):
  np.testing.assert_allclose(right_absorption_weights([2]*4),[.25,.5,.75])
 def test_independent_matrix_exponential_flux_identity(self):
  g=np.array([1.,2.,3.,4.]);c=np.array([3.,2.,5.]);p=right_absorption_weights(g)
  K=np.diag(g[:-1]+g[1:])-np.diag(g[1:-1],1)-np.diag(g[1:-1],-1);M=-K/c[:,None];u0=np.array([2.,-1.,4.])
  # Augmented state integrates physical right heat independently of adjoint.
  B=np.zeros((4,4));B[:3,:3]=M;B[3,2]=g[-1]
  for t in [.1,1.,10.]:
   state=expm(B*t)@np.r_[u0,0.]
   self.assertAlmostEqual(state[3],np.dot(c*p,u0-state[:3]),places=11)
   self.assertLessEqual(abs(state[3]),initial_heat_bound(g,c,abs(u0))['initial_state_right_heat_error_bound_J_m2']+1e-10)
 def test_invalid_passivity_refused(self):
  with self.assertRaises(ValueError):right_absorption_weights([1.,0.])
if __name__=='__main__':unittest.main()
