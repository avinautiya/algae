"""Independent semidiscrete conduction oracle, no field targets."""
import sys
import unittest
from pathlib import Path
import numpy as np
sys.path.insert(0,str(Path(__file__).parents[1]/'theory'))
from weak_graph import fit,lowest_edge_enclosure

class IndependentWeakGraphTests(unittest.TestCase):
    @staticmethod
    def fixture():
        z=np.array([0,.07,.19,.31,.5]);C=1.8e6;G=np.array([17.,4.,9.,11.])
        n=len(z);M=np.zeros((n,n));K=np.zeros((n,n))
        for j,h in enumerate(np.diff(z)):
            M[j:j+2,j:j+2]+=C*h/6*np.array([[2,1],[1,2]])
            K[j:j+2,j:j+2]+=G[j]*np.array([[1,-1],[-1,1]])
        dt=900.;t=np.arange(0,5*86400+dt,dt);T=np.zeros((len(t),n))
        T[0]=[9,7,4,2,0]
        for a in range(len(t)-1):
            Tb=np.array([9+3*np.sin(t[a+1]*2*np.pi/86400),0.])
            A=M/dt+K/2;rhs=(M/dt-K/2)@T[a]
            T[a+1,[0,-1]]=Tb
            T[a+1,1:-1]=np.linalg.solve(A[1:-1,1:-1],rhs[1:-1]-A[1:-1][:,[0,-1]]@Tb)
        return z,t,T,C,G

    def test_recovers_transient_heterogeneous_conductances(self):
        z,t,T,C,G=self.fixture()
        result=fit(z,t,T,capacity_J_m3K=C,window_s=21600.)
        np.testing.assert_allclose(result['conductances_W_m2K'],G,rtol=2e-10,atol=2e-10)
        self.assertEqual(result['rank'],len(G))

    def test_conditional_lp_contains_true_deepest_edge_after_sensor_errors(self):
        z,t,T,C,G=self.fixture()
        # Persistent offsets plus bounded time noise, shared across all windows.
        offsets=np.array([.1,-.2,.15,-.3,.05])
        perturbed=T+offsets+.05*np.sin(np.arange(len(t))[:,None]*np.array([.7,.5,.3,.9,.4]))
        result=lowest_edge_enclosure(z,t,perturbed,capacity_J_m3K=C,
                                     window_s=21600.,sensor_bound_C=.4)
        lo,hi=result['lowest_conductance_W_m2K']
        self.assertLessEqual(lo,G[-1]);self.assertGreaterEqual(hi,G[-1])
        self.assertFalse(result['icemelt_prediction'])

if __name__=='__main__':unittest.main()
