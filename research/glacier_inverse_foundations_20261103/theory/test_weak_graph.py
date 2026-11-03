import unittest
import numpy as np
from weak_graph import fit,lowest_edge_enclosure
class WeakTests(unittest.TestCase):
    def test_scaled_temperature_problem_preserves_scale_symmetry(self):
        z=[0,.1,.2,.3];t=np.arange(0,86400*3+1,3600.)
        T=np.array([[6+2*np.sin(x/86400*6.283),4+np.sin(x/86400*6.283-.3),2+.4*np.sin(x/86400*6.283-.8),0]for x in t])
        a=fit(z,t,T,capacity_J_m3K=1e6,conductivity_bounds_W_mK=(.001,100))
        b=fit(z,t,T,capacity_J_m3K=2e6,conductivity_bounds_W_mK=(.002,200))
        np.testing.assert_allclose(b['conductances_W_m2K'],2*np.array(a['conductances_W_m2K']),rtol=1e-9)
    def test_steady_series_flux_ratio(self):
        z=[0,.1,.2,.3];t=np.arange(0,86400*3+1,3600.)
        T=np.tile([10,8,3,0],(len(t),1))
        r=fit(z,t,T,capacity_J_m3K=1e6)
        flux=np.array(r['conductances_W_m2K'])*np.array([2,5,3])
        np.testing.assert_allclose(flux,flux[0],rtol=1e-10)
        self.assertLess(r['balance_RMSE_W_m2'],1e-10)
    def test_bad_input(self):
        with self.assertRaises(ValueError):fit([0,.1,.2],[0,0,1],np.zeros((3,3)),capacity_J_m3K=1e6)
        for value in (float('nan'),float('inf')):
            with self.assertRaises(ValueError):fit([0,.1,.2],[0,1,2],np.zeros((3,3)),capacity_J_m3K=1e6,sensor_bound_C=value)
    def test_necessary_lp_contains_exact_steady_conductances(self):
        z=[0,.1,.2,.3];t=np.arange(0,86400*3+1,3600.);T=np.tile([10,8,3,0],(len(t),1))
        r=lowest_edge_enclosure(z,t,T,capacity_J_m3K=1e6,sensor_bound_C=0.)
        # Exact steady flux30W/m² gives G=[15,6,10]W/m²/K.
        lo,hi=r['lowest_conductance_W_m2K'];self.assertLessEqual(lo,10);self.assertLessEqual(10,hi)
        self.assertFalse(r['exact_joint_sensor_error_feasibility_established'])
if __name__=='__main__':unittest.main()
