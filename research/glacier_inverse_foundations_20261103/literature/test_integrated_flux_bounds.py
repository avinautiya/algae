import itertools
import unittest
from integrated_flux_bounds import residual_radius,integrated_residual,basal_energy_interval

class TestIntegratedBounds(unittest.TestCase):
    def test_all_shared_endpoint_error_corners(self):
        # Same endpoints simultaneously enter M*dT and trapezoidal integral.
        mass=[3,7,2];delta=[.1,.4,.2];edges=[(2,1,0),(5,1,2)];duration=9
        radius=residual_radius(mass,delta,edges,duration)
        for signs in itertools.product((-1,1), repeat=6):
            a=[d*s for d,s in zip(delta,signs[:3])]
            b=[d*s for d,s in zip(delta,signs[3:])]
            ints=[duration*(u+v)/2 for u,v in zip(a,b)]
            err=integrated_residual(mass,a,b,ints,edges)
            self.assertLessEqual(abs(err),radius+1e-12)

    def test_arbitrary_time_nodes_obey_pointwise_bound(self):
        # 3 sensors x 3 time points ->512 coherent trajectories.
        mass=[1,4,1];delta=[.1,.2,.3];edges=[(2,1,0),(3,1,2)];duration=2
        radius=residual_radius(mass,delta,edges,duration)
        for signs in itertools.product((-1,1), repeat=9):
            t=[[delta[j]*signs[3*i+j] for j in range(3)] for i in range(3)]
            ints=[(t[0][j]+2*t[1][j]+t[2][j])/2 for j in range(3)]
            err=integrated_residual(mass,t[0],t[-1],ints,edges)
            self.assertLessEqual(abs(err),radius+1e-12)

    def test_common_capacity_conductance_scaling(self):
        m=[1,4,1];a=[1,2,1];b=[2,3,1];s=[3,8,3];g=[(2,1,0),(3,1,2)]
        r=integrated_residual(m,a,b,s,g)
        scale=7
        rs=integrated_residual([scale*x for x in m],a,b,s,[(scale*x,i,j) for x,i,j in g])
        self.assertAlmostEqual(rs,scale*r)

    def test_signed_flux_bound_not_melt_clipped(self):
        self.assertEqual(basal_energy_interval(1,4,0,2,.5),(-4,4))
        self.assertEqual(basal_energy_interval(1,4,3,2,.5),(2,16))

    def test_invalid_nonpassive_input(self):
        with self.assertRaises(ValueError):
            residual_radius([1],[.1],[(-1,0,0)],2)
        with self.assertRaises(ValueError):
            basal_energy_interval(2,1,1,1,.1)

if __name__=='__main__': unittest.main()
