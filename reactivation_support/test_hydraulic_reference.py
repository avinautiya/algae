import unittest,math
from scipy.integrate import solve_ivp
from .hydraulic_reference import upstream_head,unit_discharge,partition,Lane,_length_function

class ReferenceTests(unittest.TestCase):
    def test_frictionless_critical_energy(self):
        for q in (.001,.07,.3):
            self.assertAlmostEqual(upstream_head(q,2,0),1.5*(q*q/9.81)**(1/3),places=12)
    def test_matches_independent_backward_gvf_ode(self):
        g=9.81
        for q,n,L in ((.03,.02,2.),(.1,.03,3.),(.2,.01,.5)):
            hc=(q*q/g)**(1/3)
            def rhs(x,h):return (n*n*q*q/h**(10/3))/(1-q*q/(g*h**3))
            ode=solve_ivp(rhs,[0,L],[hc*(1+1e-7)],rtol=1e-10,atol=1e-13)
            self.assertTrue(ode.success)
            h=ode.y[0,-1];H=h+q*q/(2*g*h*h)
            self.assertAlmostEqual(upstream_head(q,L,n)/H,1.,places=8)
    def test_inverse(self):
        for q in (.001,.03,.1,.3):
            H=upstream_head(q,2.,.02)
            self.assertAlmostEqual(unit_discharge(H,2.,.02)/q,1.,places=10)
    def test_equal_lanes_split_equal(self):
        p=partition(.07,[Lane(0,.5),Lane(0,.5)])
        self.assertAlmostEqual(p['fractions'][0],.5,places=10)
    def test_flow_conservation(self):
        p=partition(.07,[Lane(-.08,.5),Lane(.12,.25),Lane(0,.5)])
        self.assertLess(abs(p['mass_residual']),1e-11)
    def test_friction_reduces_capacity_at_same_head(self):
        self.assertGreater(unit_discharge(.2,2.,0.),unit_discharge(.2,2.,.02))
    def test_stable_small_length_limit(self):
        for u in (1e-8,1e-6,1e-4):self.assertAlmostEqual(_length_function(u)/(1.5*u*u),1.,delta=2e-4)
    def test_dry_lane(self):
        self.assertEqual(unit_discharge(-.2,2.,.02),0)
if __name__=='__main__':unittest.main()
