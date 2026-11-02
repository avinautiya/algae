import unittest,math
from storage_identity import Network,exact_segment,rk4_segment,screening_bound

class IdentityTests(unittest.TestCase):
    def test_independent_rk4(self):
        n=Network(100000,10,200,50,1000)
        a=exact_segment(n,0.1,0.5,120,2000);b=rk4_segment(n,0.1,0.5,120,2000,400)
        self.assertAlmostEqual(a['temperature_C'],b['temperature_C'],places=10)
        self.assertLess(abs(a['ice_J']-b['ice_J']),2e-6)
        self.assertLess(abs(a['water_J']-b['water_J']),2e-6)
    def test_timevarying_inputs_telescope_identity(self):
        n=Network(100000,10,200,50,1000);t=.2;ice=qs=water=air=0
        for duration,tw,qa in [(7200,.5,100),(3600,.2,0),(7200,.1,200),(2000,0,0)]:
            r=exact_segment(n,t,tw,qa,duration);t=r['temperature_C'];ice+=r['ice_J'];qs+=r['instantaneous_ice_J'];water+=r['water_J'];air+=qa*duration
        du=n.capacity_J_K*(t-.2);gamma=n.diagnostics()['terminal_storage_to_ice_factor']
        self.assertLess(abs(qs-ice-gamma*du),1e-7)
        self.assertLess(abs(water+air-ice-du),1e-7)
    def test_gamma_bounded(self):
        for gb in [0,.1,100,100000]:
            for gd in [.01,10,100000]:
                for F in [.1,100,1000000]:
                    g=Network(1,1,gb,gd,F).diagnostics()['terminal_storage_to_ice_factor']
                    self.assertGreaterEqual(g,0);self.assertLessEqual(g,1+1e-14)
    def test_dry_limit(self):
        n=Network(100,10,0,0,0);r=exact_segment(n,2,0,0,1000)
        self.assertAlmostEqual(r['ice_J'],200,places=9)
        self.assertEqual(r['water_J'],0)
        self.assertAlmostEqual(n.diagnostics()['terminal_storage_to_ice_factor'],1)
    def test_bound_footprint(self):
        n=Network(100000,10,200,50,1000)
        a=screening_bound(n,2,1);b=screening_bound(n,2,2)
        self.assertAlmostEqual(a['equivalent_lowering_bound_m'],2*b['equivalent_lowering_bound_m'])
    def test_order_effect_vanishes_after_common_relaxation(self):
        n=Network(10000,10,100,100,500)
        def run(seq):
            t=0;ice=0
            for dt,tw,q in seq:
                r=exact_segment(n,t,tw,q,dt);t=r['temperature_C'];ice+=r['ice_J']
            return t,ice
        a=[(3600,.5,200),(3600,.1,0)];b=list(reversed(a))
        ta,ia=run(a);tb,ib=run(b)
        self.assertGreater(abs(ia-ib),1)
        relax=40*n.diagnostics()['relaxation_time_s']
        ta,ia=run(a+[(relax,0,0)]);tb,ib=run(b+[(relax,0,0)])
        self.assertLess(abs(ia-ib),1e-7)
    def test_invalid(self):
        for values in [(0,1,1,1,1),(1,-1,1,1,1),(1,1,1,1,0),(True,1,1,1,1)]:
            with self.assertRaises(ValueError):Network(*values)

if __name__=='__main__':unittest.main()
