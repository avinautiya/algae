import math
import unittest
import numpy as np
from scipy.integrate import solve_ivp
from recovery_budget import static_bound, deadline_certificate


class Tests(unittest.TestCase):
    def params(self, **kw):
        p = dict(area_m2=10., coefficient=.7, exponent=1.5,
                 bypass_m=0., toe_m=.2, initial_stage_m=.05,
                 peak_m3_s=.1)
        p.update(kw)
        return p

    def test_independent_ode(self):
        p = self.params()
        b = static_bound(**p)
        def event(t, y): return y[0]-p['toe_m']
        event.terminal = True
        r = solve_ivp(lambda t, y: [(p['peak_m3_s']-p['coefficient']*
            max(y[0]-p['bypass_m'], 0.)**p['exponent'])/p['area_m2']],
            (0, b.minimum_time_s*2), [p['initial_stage_m']],
            events=event, rtol=1e-11, atol=1e-13)
        self.assertAlmostEqual(r.t_events[0][0]/b.minimum_time_s, 1., places=9)

    def test_linear_exact(self):
        p = self.params(exponent=1.)
        b = static_bound(**p)
        # Increase peak above linear critical discharge.
        self.assertFalse(b.reachable_at_peak)
        p['peak_m3_s'] = .3
        b = static_bound(**p)
        exact = p['area_m2']/p['coefficient']*math.log(
            (p['peak_m3_s']-p['coefficient']*p['initial_stage_m']) /
            (p['peak_m3_s']-p['coefficient']*p['toe_m']))
        self.assertAlmostEqual(b.minimum_time_s, exact, places=10)

    def test_equal_peak_is_not_finite_contact(self):
        p = self.params()
        p['peak_m3_s'] = p['coefficient']*p['toe_m']**p['exponent']
        self.assertFalse(static_bound(**p).reachable_at_peak)

    def test_negative_datum_invariance(self):
        a = static_bound(**self.params())
        b = static_bound(**self.params(bypass_m=-8., toe_m=-7.8,
                                      initial_stage_m=-7.95))
        self.assertAlmostEqual(a.minimum_time_s, b.minimum_time_s, places=10)

    def test_below_bypass_part(self):
        p = self.params(initial_stage_m=-.05)
        a = static_bound(**p)
        b = static_bound(**self.params(initial_stage_m=0.))
        self.assertAlmostEqual(a.minimum_time_s-b.minimum_time_s, 5., places=9)

    def test_volume_time_gate(self):
        b = static_bound(**self.params())
        self.assertTrue(deadline_certificate(b, available_volume_m3=b.minimum_volume_m3,
                                              available_time_s=b.minimum_time_s))
        self.assertFalse(deadline_certificate(b, available_volume_m3=.99*b.minimum_volume_m3,
                                               available_time_s=10*b.minimum_time_s))

    def test_initial_contact(self):
        b = static_bound(**self.params(initial_stage_m=.3))
        self.assertEqual(b.minimum_time_s, 0.)

    def test_unsteady_inflows_and_eroding_bypass_obey_bound(self):
        p = self.params()
        b = static_bound(**p)
        for erosion in (0., .00002, .0001):
            def rhs(t, y):
                inflow = .1*(.8+.2*math.sin(t/9.)**2)
                outflow = .7*max(y[0]+erosion*t, 0.)**1.5
                return [(inflow-outflow)/10., inflow]
            def event(t, y): return y[0]-.2
            event.terminal = True
            event.direction = 1
            r = solve_ivp(rhs, (0, b.minimum_time_s*5), [.05, 0.],
                          events=event, rtol=1e-10, atol=1e-12)
            self.assertEqual(len(r.t_events[0]), 1)
            self.assertGreaterEqual(r.t[-1], b.minimum_time_s)
            self.assertGreaterEqual(r.y[1,-1], b.minimum_volume_m3)

    def test_validation(self):
        for kw in (dict(area_m2=0), dict(exponent=-1), dict(toe_m=np.nan),
                   dict(coefficient=True), dict(peak_m3_s=-1)):
            with self.assertRaises(ValueError): static_bound(**self.params(**kw))


if __name__ == '__main__': unittest.main()
