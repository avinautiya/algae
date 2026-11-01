import unittest
import math
from transmission import Route, route_discharge_m3_s as q, solve_common_stage, hydrostatic_diagnostics


class TransmissionTests(unittest.TestCase):
    def test_floor_passivity_and_positive_transmission(self):
        r = Route('rock', 10, 2, 3, 0, aperture_height_m=0.5, discharge_coefficient=0.6)
        self.assertEqual(q(r, 9, 0), 0)
        self.assertEqual(q(r, 10, 0), 0)
        self.assertEqual(q(r, 11, 12), 0)
        self.assertGreater(q(r, 11, 0), 0)

    def test_monotonic_and_zero_transmission(self):
        for r in (Route('aperture', 1, 2, 3, 1.7, 0.5, 0.6),
                  Route('darcy', 1, 2, 3, 1.7, darcy_conductance_m2_s=0.2)):
            flows = [q(r, i/20, 1.8) for i in range(150)]
            self.assertTrue(all(a <= b for a, b in zip(flows, flows[1:])))
        self.assertEqual(q(Route('sealed', 1, 2, 3, 0), 100, 0), 0)

    def test_conservation_and_bound(self):
        routes = [Route('a', 1, 2, 3, 1.7), Route('b', 2, 1, 2, 1.3)]
        solution = solve_common_stage(routes, 10, 0)
        self.assertAlmostEqual(sum(solution.route_flows_m3_s.values()), 10, places=8)
        self.assertAlmostEqual(solution.residual_m3_s, 0, places=8)
        with self.assertRaises(ValueError):
            solve_common_stage(routes, 10, 0, max_stage_m=3)
        with self.assertRaises(ValueError):
            solve_common_stage([Route('a', 0, 1, 1, 0), Route('b', 0, 1, 1, 0)], 1, 0)
        self.assertEqual(solve_common_stage(routes, 0, 0).residual_m3_s, 0)

    def test_unknown_curvature_and_indicators(self):
        d = hydrostatic_diagnostics(2, 100, 3, 0.02)
        self.assertIsNone(d['curvature_acceleration_over_g'])
        self.assertEqual(d['depth_to_length'], 0.02)
        self.assertTrue(any('unknown' in w for w in d['warnings']))
        known = hydrostatic_diagnostics(2, 100, 3, 0.02, 10)
        self.assertAlmostEqual(known['curvature_acceleration_over_g'], 9/(9.80665*10))
        self.assertNotIn('validated', known)

    def test_aperture_against_independent_vertical_quadrature(self):
        r=Route('aperture',0.,2.,.7,0.,.6,.8)
        for stage,receiver in ((.3,-1.),(1.,.2),(1.,.8),(1.,1.1)):
            a=min(.6,stage);n=20000;dz=a/n
            integral=sum(math.sqrt(2*9.80665*max(0.,stage-max(receiver,(i+.5)*dz))) for i in range(n))*dz
            self.assertAlmostEqual(q(r,stage,receiver),.8*.7*integral,delta=1e-6)

    def test_tiny_aperture_retains_small_flow(self):
        r=Route('tiny',0.,2.,1.,0.,1e-12,1.)
        self.assertAlmostEqual(q(r,1.,-1.)/(1e-12*math.sqrt(2*9.80665)),1.,places=10)
        with self.assertRaises(ValueError):Route('invalid',0.,True,1.,1.)
        with self.assertRaises(ValueError):solve_common_stage([r,Route('b',0.,0.,1.,1.)],1.,-1.,max_iterations=1.5)


if __name__ == '__main__':
    unittest.main()
