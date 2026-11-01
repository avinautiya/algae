import math, random, unittest
from ice_cliff_survival import scalar_model as sm
from ice_cliff_survival.toe_stock import ToeStock
def family(rnd):
    a, p, ell, c, m, D, S0 = rnd.uniform(0.5, 3), rnd.uniform(0.5, 2), rnd.uniform(0.5, 4), rnd.uniform(0.2, 2), rnd.uniform(1.0, 2.5), rnd.uniform(0.1, 1.5), rnd.uniform(0.0, 0.5)
    F = lambda T, s: a * T ** p * math.exp(-s / ell); C = lambda s, Q: c * Q * s ** m
    return sm.ScalarModel(F=F, C=C, D=D, S0=lambda t: S0), (a, p, ell, c, m, D, S0)
class TestScalar(unittest.TestCase):
    def test_equilibrium_sensitivities_match_finite_difference(self):
        rnd = random.Random(3); n = 0
        for _ in range(60):
            m, _p = family(rnd); Q = 1.0; T = rnd.uniform(0.5, 3.0); s = sm.equilibrium(m, T, Q, 0.0, 50.0)
            if s is None: continue
            n += 1; sens = sm.sensitivities(m, T, Q, s); self.assertLess(sens["g_s"], 0); h = 1e-5
            s_p, s_m = sm.equilibrium(m, T + h, Q, 0.0, 50.0), sm.equilibrium(m, T - h, Q, 0.0, 50.0); R_p, R_m = m.F(T + h, s_p), m.F(T - h, s_m)
            self.assertAlmostEqual((s_p - s_m) / (2 * h), sens["ds_dT"], 5); self.assertAlmostEqual((R_p - R_m) / (2 * h), sens["dR_dT"], 5); self.assertGreaterEqual(sens["dR_dT"], 0.0); self.assertGreater(sens["ds_dT"], 0.0)
        self.assertGreater(n, 30)
    def test_uniqueness_monotone_g(self):
        rnd = random.Random(5)
        for _ in range(20):
            m, _p = family(rnd); T = 1.7; g = [m.rhs(0.0, 0.05 * i, T, 1.0) for i in range(400)]; self.assertTrue(all(g[i + 1] < g[i] for i in range(len(g) - 1)))
    def test_hypothesis_violation_can_reverse_sign(self):
        # C_s < 0 (export decreasing with stock; violates H) with C_s - D F_s > 0: dR*/dT = F_T C_s/(C_s - D F_s) < 0
        F_T, F_s, C_s, D = 1.0, -2.0, -0.5, 1.0; den = C_s - D * F_s; self.assertGreater(den, 0); self.assertLess(F_T * C_s / den, 0)
    def test_cumulative_response_nonnegative_and_identity(self):
        rnd = random.Random(11)
        for _ in range(12):
            m, _p = family(rnd); TL = lambda t: 1.0 + 0.5 * math.sin(2 * math.pi * t / 5.0); TH = lambda t: TL(t) + 0.3 + 0.2 * math.sin(7.0 * t) ** 2; Q = lambda t: 1.0
            a = sm.integrate(m, TH, Q, 0.3, 10.0, 4000); b = sm.integrate(m, TL, Q, 0.3, 10.0, 4000)
            self.assertGreaterEqual(a["s"] - b["s"], -1e-9); lhs = m.D * (a["int_R"] - b["int_R"]); rhs = (a["s"] - b["s"]) + (a["int_C"] - b["int_C"]) - 0.0
            self.assertAlmostEqual(lhs, rhs, 6); self.assertGreaterEqual(lhs, -1e-9)
    def test_clipping_is_refused(self):
        m = sm.ScalarModel(F=lambda T, s: 0.0, C=lambda s, Q: 1.0, D=1.0, S0=lambda t: 0.0)       # C(0) != 0 violates admissibility: s goes negative
        with self.assertRaises(ValueError): sm.integrate(m, lambda t: 1.0, lambda t: 1.0, 0.0, 1.0, 100)
    def test_accumulation_bound_random_bounded_forcing(self):
        rnd = random.Random(2)
        for _ in range(30):
            Smin, Cmax = rnd.uniform(2e-5, 6e-5), rnd.uniform(0.0, 1.9e-5); Vs0 = rnd.uniform(0, 5.0); Vcrit = Vs0 + rnd.uniform(1.0, 20.0); st = ToeStock(Vs0); dt = 3600.0; hit = None
            for i in range(20000):
                st.step(dt, {"S_ext": Smin + rnd.uniform(0, 3e-5)}, rnd.uniform(0, Cmax))
                self.assertGreaterEqual(st.Vs, Vs0 + (Smin - Cmax) * st.t - 1e-9)
                if hit is None and st.Vs >= Vcrit: hit = st.t; break
            self.assertIsNotNone(hit); self.assertLessEqual(hit, sm.accumulation_hitting_bound(Vs0, Vcrit, Smin, Cmax) + dt)
