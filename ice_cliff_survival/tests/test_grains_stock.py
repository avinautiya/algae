import math, random, unittest
from ice_cliff_survival import grains, toe_stock, budgets
from ice_cliff_survival.toe_stock import ToeStock
class TestGrains(unittest.TestCase):
    def test_concentration_conversion_roundtrip_and_value(self):
        self.assertAlmostEqual(grains.c_mix_from_c_ice(0.1), 0.1 / 1.1, 14); self.assertAlmostEqual(grains.c_ice_from_c_mix(grains.c_mix_from_c_ice(0.37)), 0.37, 14)
    def test_englacial_release_equals_independent_volume_bookkeeping(self):
        # independent: sweep V_mix = A r dt of mixture; pure ice melted = (1-c_mix) V_mix; grains = c_ice * V_ice
        c_mix, A, r, dt = 0.04, 80.0, 2e-6, 86400.0; V_mix = A * r * dt; V_ice = (1 - c_mix) * V_mix; grains_ind = grains.c_ice_from_c_mix(c_mix) * V_ice
        self.assertAlmostEqual(grains.englacial_source(c_mix, A, r) * dt, grains_ind, places=15)
    def test_relative_flux_frame_equals_source(self):
        # frozen grains (u_s = 0), boundary moving INTO ice at speed r => along the inward (toward-water) normal the boundary velocity is -r
        c, A, r = 0.05, 12.0, 3e-6; self.assertAlmostEqual(grains.relative_flux(c, 0.0, -r, A), grains.englacial_source(c, A, r), 18)
    def test_solid_bulk_conversion(self):
        self.assertAlmostEqual(grains.solid_from_bulk(10.0, 0.3), 7.0); self.assertAlmostEqual(grains.bulk_from_solid(7.0, 0.3), 10.0)
    def test_mantle_source_units_and_solid_basis(self):
        # bulk throughput of a strip: W d v (m3 bulk/s); solid = (1-phi)*that; delivery fraction scales it
        W, d, v, phi, f = 20.0, 0.3, 1e-6, 0.35, 0.8; bulk = W * d * v; self.assertAlmostEqual(grains.mantle_source(phi, W, d, v, f), (1 - phi) * bulk * f, 18)
        with self.assertRaises(ValueError): grains.mantle_source(1.0, W, d, v, f)
class TestMantleSplit(unittest.TestCase):
    def test_split_conserves_released_solid(self):
        d, u = grains.mantle_split(0.3, 10.0, 0.4, 2e-6, 0.65); self.assertAlmostEqual(d + u, 0.7 * 10.0 * 0.4 * 2e-6, 18); self.assertAlmostEqual(d / (d + u), 0.65, 12)
class TestToeStock(unittest.TestCase):
    def test_constant_forcing_analytic_linear_growth(self):
        st = ToeStock(Vs=2.0); S, C, dt = 5e-5, 2e-5, 600.0
        for _ in range(1000): st.step(dt, {"S_mantle": S, "S_up": 0.0}, C)
        self.assertAlmostEqual(st.Vs, 2.0 + (S - C) * 1000 * dt, 12); self.assertAlmostEqual(st.budget_residual(2.0), 0.0, 12)
    def test_zero_stock_admissible_and_export_recorded_actual(self):
        st = ToeStock(0.0); r = st.step(100.0, {"S_ext": 1e-5}, 5e-5)         # capacity > supply: everything supplied leaves, stock stays 0
        self.assertEqual(st.Vs, 0.0); self.assertAlmostEqual(r["Q_out"], 1e-5, 18); self.assertAlmostEqual(st.cum["Q_out_capacity_unused"], 4e-5 * 100.0, 15); self.assertAlmostEqual(st.budget_residual(0.0), 0.0, 15)
        r = st.step(100.0, {}, 5e-5); self.assertEqual(r["Q_out"], 0.0); self.assertEqual(st.Vs, 0.0)
    def test_export_limited_by_availability_drains_stock_exactly(self):
        st = ToeStock(1.0); st.step(10.0, {}, 1.0); self.assertEqual(st.Vs, 0.0); self.assertAlmostEqual(st.cum["Q_out_actual"], 1.0, 15)      # capacity 1 m3/s would remove 10 m3: only 1 m3 existed
    def test_no_double_count_of_mantle_in_external(self):
        a, b = ToeStock(0.0), ToeStock(0.0); mant = 3e-5; ext_total = 4e-5      # caller's external record already contains the mantle release
        a.step(1000.0, {"S_mantle": mant, "S_ext": ext_total}, 0.0, mantle_double_counted_in_ext=mant); b.step(1000.0, {"S_mantle": mant, "S_ext": ext_total - mant}, 0.0)
        self.assertAlmostEqual(a.Vs, b.Vs, 15); self.assertAlmostEqual(a.Vs, ext_total * 1000.0, 15)
        with self.assertRaises(ValueError): ToeStock().step(1.0, {"S_ext": 1e-6}, 0.0, mantle_double_counted_in_ext=2e-6)
    def test_independent_reconstruction_random_forcing(self):
        rnd = random.Random(7); st = ToeStock(0.5); t = [0.0]; srcs = []; outs = []
        for i in range(500):
            dt = rnd.uniform(10, 900); S = {k: rnd.uniform(0, 3e-5) for k in ("S_up", "S_mantle", "S_eng", "S_ext")}; C = rnd.uniform(0, 8e-5); r = st.step(dt, S, C); t.append(t[-1] + dt); srcs.append(S); outs.append(r["Q_out"])
        self.assertAlmostEqual(budgets.grain_budget(0.5, t, srcs, outs), st.Vs, 9); self.assertGreaterEqual(st.Vs, 0.0)
    def test_unknown_source_rejected(self):
        with self.assertRaises(ValueError): ToeStock().step(1.0, {"rockfall": 1.0}, 0.0)
