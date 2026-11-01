import datetime as dt, json, math, os, sys, tempfile, unittest
import numpy as np
ROOT = os.path.abspath(os.path.join(os.path.dirname(__file__), ".."))
sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, "ice_cliff_survival", "reactivation", ".."))
from ice_integration import _paths, forcing, hydraulics, finite_heat, coupling
import predict_external as PE
from gap_reconstructed_temperature import GapReconstructedTemperatureSeries
from native_temperature import TemperatureSeries
T0 = dt.datetime(2019, 7, 20)
def secs(*h): return [T0 + dt.timedelta(hours=x) for x in h]
class TestForcing(unittest.TestCase):
    def test_irregular_sampling_integrates_over_elapsed_time(self):
        hs = [0, 0.001, 0.0035, 0.01, 0.2]; vals = [1.0, 3.0, 2.0, 4.0, 0.5]       # unequal spacing inside the 60 s support limit? use explicit seconds
        recs = [(T0 + dt.timedelta(seconds=s), v) for s, v in zip([0, 20, 50, 80, 120], vals)]; ser = TemperatureSeries(recs, max_gap_seconds=60)
        exact = sum((b - a) * (va + vb) / 2 for (a, va), (b, vb) in zip(zip([0, 20, 50, 80], vals[:4]), zip([20, 50, 80, 120], vals[1:])))
        r = ser.integrate(T0, T0 + dt.timedelta(seconds=120)); self.assertAlmostEqual(r["integral_C_s"], exact, 12); self.assertEqual(r["covered_seconds"], 120)
        r2 = ser.integrate(T0 + dt.timedelta(seconds=10), T0 + dt.timedelta(seconds=65)); self.assertAlmostEqual(r2["covered_seconds"], 55, 12)
    def test_offsets_exactly_once_on_real_source(self):
        from original_temperature_adapter import source_rows, OFFSETS
        for station in (3, 5):
            s = TemperatureSeries.from_archive(_paths.SWISS_SOURCE, station); rows = list(source_rows(__import__("pathlib").Path(_paths.SWISS_SOURCE) / f"products/temp{station}_corr_raw.csv"))
            raw = {r[2]: r[3] for r in rows}; from native_temperature import EPOCH
            for a, b, ta, tb in s.segments[:50]:
                t = EPOCH + dt.timedelta(seconds=a); self.assertAlmostEqual(ta, raw[t] + OFFSETS[station], 12, msg="offset must be applied exactly once (not 0x, not 2x)")
    def test_support_accounting_never_zero_fills(self):
        # 0-1 h supported (20 s spacing); 1-3 h gap (reconstructible, <=24 h); 3-3.01 h supported; 3.01 h .. 3.01 h+3 days gap (unsupported >24h) ; supported again
        recs = []; t = 0.0
        for k in range(0, 181): recs.append((T0 + dt.timedelta(seconds=20 * k), 2.0))
        for k in range(0, 31): recs.append((T0 + dt.timedelta(hours=3, seconds=20 * k), 4.0))
        end2 = T0 + dt.timedelta(hours=3, seconds=600)
        for k in range(0, 31): recs.append((end2 + dt.timedelta(days=3, seconds=20 * k), 6.0))
        ser = GapReconstructedTemperatureSeries(recs, 60, 86400)
        ref = forcing.FieldForcing.__new__(forcing.FieldForcing); ref.series = ser; ref.method = "linear_endpoint_bridge"; ref.max_gap = 86400
        acc = ref.integrate(T0, end2 + dt.timedelta(days=3, seconds=600))
        self.assertAlmostEqual(acc.observed_s + acc.reconstructed_s + acc.unsupported_s, acc.requested_s, 6); self.assertGreater(acc.reconstructed_s, 7000); self.assertGreater(acc.unsupported_s, 2 * 86400 * 0.99)
        self.assertFalse(acc.strict_native_only); mean = acc.mean_temperature(); self.assertTrue(2.0 < mean < 6.0)       # weighted by SUPPORTED seconds only
        strict = acc.mean_temperature(include_reconstructed=False); self.assertAlmostEqual(strict, (2.0 * 3600 + 4.0 * 600 + 6.0 * 600) / (3600 + 1200), 6)
class TestHydraulics(unittest.TestCase):
    def test_pressure_is_not_depth_without_clearance(self):
        with self.assertRaises(hydraulics.HydraulicRefusal): hydraulics.depth_from_pressure(0.2)
        d, b = hydraulics.depth_from_pressure(0.2, allow_proxy=True); self.assertEqual((d, b), (0.2, "pressure_proxy_zero_clearance_assumed"))
        d, b = hydraulics.depth_from_pressure(0.2, sensor_clearance_m=0.07); self.assertAlmostEqual(d, 0.27); self.assertEqual(b, "measured_bed_clearance")
        for bad in (0.0, float("nan"), -1.0): 
            with self.assertRaises(hydraulics.HydraulicRefusal): hydraulics.depth_from_pressure(bad, 0.0)
    def test_consistent_re_and_labels_and_no_station_substitution(self):
        s = hydraulics.SectionHydraulics("P5", 2.0, 0.18, 0.15, "pressure_proxy_zero_clearance_assumed", "reconstructed", "daily stage proxy"); d = s.derived()
        self.assertAlmostEqual(d["Re"] / d["Re_alt"], 1.0, 12); self.assertTrue(d["labels"]["reconstructed"])
        with self.assertRaises(hydraulics.HydraulicRefusal): hydraulics.require_station(s, "P3")
        with self.assertRaises(hydraulics.HydraulicRefusal): hydraulics.SectionHydraulics("P3", 2.0, 0.18, 0.15, "unlabelled", "measured", "x")
def sect(n, L, depth=0.2, width=0.5): return PE.Section(n, width, depth, L, False, "synthetic fixture", "synthetic_fixture")
def acc(T, s=3600.0): return forcing.SupportAccount(s, s, 0.0, 0.0, T * s, 0.0, "linear_endpoint_bridge", 0.0)
class TestFiniteHeat(unittest.TestCase):
    M = "gnielinski_smooth_unfitted"
    def test_conservation_and_separation(self):
        m = finite_heat.IntegratedHeatModel([sect("a", 5.0), sect("b", 5.0)], self.M, -0.0); m.advance(acc(2.0), 0.07); r = m.residuals()
        self.assertLess(abs(r["water_side_J"]), 1e-6 * m.ledger["incoming_water_heat_J"]); self.assertLess(abs(r["interface_J"]), 1e-6 * m.ledger["latent_J"]); self.assertIn("atmospheric exchange", r["not_implemented"])
        self.assertGreater(m.ledger["floor_heat_J"], m.ledger["wall_heat_J"] * 0)
    def test_zero_heat(self):
        m = finite_heat.IntegratedHeatModel([sect("a", 5.0)], self.M); m.advance(acc(0.0), 0.07)
        self.assertEqual(m.ledger["floor_heat_J"] + m.ledger["wall_heat_J"] + m.ledger["latent_J"], 0.0)
    def test_strong_depletion_cannot_exceed_supply(self):
        m = finite_heat.IntegratedHeatModel([sect("long", 1e6)], self.M); m.advance(acc(2.0), 0.07); L = m.ledger
        self.assertLessEqual(L["floor_heat_J"] + L["wall_heat_J"], L["incoming_water_heat_J"] * (1 + 1e-12)); self.assertLess(L["outgoing_water_heat_J"], 1e-9 * L["incoming_water_heat_J"])
    def test_section_subdivision(self):
        one = finite_heat.IntegratedHeatModel([sect("a", 12.0)], self.M); one.advance(acc(1.5), 0.05); three = finite_heat.IntegratedHeatModel([sect("a", 4.0), sect("b", 4.0), sect("c", 4.0)], self.M); three.advance(acc(1.5), 0.05)
        for k in ("floor_heat_J", "wall_heat_J", "outgoing_water_heat_J", "latent_J"): self.assertAlmostEqual(one.ledger[k] / three.ledger[k], 1.0, 12)
    def test_time_mean_temperature_exact_for_linear_law(self):
        a = finite_heat.IntegratedHeatModel([sect("a", 8.0)], self.M); a.advance(acc(1.0), 0.06); a.advance(acc(3.0), 0.06); b = finite_heat.IntegratedHeatModel([sect("a", 8.0)], self.M); b.advance(acc(2.0, 7200.0), 0.06)
        self.assertAlmostEqual(a.ledger["floor_heat_J"] / b.ledger["floor_heat_J"], 1.0, 12)
    def test_unsupported_time_gets_no_heat_and_is_not_zero_contact(self):
        m = finite_heat.IntegratedHeatModel([sect("a", 5.0)], self.M); m.advance(forcing.SupportAccount(3600.0, 0.0, 0.0, 3600.0, 0.0, 0.0, "x", 0.0), 0.07)
        self.assertEqual(m.ledger["floor_heat_J"], 0.0); self.assertEqual(m.time["unsupported_s"], 3600.0); self.assertIn("NOT inferred", m.history[-1]["note"])
    def test_strict_mode_treats_reconstructed_as_unsupported(self):
        a = forcing.SupportAccount(3600.0, 1800.0, 1800.0, 0.0, 1800 * 2.0, 1800 * 4.0, "linear_endpoint_bridge", 0.0)
        s = finite_heat.IntegratedHeatModel([sect("a", 5.0)], self.M); s.advance(a, 0.07, include_reconstructed=False); r = finite_heat.IntegratedHeatModel([sect("a", 5.0)], self.M); r.advance(a, 0.07, include_reconstructed=True)
        self.assertEqual(s.time["reconstructed_s"], 0.0); self.assertEqual(s.time["unsupported_s"], 1800.0); self.assertEqual(r.time["reconstructed_s"], 1800.0)
        self.assertLess(s.ledger["floor_heat_J"], r.ledger["floor_heat_J"])
    def test_restart_equivalence_and_fingerprint_refusal(self):
        secs_ = [sect("a", 6.0), sect("b", 6.0)]; accs = [(acc(1.0), 0.05), (acc(2.5, 1800.0), 0.08), (acc(0.5, 5400.0), 0.06)]
        full = finite_heat.IntegratedHeatModel(secs_, self.M, metadata=dict(src="x"))
        for a, q in accs: full.advance(a, q)
        part = finite_heat.IntegratedHeatModel(secs_, self.M, metadata=dict(src="x")); part.advance(*accs[0]); st = json.loads(json.dumps(part.state()))
        res = finite_heat.IntegratedHeatModel.restore(secs_, self.M, st, metadata=dict(src="x"))
        for a, q in accs[1:]: res.advance(a, q)
        for k in full.ledger: self.assertAlmostEqual(full.ledger[k], res.ledger[k], delta=1e-9 * max(1.0, abs(full.ledger[k])))
        with self.assertRaises(ValueError): finite_heat.IntegratedHeatModel.restore([sect("a", 7.0), sect("b", 6.0)], self.M, st, metadata=dict(src="x"))
        with self.assertRaises(ValueError): finite_heat.IntegratedHeatModel.restore(secs_, self.M, st, metadata=dict(src="DIFFERENT"))
class TestAgainstExistingCoupledSolver(unittest.TestCase):
    def test_predictor_equals_swe2d_in_the_plug_limit(self):
        """identical NTU => identical delivered heat: predictor P = rho cp Q T (1-exp(-NTU)) vs the coupled solver's bed-exchange sink on a flat lane (supported limiting case; not moving geometry)"""
        from ice_cliff_survival.reactivation.swe2d import Sim, RHO_W, CP
        s = sect("lane", 2.0, depth=0.3, width=0.5); p = s.width_m + 2 * s.depth_m; h = PE.H_FIXED; heff = h * p / s.width_m       # same h*p*L as h_eff*W*L
        sim = Sim(dx=0.1, z_toe=0.0, z_byp=0.0, heff=heff, Tw=0.5, tail_slope=0.5, ramp_ratio=0.0, pool_floor=-0.3); sim.hydraulics(tmax=300.0); qj = sim.Qsec(sim.toe_cols()); d = sim.diagnostics()
        measured = float((sim.last_qacc[sim.i0:sim.i1 + 1, sim.toe_cols()]).sum() * sim.dx ** 2)
        r = PE.predict_interval([PE.Section("lane", 0.5, 0.3, 2.0, False, "fixture", "synthetic_fixture")], qj, 0.5, 3600.0, 0.0, "constant_h_fitted")
        self.assertAlmostEqual(measured / (r["removed_water_heat_J"] / 3600.0), 1.0, delta=0.01)
class TestCoupledRestart(unittest.TestCase):
    def test_restart_equivalence_and_unsupported_interval(self):
        from ice_cliff_survival.reactivation import cases
        def drv(): return coupling.CoupledDriver(cases.build("strip", dx=0.1), dict(src="fixture", offsets_once=True), hyd_max=20.0)
        plan = [(acc(1.0, 6 * 3600.0), 0.07), (forcing.SupportAccount(6 * 3600.0, 0.0, 0.0, 6 * 3600.0, 0.0, 0.0, "x", 0.0), 0.07), (acc(1.5, 6 * 3600.0), 0.07), (acc(2.0, 6 * 3600.0), 0.07)]
        a = drv()
        for x in plan: a.advance(*x)
        b = drv(); b.advance(*plan[0]); b.advance(*plan[1]); path = os.path.join(tempfile.mkdtemp(), "ck.npz"); b.checkpoint(path)
        c = drv(); c.restore(path)
        for x in plan[2:]: c.advance(*x)
        self.assertTrue(np.allclose(a.sim.zb, c.sim.zb, atol=1e-12, rtol=0)); self.assertEqual(a.time["unsupported_s"], c.time["unsupported_s"]); self.assertEqual(a.log[1]["advanced"], False)
        z0 = drv(); zb0 = z0.sim.zb.copy(); z0.advance(forcing.SupportAccount(3600.0, 0, 0, 3600.0, 0, 0, "x", 0.0), 0.07); self.assertTrue(np.array_equal(zb0, z0.sim.zb), "unsupported interval must not alter the bed")
        d = drv()
        with self.assertRaises(ValueError): 
            d2 = coupling.CoupledDriver(cases.build("strip", dx=0.1), dict(src="DIFFERENT"), hyd_max=20.0); d2.restore(path)
if __name__ == "__main__": unittest.main()
