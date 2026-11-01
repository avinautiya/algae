"""M2 regression tests for the five reproduced M1 defects."""
import math, os, tempfile, unittest, json
from ice_cliff_survival import geometry, cfd_contract as C, UnsupportedRun
from ice_cliff_survival.config import ToeConfig, Provenance
from ice_cliff_survival.runner import ToeRun
from ice_cliff_survival.closures import ConstantSources, ConstantCapacity, LinearCapacity, TableSources
from ice_cliff_survival.units import RHO_W, RHO_I, LF, K_W
def ray_thickness(z, zf, th_deg, be_deg):
    """independent explicit ray-line intersection in coordinates (x, z). Ice cliff, open side to the right: face line through (-z cot th, z), i.e. rising to the LEFT (back-tilted, ice behind it); ground z=0.
    Debris pile: surface line through the face point at height zf, descending to the right at angle beta: Z = zf - tan(be) (X - xf), xf = -zf cot th. Unit normal into the open side n=(sin th, cos th).
    d solves P + d n on the pile line (linear system, no trigonometric identities used)."""
    th, be = math.radians(th_deg), math.radians(be_deg)
    if z >= zf: return 0.0
    cot = 0.0 if th_deg == 90.0 else math.cos(th) / math.sin(th); xf = -zf * cot; P = (-z * cot, z); n = (math.sin(th), math.cos(th)); m = math.tan(be)
    # Z - zf + m (X - xf) = 0 with (X, Z) = P + d n  ->  d (n_z + m n_x) = zf - P_z - m (P_x - xf)
    return (zf - P[1] - m * (P[0] - xf)) / (n[1] + m * n[0])
class TestWedgeThickness(unittest.TestCase):
    def test_reported_repro(self):
        self.assertAlmostEqual(geometry.debris_normal_thickness(0.003, 0.012, 62.0, 33.0), 0.00565014447451, 12)
        self.assertAlmostEqual(geometry.debris_normal_thickness_v1_projected_horizontal(0.003, 0.012, 62.0, 33.0), 0.00801133650308, 12)   # old value kept for provenance
    def test_matches_ray_line_intersection(self):
        for th, be in ((62, 33), (45, 20), (80, 35), (89, 10), (30, 29), (90, 33), (90, 5)):
            for z in (0.0, 0.001, 0.003, 0.0119):
                self.assertAlmostEqual(geometry.debris_normal_thickness(z, 0.012, th, be), ray_thickness(z, 0.012, th, be), 12, f"{th},{be},{z}")
    def test_vertical_face_closed_form(self):
        self.assertAlmostEqual(geometry.debris_normal_thickness(0.002, 0.012, 90.0, 33.0), 0.010 / math.tan(math.radians(33)), 14)
    def test_crest_and_above_zero_and_linear(self):
        self.assertEqual(geometry.debris_normal_thickness(0.012, 0.012, 62, 33), 0.0); self.assertEqual(geometry.debris_normal_thickness(0.5, 0.012, 62, 33), 0.0)
        d = [geometry.debris_normal_thickness(z, 0.012, 62, 33) for z in (0.0, 0.004, 0.008)]; self.assertAlmostEqual(d[0] - d[1], d[1] - d[2], 14)
    def test_invalid_inputs(self):
        for th, be in ((33, 62), (62, 62), (62, 0), (62, -1), (91, 33), (0, 0), (math.nan, 33), (62, math.inf)):
            with self.assertRaises(ValueError): geometry.debris_normal_thickness(0.003, 0.012, th, be)
        for z, zf in ((math.nan, 0.01), (0.0, math.inf), (0.0, -1.0)):
            with self.assertRaises(ValueError): geometry.debris_normal_thickness(z, zf, 62, 33)
    def test_shield_conductance_uses_corrected_thickness(self):
        self.assertAlmostEqual(geometry.shield_conductance(0.003, 0.012, 62, 33, 0.5), 0.5 / 0.00565014447451, 8)
class TestFiniteValidation(unittest.TestCase):
    def rec(self, **kw):
        base = dict(t0=0.0, t1=1.0, q_wall_W=0.0, wetted_area_m2=1.0, melt_volume_ice_m3=0.0, meltwater_volume_m3=0.0, sensible_in_J=0.0, sensible_out_J=0.0, storage_change_J=0.0, dissipation_J=0.0); base.update(kw); return C.ExchangeRecord(**base)
    def test_nonfinite_rejected_in_every_field(self):
        for f in ("t0", "t1", "q_wall_W", "wetted_area_m2", "melt_volume_ice_m3", "meltwater_volume_m3", "sensible_in_J", "sensible_out_J", "storage_change_J", "dissipation_J"):
            for bad in (math.nan, math.inf, -math.inf):
                with self.assertRaises(ValueError, msg=f"{f}={bad}"): C.validate_records([self.rec(**{f: bad})])
                with self.assertRaises(ValueError): C.diagnose_records([self.rec(**{f: bad})])
    def test_bad_tolerances_rejected(self):
        for kw in (dict(rtol=math.nan), dict(latent_rtol=-1.0), dict(energy_rtol=math.inf), dict(latent_rtol="x")):
            with self.assertRaises(ValueError): C.validate_records([self.rec()], **kw)
    def test_signed_storage_and_flux_allowed_diagnostics_do_not_enforce(self):
        r = self.rec(storage_change_J=-5.0, sensible_in_J=0.0, sensible_out_J=5.0); d = C.diagnose_records([r]); self.assertAlmostEqual(d[0]["water_energy_residual_J"], 0.0)
        bad = self.rec(sensible_in_J=7.0); d = C.diagnose_records([bad]); self.assertEqual(d[0]["water_energy_residual_J"], 7.0)       # diagnostic pass reports, does not raise
        with self.assertRaises(ValueError): C.validate_records([bad])                                                                       # gate raises
class TestNeumannFixture(unittest.TestCase):
    def test_incompatible_densities_rejected(self):
        with self.assertRaises(ValueError): C.neumann_fixture(0.2, rho_i=917.0)
        with self.assertRaises(ValueError): C.neumann_fixture(0.2, rho_w=1000.0, rho_i=917.0)
    def test_records_equal_independent_flux_integral_and_latent(self):
        Ste = 0.2; lam, Tw, recs = C.neumann_fixture(Ste, n=7, t_end=3600.0); rho, cp = RHO_W, 4217.0; kap = K_W / (rho * cp)
        for r in recs:          # independent Gauss-Legendre-free midpoint quadrature of the analytic conductive flux at the front, substitution t=t0+(t1-t0)u^2 not needed when t0>0; handle t0=0 with t=t1 u^2
            N = 200000; q = 0.0
            for i in range(N):
                u = (i + 0.5) / N; t = r.t0 + (r.t1 - r.t0) * u if r.t0 > 0 else r.t1 * u * u; w = (r.t1 - r.t0) / N if r.t0 > 0 else 2 * r.t1 * u / N
                q += K_W * Tw * math.exp(-lam ** 2) / (math.erf(lam) * math.sqrt(math.pi * kap * t)) * w
            self.assertAlmostEqual(r.heat_J() / q, 1.0, 5); self.assertAlmostEqual(rho * LF * r.melt_volume_ice_m3 / q, 1.0, 5)
        self.assertTrue(C.validate_records(recs, latent_rtol=1e-9, energy_rtol=1e-9))
    def test_no_rescale_emits_own_analytic_heat(self):
        lam, Tw, recs = C.neumann_fixture(0.05); self.assertAlmostEqual(sum(r.heat_J() for r in recs) / (RHO_W * LF * sum(r.melt_volume_ice_m3 for r in recs)), 1.0, 12)
    def test_hot_wall_heat_and_liquid_storage_independent(self):
        Ste = 0.2; lam, Tw, recs = C.neumann_fixture(Ste, n=4, t_end=3600.0); rho, cp = RHO_W, 4217.0; kap = K_W / (rho * cp)
        T = lambda x, t: Tw * (1 - math.erf(x / (2 * math.sqrt(kap * t))) / math.erf(lam))
        def U(t):       # independent midpoint integration of rho cp T over the liquid 0<x<s(t)
            if t == 0: return 0.0
            s_ = 2 * lam * math.sqrt(kap * t); N = 200000; return rho * cp * sum(T((i + 0.5) / N * s_, t) for i in range(N)) * s_ / N
        def hot(r):     # integrate analytic wall conductive flux  -k dT/dx|_0 = k Tw/(erf(lam) sqrt(pi kap t)), substitution to remove the t^-1/2 singularity
            N = 200000; q = 0.0
            for i in range(N):
                u = (i + 0.5) / N; t = r.t0 + (r.t1 - r.t0) * u if r.t0 > 0 else r.t1 * u * u; w = (r.t1 - r.t0) / N if r.t0 > 0 else 2 * r.t1 * u / N
                q += K_W * Tw / (math.erf(lam) * math.sqrt(math.pi * kap * t)) * w
            return q
        for r in recs:
            self.assertAlmostEqual(r.sensible_in_J / hot(r), 1.0, 5); self.assertAlmostEqual(r.storage_change_J / (U(r.t1) - U(r.t0)), 1.0, 5)
            self.assertNotEqual(r.sensible_in_J, r.heat_J()); self.assertAlmostEqual(r.water_energy_residual() / r.sensible_in_J, 0.0, 12)
class TestSolveLambda(unittest.TestCase):
    """independent check of cfd/stefan2d.solve_lambda: build the similarity solution T(x,t) and test the PDE (with liquid velocity u=(1-r) sdot), both boundary values and the Stefan condition by finite differences."""
    def test_solution_satisfies_pde_and_stefan_condition(self):
        from ice_cliff_survival.cfd.stefan2d import solve_lambda
        cp, kap = 4217.0, K_W / (RHO_W * 4217.0)
        for Ste in (0.01, 0.2, 1.0):
            for r in (1.0, RHO_I / RHO_W):
                lam, b = solve_lambda(Ste, r, True); self.assertAlmostEqual(b, (1 - r) * lam, 14); Tw = Ste * LF / cp; E = lambda t, x: Tw * (1 - (math.erf(x / (2 * math.sqrt(kap * t)) - b) + math.erf(b)) / (math.erf(lam - b) + math.erf(b)))
                t = 100.0; s = lambda t: 2 * lam * math.sqrt(kap * t); x = 0.5 * s(t); h = 1e-2 * s(t); dt = 1e-2
                Tt = (E(t + dt, x) - E(t - dt, x)) / (2 * dt); Tx = (E(t, x + h) - E(t, x - h)) / (2 * h); Txx = (E(t, x + h) - 2 * E(t, x) + E(t, x - h)) / h ** 2; sd = lam * math.sqrt(kap / t); u = (1 - r) * sd
                self.assertLess(abs(Tt + u * Tx - kap * Txx) / abs(kap * Txx), 2e-4)
                self.assertAlmostEqual(E(t, 0.0) / Tw, 1.0, 12); self.assertAlmostEqual(E(t, s(t)), 0.0, 9)
                Tx_s = (3 * E(t, s(t)) - 4 * E(t, s(t) - h) + E(t, s(t) - 2 * h)) / (2 * h); rho_i = r * RHO_W; self.assertLess(abs(-K_W * Tx_s / (rho_i * LF * sd) - 1.0), 5e-5)
class TestToeRunSafeguards(unittest.TestCase):
    cfg = ToeConfig(W=10.0, porosity=0.3, H_cliff=15.0, theta_deg=60.0, beta_deg=35.0); P = Provenance("c", "synthetic")
    def mk(self, S=3e-5, Q=1e-5, dt=600.0, **kw): return ToeRun(self.cfg, ConstantSources({"S_mantle": S}), ConstantCapacity(Q), dict(sources=self.P, capacity=self.P), dt, **kw)
    def test_dt_rejected(self):
        for dt in (0.0, -1.0, math.nan, math.inf, -math.inf, True, "1"):
            with self.assertRaises(ValueError): self.mk(dt=dt)
    def test_horizons(self):
        r = self.mk(); r.advance(1200.0)
        for h in (math.nan, math.inf, -math.inf): 
            with self.assertRaises(ValueError): r.advance(h)
        with self.assertRaises(ValueError): r.advance(600.0)          # backward
        n = len(r.history); r.advance(1200.0); self.assertEqual(len(r.history), n)      # unchanged horizon: no-op
    def test_opaque_callables_refused_unless_unsafe(self):
        with self.assertRaises(UnsupportedRun): ToeRun(self.cfg, lambda t: {}, ConstantCapacity(0.0), dict(sources=self.P, capacity=self.P), 600.0)
        r = ToeRun(self.cfg, lambda t: {"S_mantle": 1e-5}, lambda t, V: 1e-5, dict(sources=self.P, capacity=self.P), 600.0, unsafe_opaque=True); r.advance(1200.0); self.assertFalse(r.restart_verified)
        p = os.path.join(tempfile.mkdtemp(), "ck.json")
        with self.assertRaises(UnsupportedRun): r.save(p)
        r.save(p, unsafe_restart=True); self.assertFalse(json.load(open(p))["restart_verified"])
        r2 = ToeRun(self.cfg, lambda t: {"S_mantle": 9e-5}, lambda t, V: 1e-5, dict(sources=self.P, capacity=self.P), 600.0, unsafe_opaque=True)
        with self.assertRaises(UnsupportedRun): r2.load(p)
        r2.load(p, unsafe_restart=True); self.assertFalse(r2.restart_verified)      # accepted but flagged; a changed callable is undetectable, which is exactly why it is unverified
        reg = self.mk()
        with self.assertRaises(ValueError): reg.load(p, unsafe_restart=True)         # registry run cannot adopt an opaque checkpoint
    def ck(self, run, t=3000.0):
        run.advance(t); p = os.path.join(tempfile.mkdtemp(), "ck.json"); run.save(p); return p
    def test_changed_coefficient_refused(self):
        p = self.ck(self.mk())
        with self.assertRaises(ValueError): self.mk(S=3.0000001e-5).load(p)
        with self.assertRaises(ValueError): self.mk(Q=2e-5).load(p)
    def test_changed_closure_type_refused(self):
        p = self.ck(self.mk()); r = ToeRun(self.cfg, ConstantSources({"S_mantle": 3e-5}), LinearCapacity(1e-5, 0.0), dict(sources=self.P, capacity=self.P), 600.0)
        with self.assertRaises(ValueError): r.load(p)
    def test_changed_forcing_before_checkpoint_refused_but_extension_allowed(self):
        e = [0.0, 3600.0, 7200.0]; rates = [{"S_mantle": 3e-5}, {"S_mantle": 5e-5}]
        mk = lambda e, rates: ToeRun(self.cfg, TableSources(e, rates), ConstantCapacity(1e-5), dict(sources=self.P, capacity=self.P), 600.0)
        a = mk(e, rates); p = self.ck(a, 5400.0)
        with self.assertRaises(ValueError): mk(e, [{"S_mantle": 3.1e-5}, {"S_mantle": 5e-5}]).load(p)            # change in the already-simulated interval 0
        with self.assertRaises(ValueError): mk(e, [{"S_mantle": 3e-5}, {"S_mantle": 5.5e-5}]).load(p)            # change in the interval containing the checkpoint time
        b = mk([0.0, 3600.0, 7200.0, 10800.0], rates + [{"S_mantle": 1e-5}]); b.load(p); b.advance(9000.0)          # permitted: extension after the checkpoint
        ref = mk([0.0, 3600.0, 7200.0, 10800.0], rates + [{"S_mantle": 1e-5}]); ref.advance(9000.0)
        self.assertEqual(b.stock.Vs, ref.stock.Vs); self.assertEqual(b.stock.cum, ref.stock.cum)
        with self.assertRaises(ValueError): mk(e, rates).load(p) if False else mk([0.0, 3600.0], rates[:1]).load(p)    # table that cannot cover the checkpoint time
    def test_missing_and_malformed_identity(self):
        p = self.ck(self.mk()); d = json.load(open(p)); q = p + ".bad"
        for mut in (lambda d: d["ident"].pop("cfg"), lambda d: d["ident"].pop("closures"), lambda d: d["ident"]["closures"].pop("sources"), lambda d: d.pop("ident"), lambda d: d["stock"].pop("t"), lambda d: d.pop("restart_verified"), lambda d: d["stock"].__setitem__("t", math.nan), lambda d: d["stock"].__setitem__("t", "x"), lambda d: d.pop("Vs0")):
            dd = json.loads(json.dumps(d)); mut(dd); json.dump(dd, open(q, "w"))
            with self.assertRaises(ValueError): self.mk().load(q)
        open(q, "w").write("[1,2]")
        with self.assertRaises(ValueError): self.mk().load(q)
    def test_unchanged_restart_equivalence(self):
        a = self.mk(); a.advance(36000.0); p = self.ck(self.mk(), 18000.0); c = self.mk(); c.load(p); c.advance(36000.0)
        self.assertEqual(a.stock.Vs, c.stock.Vs); self.assertEqual(a.stock.cum, c.stock.cum)
    def test_registry_roundtrip_and_unknown_spec(self):
        from ice_cliff_survival import closures
        t = TableSources([0.0, 1.0, 2.0], [{"S_up": 1e-6}, {"S_up": 2e-6}]); self.assertEqual(closures.from_spec(t.spec()).spec(), t.spec())
        for bad in ({}, {"type": "nope", "params": {}}, {"type": "constant_capacity"}): 
            with self.assertRaises(ValueError): closures.from_spec(bad)
if __name__ == "__main__": unittest.main()
class TestShieldPlanar(unittest.TestCase):
    def test_series_matches_independent_fv_and_late_limit(self):
        from ice_cliff_survival import shield_planar as S
        d, k, rc, Tw = 0.007, 0.5, 1.5e6, 1.2
        for t in (60.0, 600.0, 3600.0): self.assertAlmostEqual(float(S.Q_series(t, d, k, rc, Tw)) / S.fv_numeric(t, d, k, rc, Tw, N=200, nsteps=2000), 1.0, 4)
        t = 1e5; self.assertAlmostEqual((float(S.Q_massless(t, d, k, Tw)) - float(S.Q_series(t, d, k, rc, Tw))) / (rc * d * d * Tw / 6.0), 1.0, 9)
        self.assertGreater(float(S.Q_massless(60.0, d, k, Tw)), float(S.Q_series(60.0, d, k, rc, Tw)))      # massless always over-delivers heat
