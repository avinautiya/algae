import math, unittest
from ice_cliff_survival import geometry, interface, cfd_contract
from ice_cliff_survival.units import RHO_I, RHO_W, LF, K_W
def shoelace(P): return 0.5 * abs(sum(P[i][0] * P[(i + 1) % len(P)][1] - P[(i + 1) % len(P)][0] * P[i][1] for i in range(len(P))))
class TestGeometry(unittest.TestCase):
    def test_wedge_area_matches_polygon(self):
        for th, be, zf in ((60, 35, 3.0), (80, 30, 7.5), (45, 20, 1.2)):
            t, b = math.radians(th), math.radians(be)
            # vertices: toe (0,0); pile outer edge (-L,0) where the repose line meets the ground; apex on the face (zf/tan th, zf) with the repose line passing through it: zf = tan(b) (x_f + L)
            xf = zf / math.tan(t); L = zf / math.tan(b) - xf; poly = [(-L, 0.0), (0.0, 0.0), (xf, zf)]
            self.assertAlmostEqual(geometry.wedge_area(zf, th, be), shoelace(poly), 10)
            self.assertAlmostEqual(geometry.wedge_height(shoelace(poly), th, be), zf, 10)
    def test_exposure_partition_and_critical_volume(self):
        W, phi, H, th, be = 10.0, 0.3, 15.0, 60.0, 35.0; Vc = geometry.critical_solid_volume(W, phi, H, th, be)
        e = geometry.exposure(0.999 * Vc, W, phi, H, th, be); self.assertFalse(e["fully_buried"]); e2 = geometry.exposure(1.0001 * Vc, W, phi, H, th, be); self.assertTrue(e2["fully_buried"]); self.assertAlmostEqual(e2["exposed_area"], 0.0, 12)
        self.assertAlmostEqual(e["exposed_area"] + e["covered_area"], W * H / math.sin(math.radians(th)), 10)
        for V in (0.1 * Vc, 0.5 * Vc): ex = geometry.exposure(V, W, phi, H, th, be, water_depth_at_toe=2.0); self.assertGreaterEqual(ex["exposed_area"], 0)
    def test_wetted_exposed_area_requires_water_above_pile(self):
        W, phi, H, th, be = 10.0, 0.3, 15.0, 60.0, 35.0; V = 0.3 * geometry.critical_solid_volume(W, phi, H, th, be); zf = geometry.exposure(V, W, phi, H, th, be)["pile_height_on_face"]
        self.assertEqual(geometry.exposure(V, W, phi, H, th, be, water_depth_at_toe=0.9 * zf)["wetted_exposed_area"], 0.0); self.assertGreater(geometry.exposure(V, W, phi, H, th, be, water_depth_at_toe=zf + 0.5)["wetted_exposed_area"], 0.0)
    def test_stock_threshold_alone_is_not_burial_with_changed_repose_angle(self):
        # same solid volume covers different fractions for different repose angles: a volume threshold is geometry dependent
        W, phi, H, th = 10.0, 0.3, 15.0, 70.0; V = 0.4 * geometry.critical_solid_volume(W, phi, H, th, 35.0)
        self.assertNotAlmostEqual(geometry.exposure(V, W, phi, H, th, 35.0)["covered_fraction"], geometry.exposure(V, W, phi, H, th, 25.0)["covered_fraction"], places=3)
class TestInterface(unittest.TestCase):
    def test_stefan_speed_matches_exact_neumann_solution_equal_density(self):
        Ste = 0.2; lam, Tw, _ = cfd_contract.neumann_fixture(Ste, n=2); kap = K_W / (RHO_W * 4217.0)
        for t in (100.0, 1000.0, 3600.0):
            dTdn = Tw / math.erf(lam) * math.exp(-lam ** 2) / math.sqrt(math.pi * kap * t)      # d/dn with n = -x: gradient of the analytic profile at x = s
            m = interface.melt_speed(dTdn, 0.0, rho_i=RHO_W); self.assertAlmostEqual(m / (lam * math.sqrt(kap / t)), 1.0, 9)
    def test_latent_energy_equals_integrated_conductive_flux(self):
        # independent quadrature of the analytic flux: int_0^t q dt = k Tw/(erf lam sqrt(pi kap)) e^{-lam^2} * 2 sqrt(t)  vs rho L s(t)
        Ste = 0.00631; lam, Tw, _ = cfd_contract.neumann_fixture(Ste, n=2); kap = K_W / (RHO_W * 4217.0); t = 3600.0
        N = 200000; q_int = 0.0
        for i in range(N): u = (i + 0.5) / N; ti = t * u * u; q_int += K_W * Tw / math.erf(lam) * math.exp(-lam ** 2) / math.sqrt(math.pi * kap * ti) * (2 * t * u / N)        # substitution t = t_end u^2 removes the integrable singularity
        self.assertAlmostEqual(q_int / (RHO_W * LF * 2 * lam * math.sqrt(kap * t)), 1.0, 6)
    def test_density_jump_mass_balance(self):
        m = 3.4e-6; vI = -m; uw = interface.water_normal_velocity(m)
        self.assertAlmostEqual(RHO_I * (0.0 - vI) - RHO_W * (uw - vI), 0.0, 15); self.assertAlmostEqual(uw, (RHO_I / RHO_W - 1) * m, 18); self.assertLess(uw, 0)
        self.assertAlmostEqual(interface.water_normal_velocity(m, rho_i=RHO_W), 0.0, 18)      # equal-density limit: liquid at rest for stationary ice
    def test_exchange_record_phase_change_mass_check(self):
        rec = cfd_contract.ExchangeRecord(0.0, 100.0, 50.0, 2.0, 1e-3, RHO_I / RHO_W * 1e-3, 1e5, 9e4, 0.0); rec.validate()
        bad = cfd_contract.ExchangeRecord(0.0, 100.0, 50.0, 2.0, 1e-3, 1e-3, 1e5, 9e4, 0.0)
        with self.assertRaises(ValueError): bad.validate()
    def test_water_energy_residual_zero_for_consistent_step(self):
        q = 200.0; dt = 50.0; rec = cfd_contract.ExchangeRecord(0.0, dt, q, 1.0, q * dt / (RHO_I * LF), RHO_I / RHO_W * q * dt / (RHO_I * LF), 1.0e5, 1.0e5 - q * dt - 300.0, 300.0)
        self.assertAlmostEqual(rec.water_energy_residual(), 0.0, 9); self.assertAlmostEqual(rec.validate()["latent_over_heat"], 0.0, 12)
class TestContractIO(unittest.TestCase):
    def test_csv_roundtrip_and_grain_inputs_match_englacial_source(self):
        import os, tempfile
        from ice_cliff_survival import grains
        rec = cfd_contract.ExchangeRecord(0.0, 3600.0, 120.0, 4.0, 120.0 * 3600 / (RHO_I * LF), RHO_I / RHO_W * 120.0 * 3600 / (RHO_I * LF), 5e5, 4.4e5, 1e3, 0.0)
        p = os.path.join(tempfile.mkdtemp(), "x.csv"); open(p, "w").write("t0,t1,q_wall_W,wetted_area_m2,melt_volume_ice_m3,meltwater_volume_m3,sensible_in_J,sensible_out_J,storage_change_J\n" + ",".join(repr(v) for v in (rec.t0, rec.t1, rec.q_wall_W, rec.wetted_area_m2, rec.melt_volume_ice_m3, rec.meltwater_volume_m3, rec.sensible_in_J, rec.sensible_out_J, rec.storage_change_J)) + "\n")
        r2 = cfd_contract.parse_exchange_csv(p)[0]; self.assertEqual(r2, rec)
        # englacial release from PURE-ICE volume lost vs mixture-sweep bookkeeping: pure-ice volume V_ice corresponds to mixture volume V_ice/(1-c_mix)
        c = 0.06; V_mix = rec.melt_volume_ice_m3 / (1 - c); self.assertAlmostEqual(cfd_contract.volume_to_grain_inputs(rec, c), c * V_mix, 15)
