import datetime as dt, math, os, unittest
from ice_cliff_survival import transport, observation, forcing, UnsupportedRun
from ice_cliff_survival.config import ToeConfig, Provenance
from ice_cliff_survival.runner import ToeRun
HAVE = os.path.exists(os.path.join(forcing.RAW, forcing.KENN)) and os.path.exists(os.path.join(forcing.RAW, forcing.STREAMS))
class TestTransport(unittest.TestCase):
    def test_reproduces_published_petersen_clast_sizes(self):
        # Petersen et al. 2024 Sect. 5.2: v = 0.6-1.13 m/s mobilises half-spheres 9-30 cm and cubes 14-51 cm (rho 2700, mu 0.05, no buoyancy)
        self.assertAlmostEqual(100 * transport.max_mobile_diameter(0.6, "half_sphere"), 8.6, delta=0.4); self.assertAlmostEqual(100 * transport.max_mobile_diameter(1.13, "half_sphere"), 30.4, delta=0.7)
        self.assertAlmostEqual(100 * transport.max_mobile_diameter(0.6, "cube"), 14.3, delta=0.5); self.assertAlmostEqual(100 * transport.max_mobile_diameter(1.13, "cube"), 50.6, delta=0.8)
    def test_scaling_and_buoyancy_direction(self):
        self.assertAlmostEqual(transport.max_mobile_diameter(2.0) / transport.max_mobile_diameter(1.0), 4.0, 12)
        self.assertGreater(transport.max_mobile_diameter(1.0, buoyancy=True), transport.max_mobile_diameter(1.0, buoyancy=False))      # buoyancy lowers the normal force on the ice, so larger clasts are mobilised
    def test_clearing_capacity_zero_below_threshold(self):
        c = transport.ClearingModel(0.3, 0.8, 0.3, 5.0); self.assertEqual(c.capacity(10.0, 0.7), 0.0); self.assertAlmostEqual(c.capacity(10.0, 1.0), 10.0 * 0.3 * 1.0 / 5.0, 12)
    def test_partial_submergence_factor(self):
        self.assertEqual(transport.partial_submergence_factor(0.2, 0.5), 1.0); self.assertAlmostEqual(transport.partial_submergence_factor(0.4, 0.1), 0.25)
class TestObservation(unittest.TestCase):
    def test_horizontal_to_normal_recession_planar_face(self):
        for th in (30.0, 55.0, 90.0):
            M = 4e-6; t = math.radians(th); n = (math.sin(t), -math.cos(t))            # unit normal of a face with direction (cos th, sin th); stake along +x
            def D_of(time):                                                              # face line n.(x,y) = c0 + M time, stake on y = 0: x = (c0 + M time)/n_x
                return (0.0 + M * time) / n[0]
            dD = D_of(86400.0) - D_of(0.0); self.assertAlmostEqual(observation.normal_recession_from_horizontal(dD / 86400.0, th) / M, 1.0, 12); self.assertAlmostEqual(observation.horizontal_from_normal(M, th) * 86400.0, dD, 9)
    def test_interval_rates_irregular_and_uncertainty(self):
        t = [dt.datetime(2022, 6, 1, 0), dt.datetime(2022, 6, 2, 6), dt.datetime(2022, 6, 5, 6)]; D = [10.0, 40.0, 70.0]; r = observation.interval_rates(t, D)
        self.assertEqual(len(r), 2); self.assertAlmostEqual(r[0]["rate_cm_day"], 30.0 / 1.25, 12); self.assertAlmostEqual(r[1]["rate_cm_day"], 10.0, 12); self.assertGreater(r[0]["sigma_cm_day"], r[1]["sigma_cm_day"])
    def test_bad_slope_rejected(self):
        with self.assertRaises(ValueError): observation.normal_recession_from_horizontal(1.0, 0.0)
    @unittest.skipUnless(HAVE, "pinned archives not downloaded")
    def test_real_backwasting_loader(self):
        L = forcing.load_backwasting(); self.assertEqual(len(L), 13); self.assertEqual(sum(len(x["D_cm"]) for x in L), 304)
        for x in L: self.assertTrue(all(b >= a for a, b in zip(x["D_cm"], x["D_cm"][1:])))
    @unittest.skipUnless(HAVE, "pinned archives not downloaded")
    def test_gauge_discharge_in_published_range(self):
        import zipfile, csv, io
        z = forcing.open_zip(forcing.STREAMS); N = [n for n in z.namelist() if "gauge_1_2021_07_13/" in n and "__MACOSX" not in n]
        profs = [[tuple(map(float, l.split(","))) for l in z.read(n).decode("utf-8-sig").strip().splitlines()] for n in N if "prof" in n and n.endswith(".csv")]
        tm = [r for r in csv.DictReader(io.StringIO(z.read([n for n in N if n.endswith("time.csv")][0]).decode("utf-8-sig"))) if r["Use"].strip() == "1"]; T = [float(r["Time"]) for r in tm]
        q = forcing.gauge_discharge(profs, 8.0, T); self.assertTrue(0.06 <= q["Q_m3_s"] <= 0.21)
class TestRunner(unittest.TestCase):
    def mk(self):
        cfg = ToeConfig(W=10.0, porosity=0.3, H_cliff=15.0, theta_deg=60.0, beta_deg=35.0); P = Provenance("c", "synthetic")
        return ToeRun(cfg, lambda t: {"S_mantle": 3e-5}, lambda t, V: 1e-5, dict(sources=P, capacity=P), "const", 600.0)
    def test_deterministic_restart_equivalence(self):
        a = self.mk(); a.advance(60 * 600.0); b = self.mk(); b.advance(30 * 600.0); p = "/tmp/_icsv_ck.json"; b.save(p); c = self.mk(); c.load(p); c.advance(60 * 600.0)
        self.assertEqual(a.stock.Vs, c.stock.Vs); self.assertEqual(a.stock.cum, c.stock.cum)
    def test_restart_refused_on_mismatch(self):
        b = self.mk(); b.advance(1200.0); p = "/tmp/_icsv_ck2.json"; b.save(p); cfg = ToeConfig(W=11.0, porosity=0.3, H_cliff=15.0, theta_deg=60.0, beta_deg=35.0); P = Provenance("c", "synthetic")
        c = ToeRun(cfg, lambda t: {"S_mantle": 3e-5}, lambda t, V: 1e-5, dict(sources=P, capacity=P), "const", 600.0)
        with self.assertRaises(ValueError): c.load(p)
    def test_environmental_run_refused(self):
        cfg = ToeConfig(W=10.0, porosity=0.3, H_cliff=15.0, theta_deg=60.0, beta_deg=35.0); P = Provenance("c", "measured")
        with self.assertRaises(UnsupportedRun): ToeRun(cfg, lambda t: {}, lambda t, V: 0.0, dict(sources=P, capacity=P), "x", 600.0, label="environmental")
    def test_config_validation(self):
        with self.assertRaises(ValueError): ToeConfig(W=10.0, porosity=1.0, H_cliff=15.0, theta_deg=60.0, beta_deg=35.0)
        with self.assertRaises(ValueError): ToeConfig(W=10.0, porosity=0.3, H_cliff=15.0, theta_deg=30.0, beta_deg=35.0)
