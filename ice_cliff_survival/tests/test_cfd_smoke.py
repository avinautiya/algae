import math, os, tempfile, unittest
HAVE = os.path.exists("/tmp/claude-0/ofbuild/bin/stefan2DFoam")
@unittest.skipUnless(HAVE, "stefan2DFoam not built (see ice_network_v2/foam/BUILD.md)")
class TestHeatMovingIce2D(unittest.TestCase):
    def test_one_d_limit_ledgers_and_shield_limits(self):
        from ice_cliff_survival.cfd import stefan2d as S
        wd = tempfile.mkdtemp(); r = S.run(wd + "/a", 0.00631, 40, 4, 10.0, 3600.0, jump=True); L = r["L"]; sex = float(S.neumann_s(r, 3600.0)[0])
        self.assertLess(abs(L["s_mean"][-1] - sex) / sex, 3e-5); self.assertEqual(float(L["s_max"][-1] - L["s_min"][-1]), 0.0)                              # 1-D limit: exact Neumann solution, rows identical
        self.assertLess(abs(L["resid_disc"][-1]) / (L["E_hot_cond"][-1] + L["E_hot_adv"][-1]), 1e-12); self.assertLess(abs(L["gcl_resid"][-1]) / L["Vsw_geom"][-1], 1e-12)  # discrete energy ledger and geometric conservation
        self.assertAlmostEqual(L["resid_phys"][-1] / L["E_latent"][-1], 1.66e-3, delta=2e-4)                                                                # sensible+latent residual: first-order in dt (reported, not hidden)
        rs = S.run(wd + "/b", 0.0005, 40, 4, 10.0, 3600.0, jump=False, rho_i=S.RW, hShield=[100.0] * 4); sq = S.shielded_s(rs["Tw"][0], 100.0, rs["params"]["s0"], 3600.0, rho_i=S.RW)
        self.assertLess(abs(rs["L"]["s_mean"][-1] - sq) / sq, 2e-4); self.assertLess(abs(rs["L"]["resid_disc"][-1]) / (rs["L"]["E_hot_cond"][-1] + rs["L"]["E_hot_adv"][-1]), 1e-12)
HAVE2 = os.path.exists("/tmp/claude-0/ofbuild/bin/flowStefanFoam")
@unittest.skipUnless(HAVE2, "flowStefanFoam not built")
class TestFlowStefan(unittest.TestCase):
    def test_small_coupled_case_closes_and_local_mass_holds(self):
        from ice_cliff_survival.cfd import flowstefan as F
        for rho in (1000.0, 917.0):
            R = F.run_case(tempfile.mkdtemp() + "/c", H=0.02, dz=0.01, Lu=0.04, Li=0.08, Ld=0.04, dx=0.008, ny=12, Ubar=0.005, Tin=2.0, dt=10.0, tend=200.0, rho_i=rho); L = R["L"]
            self.assertLess(abs(L["R_disc"][-1]) / L["Ein"][-1], 1e-9); self.assertLess(abs(L["R_phys"][-1]) / L["Eif"][-1], 1e-5); self.assertLess(float(L["massLocal"].max()), 1e-12); self.assertLess(float(abs(L["gcl"]).max()), 1e-15)
            self.assertLess(abs(L["Eif"][-1] - L["Elat"][-1]) / L["Elat"][-1], 1e-5); self.assertLess(float(L["umax"].max()), 0.0151)      # no instability
            if rho != 1000.0: self.assertAlmostEqual(float(L["Qout"][-1] - L["Qin"][-1]), -(1 - rho / 1000.0) * float(L["Vsw"][-1] - L["Vsw"][-2]) / 10.0, delta=1e-13)      # net outflow deficit = (1-r) x melt rate
    def test_unsupported_shield_with_density_jump_refused(self):
        from ice_cliff_survival.cfd.cli_flow import main
        self.assertEqual(main(["--workdir", tempfile.mkdtemp(), "--shield", "--rho-i", "917"]), 2)
class TestReactivationMetrics(unittest.TestCase):
    def test_event_definitions(self):
        import numpy as np
        from ice_cliff_survival.reactivation import metrics as M
        d = np.arange(0, 10.25, 0.25); f = np.where(d < 2, 0.3, 0.0); f = np.where(d >= 8, 0.1, f); e = M.events(d, f, 1e-3, 8.0)
        self.assertAlmostEqual(e["loss"], 1.75 + 0.25 * (0.3 - 1e-3) / 0.3, 6); self.assertEqual(e["sustained_reconnection"], 8.0); self.assertEqual(e["rewet_episodes"], 0)
        f2 = np.where((d > 3) & (d < 3.5), 0.01, np.where(d < 2, 0.3, 0.0)); e2 = M.events(d, f2, 1e-3, None); self.assertEqual(e2["rewet_episodes"], 1); self.assertIsNone(e2["sustained_reconnection"])      # brief rewetting is not reconnection
