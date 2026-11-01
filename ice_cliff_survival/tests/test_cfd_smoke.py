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
