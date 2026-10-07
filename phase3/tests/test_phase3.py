"""
Validation tests for Phase 3.   Run:  python phase3/tests/test_phase3.py   (or pytest)

1. Sobol' machinery (ParameterSpace transform + Saltelli design + SALib analysis)
   recovers the analytic indices of the Ishigami benchmark function.
2. Latin Hypercube sample is exactly stratified in every dimension.
3. Q* look-up table agrees with direct Monte Carlo chord averaging (<0.5 %).
4. The fast forward model reproduces the Phase 2 physics (direct Q*, uncached ice optics).
5. Bootstrap CIs of the mean cover the true mean for a known distribution.
"""

import os
import sys

import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))
sys.path.insert(0, os.path.join(HERE, "..", "..", "phase2"))

from param_space import Param, ParameterSpace  # noqa: E402
import stats_tools as st  # noqa: E402


def test_ishigami_sobol():
    a, b = 7.0, 0.1
    U = stats.uniform(-np.pi, 2 * np.pi)
    space = ParameterSpace([Param(f"x{i}", f"x{i}", "-", "molecular", U, "") for i in range(1, 4)])
    X, prob = space.saltelli(4096, second_order=True, seed=1)
    Y = np.sin(X[:, 0]) + a * np.sin(X[:, 1]) ** 2 + b * X[:, 2] ** 4 * np.sin(X[:, 0])
    tab, s2 = st.sobol(prob, Y, second_order=True, n_resamples=100)
    # analytic values (Ishigami & Homma 1990; Sobol' & Levitan 1999)
    V = a**2 / 8 + b * np.pi**4 / 5 + b**2 * np.pi**8 / 18 + 0.5
    V1 = 0.5 * (1 + b * np.pi**4 / 5) ** 2
    V2 = a**2 / 8
    V13 = b**2 * np.pi**8 / 18 - b**2 * np.pi**8 / 50
    S1 = np.array([V1, V2, 0.0]) / V
    ST = np.array([V1 + V13, V2, V13]) / V
    assert np.allclose(tab.S1.to_numpy(), S1, atol=0.03), tab
    assert np.allclose(tab.ST.to_numpy(), ST, atol=0.03), tab
    assert abs(s2.loc["x1", "x3"] - V13 / V) < 0.05


def test_lhs_stratified():
    space = ParameterSpace([Param(f"x{i}", "", "-", "molecular", stats.uniform(0, 1), "") for i in range(5)])
    X = space.lhs(200, seed=3)
    for j in range(5):
        counts = np.bincount(np.floor(X[:, j] * 200).astype(int), minlength=200)
        assert np.all(counts == 1)


def test_qstar_table():
    import pigment_packaging as P
    from qstar_table import QStarTable
    T = QStarTable(n_aspect=21, n_tau=241)
    rng = np.random.default_rng(5)
    for k in range(6):
        g = P.CellGeometry("cylinder", rng.uniform(2.5, 7.5), rng.uniform(10, 30), n_chords=200_000, seed=50 + k)
        a = np.logspace(2, 7, 40)
        assert np.max(np.abs(T(a, g) / P.q_star(a, g) - 1)) < 5e-3


def _root():
    import biosnicar_bridge as bb
    try:
        return bb.locate_biosnicar()
    except ImportError:
        return None


def test_forward_model_matches_phase2():
    root = _root()
    if root is None:
        print("BioSNICAR not found - skipping")
        return
    import biosnicar_bridge as bb
    import cell_optics as co
    from pigment_packaging import CellGeometry
    from forward_model import ForwardModel
    fm = ForwardModel(demo=True)
    p = dict(dE_ev=0.0, f_scale=1.0, fwhm_ev=0.30, cell_length_um=20.0, cell_diameter_um=10.0,
             c_internal=50.0, grain_um=1500.0, rho_top=650.0, conc_cells_ml=1e4,
             lmct_eps=4000.0, lmct_center_nm=570.0)
    out = fm.evaluate(p)
    # Phase 2 reference: direct chord Q*, fresh runner
    mac = co.to_480(co.demo_spectrum(root).mac_at)
    cell = co.CellModel(CellGeometry("cylinder", 5.0, 20.0), 50.0)
    oC = cell.optics(mac, co.water_k_480(root), packaged=True)
    r = bb.BioSNICARRunner(root)
    spec = bb.IceSpec(1500, 650)
    a0, flx, _ = r.run(spec, 55)
    a1, _, _ = r.run(spec, 55, bb.CustomImpurity("C", oC["ext_xsc"], oC["ss_alb"], oC["asm_prm"]), 1e4)
    rf = r.forcing(a0, a1, flx, bb.sw_down_clear_sky(55))
    assert abs(out["rf_C"] / rf - 1) < 5e-3, (out["rf_C"], rf)


def test_bootstrap_coverage():
    rng = np.random.default_rng(11)
    hits = 0
    for k in range(60):
        x = rng.lognormal(0.0, 1.0, size=400)
        d = st.describe(x, n_boot=500, seed=k)
        hits += d["mean_ci_lo"] <= np.exp(0.5) <= d["mean_ci_hi"]
    assert hits >= 50          # nominal 95 % -> expect ~57/60; percentile bootstrap undercovers slightly


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print(f"PASS {name}")
