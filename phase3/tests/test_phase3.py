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
    import tddft_calibration as TC
    fm = ForwardModel(demo=True)                     # defaults: bubbly ice, chl/carotenoids, calibrated pigment
    V, ar = 2000.0, 2.0
    c = fm.cal
    p = dict(dE_ev=c["dE"], f_scale=c["f"], fwhm_ev=c["w"], cell_volume_um3=V, cell_aspect=ar,
             c_internal=22.0, conc_size_exponent=0.0, ice_ssa=0.38, rho_top=450.0, conc_cells_ml=1e4,
             fe_fraction=c["phi"], transmissivity=0.92)
    out = fm.evaluate(p)
    d = (4 * V / (np.pi * ar)) ** (1 / 3)
    assert np.isclose(out["cell_diameter_um"], d) and np.isclose(out["cell_length_um"], ar * d)
    r_b = 3 * (1 - 690 / 917) / (690 * 0.38) * 1e6                 # SSA -> BioSNICAR bubble radius
    assert np.isclose(out["ice_radius_um"], r_b)
    # Phase 2 reference: direct chord Q*, fresh runner, same calibrated MAC
    root = bb.locate_biosnicar()
    lig = co.demo_spectrum(root)
    mac = co.to_480(TC.perturbed_mac(lig, c["dE"], c["f"], c["w"]))
    cell = co.CellModel(CellGeometry("cylinder", d / 2, ar * d), 22.0, extra_pigments=co.empirical_pigments())
    oC = cell.optics(mac, co.water_k_480(root), packaged=True)
    r = bb.BioSNICARRunner(root)
    spec = bb.IceSpec(r_b, 450, rho_bottom=690.0, mode="bubbly")
    a0, flx, _ = r.run(spec, 45)
    a1, _, _ = r.run(spec, 45, bb.CustomImpurity("C", oC["ext_xsc"], oC["ss_alb"], oC["asm_prm"]), 1e4)
    rf = r.forcing(a0, a1, flx, bb.sw_down_clear_sky(45, 0.92))
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


def test_joint_calibration_draw_keeps_correlations_and_replicate_ci():
    import stats_tools as st
    from param_space import default_parameters
    cal = {k: dict(mean=m, sd=s) for k, m, s in (("dE", 0.04, 0.01), ("w", 0.6, 0.01), ("f", 60.0, 20.0),
                                                 ("phi", 0.12, 0.03))}
    cal["log_f"] = dict(mean=np.log(60.0), sd=0.3)
    pj = default_parameters(True, cal, joint_calibration=True)
    names = [p.name for p in pj]
    assert "cal_draw" in names and not {"dE_ev", "f_scale", "fwhm_ev", "fe_fraction"} & set(names)
    pi = [p.name for p in default_parameters(True, cal, joint_calibration=False)]
    assert {"dE_ev", "f_scale", "fwhm_ev", "fe_fraction"} <= set(pi)
    # the forward model maps the draw onto ONE joint posterior row
    from forward_model import ForwardModel
    rng = np.random.default_rng(0)
    z = rng.normal(size=1000)
    samples = np.column_stack([0.04 + 0.01 * rng.normal(size=1000), 0.6 + 0.01 * rng.normal(size=1000),
                               60 + 20 * z, 0.12 - 0.03 * z])                    # f and phi anti-correlated
    fm = ForwardModel.__new__(ForwardModel)
    fm.cal, fm.cal_samples = dict(dE=0, w=0.6, f=1, phi=0), samples
    seen = []
    fm.co = type("co", (), {"to_480": staticmethod(lambda f, *a: f)})
    fm.TC = type("tc", (), {"perturbed_mac": staticmethod(lambda lig, dE, f, w: (dE, w, f)),
                            "complexed_mac": staticmethod(lambda lig, dE, f, w, phi: seen.append((f, phi)) or 0)})
    fm.ligand, fm.window, fm.uv_mode = None, None, None
    for u in rng.random(300):
        fm.mac480(dict(cal_draw=u))
    f_, ph = np.array(seen).T
    assert np.corrcoef(f_, ph)[0, 1] < -0.9
    # LHS replicate CI covers the true mean, and is computed from replicate spread
    x = rng.normal(5.0, 1.0, 2000)
    lo, hi = st.replicate_ci(x, np.repeat(np.arange(10), 200))
    assert lo < 5.0 < hi and hi - lo < 0.3
    d = st.describe(x, replicate=np.repeat(np.arange(10), 200))
    assert d["ci_method"].startswith("10 LHS")
