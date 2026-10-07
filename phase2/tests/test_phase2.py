"""
Validation tests for Phase 2.   Run:  python -m pytest phase2/tests -q   (or python phase2/tests/test_phase2.py)

1. Monte-Carlo chord packaging reproduces the analytic Duysens / Morel-Bricaud sphere.
2. Q* has the correct dilute and black-body limits; cylinder mean chord obeys Cauchy (4V/S).
3. Independent physics check: exact Mie theory (miepython) for a weakly refracting absorbing
   sphere agrees with Duysens/ADA in shape (<2.5 %) and in magnitude to within the ~n_rel^2
   refractive enhancement ADA neglects (<7 %).
4. Cell optics conserve energy: sigma_abs(packaged) <= projected area, 0 <= SSA <= 1, and
   packaging can only reduce absorption (tier C <= tier B).
5. The in-memory BioSNICAR bridge reproduces BioSNICAR's own run_model() exactly.
"""

import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))

import pigment_packaging as P  # noqa: E402


def test_sphere_mc_matches_analytic():
    g = P.CellGeometry("sphere", 10.0)
    a = np.logspace(2, 7, 21)                         # m^-1
    qa = P.q_star(a, g, "analytic")
    qm = P.q_star(a, g, "chord")
    assert np.max(np.abs(qm / qa - 1.0)) < 3e-3


def test_limits_and_cauchy():
    x = np.array([1e-6, 1e-3, 5e-3])
    assert np.allclose(P.q_star_sphere(x), 1.0 - 3.0 * x / 8.0, atol=1e-5)
    big = np.array([1e3, 1e4])
    assert np.allclose(P.q_star_sphere(big) * big, 1.5, rtol=2e-3)      # black sphere
    # continuity across the series / closed-form switch
    assert abs(P.q_star_sphere(np.array([0.0099999]))[0] - P.q_star_sphere(np.array([0.0100001]))[0]) < 1e-6
    c = P.CellGeometry("cylinder", 5.0, 20.0)
    assert abs(c.chords_um().mean() / c.mean_chord_um - 1.0) < 5e-3
    q = P.q_star(np.array([1e-2, 1e3, 1e5, 1e6]), c)
    assert abs(q[0] - 1.0) < 1e-6 and np.all(np.diff(q) < 0) and np.all(q > 0)


def test_mie_vs_duysens():
    try:
        import miepython
    except ImportError:
        print("miepython not installed - skipping Mie check")
        return
    lam0, nmed, npart, d = 0.45, 1.33, 1.37, 10.0     # um; cell (n=1.37) in meltwater
    x = np.pi * d * nmed / lam0
    rho = np.array([0.05, 0.2, 1.0, 3.0, 10.0])        # a * d
    ratio = []
    for r in rho:
        a_um = r / d
        kp = a_um * lam0 / (4.0 * np.pi)
        m = (npart + 1j * kp) / nmed
        qext, qsca = miepython.efficiencies_mx(m, x)[:2]
        q_ada = P.q_star_sphere(np.array([r]))[0] * 2.0 * r / 3.0
        ratio.append((qext - qsca) / q_ada)
    ratio = np.array(ratio)
    assert np.all((ratio > 1.0) & (ratio < 1.07)), ratio
    assert np.max(np.abs(ratio / ratio[0] - 1.0)) < 0.025, ratio


def _biosnicar_root():
    import biosnicar_bridge as bb
    try:
        return bb.locate_biosnicar()
    except ImportError:
        return None


def test_cell_optics_energy_conservation():
    root = _biosnicar_root()
    if root is None:
        print("BioSNICAR not found - skipping")
        return
    import cell_optics as co
    mac = co.to_480(co.demo_spectrum(root).mac_at)
    kw = co.water_k_480(root)
    for c_int in (5.0, 50.0, 500.0):
        cell = co.CellModel(P.CellGeometry("cylinder", 5.0, 20.0), c_int)
        C = cell.optics(mac, kw, packaged=True)
        B = cell.optics(mac, kw, packaged=False, scatter_from=C)
        A_proj = cell.geom.projected_area_um2 * 1e-12
        assert np.all(C["abs_xsc"] <= A_proj * (1 + 2e-3))
        assert np.all(C["abs_xsc"] <= B["abs_xsc"] * (1 + 1e-12))
        for o in (B, C):
            assert np.all((o["ss_alb"] > 0) & (o["ss_alb"] < 1))
            assert np.allclose(o["ext_xsc"], o["abs_xsc"] + o["sca_xsc"])


def test_bridge_matches_run_model():
    root = _biosnicar_root()
    if root is None:
        print("BioSNICAR not found - skipping")
        return
    import biosnicar_bridge as bb
    from biosnicar.drivers.run_model import run_model
    r = bb.BioSNICARRunner(root, incoming=3)
    spec = bb.IceSpec(1500, 650, 850, 0.02, 2.0, "grains")
    alb, _, bba = r.run(spec, 55, r.default_impurity(), 1e4)
    o = run_model(solzen=55, incoming=3, layer_type=[0, 0], rds=[1500, 1500], rho=[650, 850],
                  dz=[0.02, 2.0], glacier_algae=1e4)
    assert np.max(np.abs(o.albedo - alb)) < 1e-12 and abs(o.BBA - bba) < 1e-12


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print(f"PASS {name}")
