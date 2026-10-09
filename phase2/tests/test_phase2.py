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


def test_photosynthetic_pigments_and_empirical_inputs():
    """Chl a/b + carotenoids add absorption within the packaging limit; the empirical inputs
    reproduce the published numbers they are built from."""
    import empirical_data as ED
    # Williamson et al. (2020): 0.04322 ng phenolics per cell; S6 pooled biovolume -> ~22 kg m^-3
    assert abs(ED.pigments_per_cell()["phenolics"][0] - 0.04322) < 1e-4
    assert 21.5 < ED.intracellular_concentration_kg_m3("phenolics") < 22.5
    g = ED.species_geometry()
    for name, V in (("nordenskioeldii", 2307.0), ("alaskanum", 822.0)):   # Halbach et al. (2022)
        d, L = g[name]["diameter_um"], g[name]["length_um"]
        assert abs(np.pi * d ** 2 / 4 * L - V) < 1e-6 * V
    srf = ED.s2_srf_480("S2A")
    wl = np.arange(205, 4996, 10.0)
    for b, centre in (("B2", 492), ("B3", 560), ("B4", 665), ("B8", 833)):   # ESA S2A central wavelengths
        assert abs(np.sum(srf[b] * wl) / srf[b].sum() - centre) < 10
    root = _biosnicar_root()
    if root is None:
        print("BioSNICAR not found - skipping cell part")
        return
    import cell_optics as co
    mac = co.to_480(co.demo_spectrum(root).mac_at)
    kw = co.water_k_480(root)
    geom = P.CellGeometry("cylinder", g["nordenskioeldii"]["diameter_um"] / 2, g["nordenskioeldii"]["length_um"])
    ci = co.empirical_phenolic_concentration()
    base = co.CellModel(geom, ci).optics(mac, kw, packaged=True)
    full = co.CellModel(geom, ci, extra_pigments=co.empirical_pigments()).optics(mac, kw, packaged=True)
    red = (co.WVL_480_NM >= 660) & (co.WVL_480_NM <= 690)                   # chl a red band
    assert np.all(full["abs_xsc"] >= base["abs_xsc"] * (1 - 1e-9))
    assert full["abs_xsc"][red].mean() > 1.5 * base["abs_xsc"][red].mean()
    assert np.all(full["abs_xsc"] <= geom.projected_area_um2 * 1e-12 * (1 + 2e-3))


def test_new_empirical_components():
    """Transmissivity, ice SSA prior, field samples, size scaling, Mie g and the TD-DFT calibration."""
    import empirical_data as ED
    tau, sd, n = ED.clear_sky_transmissivity()
    assert 0.85 < tau < 0.97 and n > 1000                     # PROMICE KAN_M clear-sky hours
    m, sdl, nn = ED.ice_ssa_prior()
    assert nn == 19 and 0.2 < np.exp(m) < 0.6
    r = ED.bubble_radius_um(0.35, 835.0)                        # Cooper et al. (2021) SSA at its density
    assert abs(ED.ssa_from_bubble_radius(r, 835.0) - 0.35) < 1e-12 and 800 < r < 1000
    tab, spectra = ED.field_samples()
    assert (tab.dataset == "s6_2017").sum() == 47 and (tab.dataset == "sgris_2021").sum() == 18
    assert "22_7_SB6" not in set(tab["sample"])                 # not in the primary count workbook
    g = ED.phenolic_size_scaling()
    assert g["n"] == 64 and g["gamma_se"] > 0
    import cell_optics as co
    gg = co.mie_g(7.0, np.zeros(co.WVL_480_NM.size))
    vis = (co.WVL_480_NM >= 400) & (co.WVL_480_NM <= 700)
    assert np.all((gg[vis] > 0.95) & (gg[vis] < 0.999))         # measured microalgae: g > 0.95
    root = _biosnicar_root()
    if root is None:
        print("BioSNICAR not found - skipping calibration part")
        return
    import tddft_calibration as TC
    cal = TC.calibrate(co.demo_spectrum(root), n_steps=800, burn=300, n_stage2=40, draws_per=10, verbose=False)
    d = cal.summary()
    assert d["shape_r2"] > 0.8 and 0 <= d["phi"]["mean"] <= 1
    wl, E, _ = ED.phenolic_extract_mac()
    k = (wl >= 450) & (wl <= 700)
    errC = np.mean(np.log(cal.mac_C(wl[k]) / E[k]) ** 2)
    errD = np.mean(np.log(cal.mac_D(wl[k]) / E[k]) ** 2)
    assert errD <= errC                                          # the Fe term explains the visible tail


def test_stibal2017_manual_ingestion():
    """The manual-ingestion path for Stibal et al. (2017): mapping-driven reading of a counts workbook and
    a spectra table (FIXTURE values, not data), sample matching, solar zenith, and clear errors."""
    import json
    import tempfile
    import pandas as pd
    import stibal2017 as SB
    with tempfile.TemporaryDirectory() as d:
        pd.DataFrame({"ID": ["A1", "A2", "A3"], "cells": [1000.0, 5000.0, 20000.0],
                      "date": ["2014-07-10", "2014-07-20", "2014-08-01"]}).to_excel(
            os.path.join(d, "si.xlsx"), sheet_name="algal cells time series data", index=False)
        wl = np.arange(350, 2501)
        sp = pd.DataFrame({"nm": wl, "A1": 0.6, "A2": 0.4, "A9": 0.3})
        sp.to_csv(os.path.join(d, "spec.csv"), index=False)
        mp = {"counts": {"file": "si.xlsx", "sheet": "algal cells time series data", "sample_col": "ID",
                         "cells_col": "cells", "datetime_col": "date"},
              "spectra": {"file": "spec.csv", "orientation": "wavelength_rows", "wavelength_col": "nm",
                          "quantity": "albedo"},
              "albedo_k_sd": 0.02}
        with open(os.path.join(d, "mapping.json"), "w") as fh:
            json.dump(mp, fh)
        tab, spec, meta = SB.load(d)
        assert list(tab["sample"]) == ["A1", "A2"]                    # A3 has no spectrum, A9 no count
        assert np.allclose(tab.cells, [1000, 5000]) and (tab.quantity == "albedo").all()
        assert np.all((tab.sza > 40) & (tab.sza < 50))               # noon at 67 N in July
        assert meta["albedo_k_sd"] == 0.02 and list(spec.columns) == ["wavelength_nm", "A1", "A2"]
        mp["counts"]["cells_col"] = "wrong"
        with open(os.path.join(d, "mapping.json"), "w") as fh:
            json.dump(mp, fh)
        try:
            SB.load(d)
            raise AssertionError("expected StibalDataError")
        except SB.StibalDataError as e:
            assert "wrong" in str(e) and "cells" in str(e)          # names the missing and available columns


def test_fe_scale_invariance():
    """Assumption 1 of SOURCES.md: if the PG and PG-Fe solutions of Prochazkova et al. (2025) Fig. 4
    did not hold equal amounts, the measured increment D is rescaled by an unknown factor. The fit
    to the extract MAC then rescales phi inversely and the tier D MAC is unchanged."""
    root = _biosnicar_root()
    if root is None:
        print("BioSNICAR not found - skipping")
        return
    import cell_optics as co
    import tddft_calibration as TC
    spec = co.demo_spectrum(root)
    wl_h, S, sS, wl_m, E, sE, D = TC._data()
    t1, _ = TC._fit_magnitude(spec, 0.0, 0.5, wl_m, E, sE, D)
    t2, _ = TC._fit_magnitude(spec, 0.0, 0.5, wl_m, E, sE, 2.0 * D)       # D twice as large
    assert 0 < t1[1] < 1 and abs(t2[1] - t1[1] / 2.0) < 2e-3 * max(1.0, t1[1])
    m1 = np.exp(t1[0]) * (TC.perturbed_mac(spec, 0.0, 1.0, 0.5)(wl_m) + t1[1] * D * 1.0)
    m2 = np.exp(t2[0]) * (TC.perturbed_mac(spec, 0.0, 1.0, 0.5)(wl_m) + t2[1] * 2.0 * D)
    assert np.max(np.abs(m1 / m2 - 1)) < 5e-3


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
    # the bridge takes counts per mL of meltwater; BioSNICAR's own input is per mL of solid ice
    o = run_model(solzen=55, incoming=3, layer_type=[0, 0], rds=[1500, 1500], rho=[650, 850],
                  dz=[0.02, 2.0], glacier_algae=1e4 * bb.MELTWATER_TO_BIOSNICAR)
    assert np.max(np.abs(o.albedo - alb)) < 1e-12 and abs(o.BBA - bba) < 1e-12
    # film geometry: algae uniform over a split crust differs from the 2-layer column only by the
    # delta-Eddington discretisation; with the cells in a thin film they are less shaded by the crust,
    # so the darkening is slightly stronger (a few per cent)
    a0, f, _ = r.run(spec, 55)
    u = bb.IceSpec(1500, 650, 850, 0.02, 2.0, "grains", film_dz=0.002, film_only=False)
    fl = bb.IceSpec(1500, 650, 850, 0.02, 2.0, "grains", film_dz=0.002, film_only=True)
    au, _, _ = r.run(u, 55, r.default_impurity(), 1e4)
    af, _, _ = r.run(fl, 55, r.default_impurity(), 1e4)
    a0u, _, _ = r.run(u, 55)
    du, df_ = r.broadband(a0u, f) - r.broadband(au, f), r.broadband(a0u, f) - r.broadband(af, f)
    assert du > 0 and 1.0 <= df_ / du < 1.12


def test_fast_lut_matches_biosnicar():
    root = _biosnicar_root()
    if root is None:
        print("BioSNICAR not found - skipping")
        return
    import biosnicar_bridge as bb
    from biosnicar.optical_properties.column_OPs import get_layer_OPs
    r = bb.BioSNICARRunner(root)
    for rds in (1000, 1740, 3000):
        ice, ssa, g, mac = r.ice(bb.IceSpec(rds, 650))
        ssa0, g0, mac0 = get_layer_OPs(ice, r.model_config)
        assert np.array_equal(ssa, ssa0) and np.array_equal(g, g0) and np.array_equal(mac, mac0)
    assert r.snap_radius(1509) == 1500 and r.snap_radius(1511) == 1520


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print(f"PASS {name}")


def test_calibration_root_count_sensitivity_flags_truncated_window():
    import cell_optics as co
    import tddft_calibration as tc
    f = np.array([0.1, 0.2, 0.5, 0.4])
    near = co.MolecularSpectrum("x", 300.0, 0.3, np.array([2.5, 3.0, 3.5, 4.6]), f)
    far = co.MolecularSpectrum("x", 300.0, 0.3, np.array([2.5, 3.0, 3.5, 9.0]), f)
    a = tc.root_count_sensitivity(near, 0.0, 0.6, drop=1)
    b = tc.root_count_sensitivity(far, 0.0, 0.6, drop=1)
    assert a["max_rel_change_in_window"] > 0.1 and b["max_rel_change_in_window"] < 1e-6
    assert tc.root_count_sensitivity(near, 0.0, 0.6, drop=4) is None


def test_film_only_places_algae_not_dust_in_film():
    import biosnicar_bridge as bb
    split = bb.IceSpec(1500, 650, 850, 0.02, 2.0, "grains", film_dz=0.002, film_only=True)
    uni = bb.IceSpec(1500, 650, 850, 0.02, 2.0, "grains", film_dz=0.002, film_only=False)
    c_alg = 1e4 * bb.MELTWATER_TO_BIOSNICAR
    assert np.allclose(bb._layer_concs(split, 1e4, 1), [c_alg * 10, 0, 0])
    assert np.allclose(bb._layer_concs(uni, 1e4, 1), [c_alg, c_alg, 0])
    assert bb._layer_concs(split, 5e5, 0) == [5e5, 5e5, 0.0]          # dust: bulk sample, uniform over crust
    assert bb._layer_concs(split, 5e5, 0, in_film=True) == [5e5 * 10, 0.0, 0.0]
    # column amount is conserved by the film rescaling: sum(c_i dz_i)
    dz = [0.002, 0.018]
    for u in (0, 1):
        for sp in (split, uni):
            assert np.isclose(np.dot(bb._layer_concs(sp, 7.0, u)[:2], dz), np.dot(bb._layer_concs(uni, 7.0, u)[:2], dz))
    root = _biosnicar_root()
    if root is None:
        return
    r = bb.BioSNICARRunner(root, incoming=3)
    import copy
    d = np.load(os.path.join(root, "data", "OP_data", "480band", "lap.npz"))
    st = "dust_greenland_Cook_CENTRAL_20190911"
    dust = copy.deepcopy(r.default_impurity())
    dust.name, dust.unit = "dust", 0
    dust.mac, dust.ssa, dust.g = d[st + "__ext_cff_mss"], d[st + "__ss_alb"], d[st + "__asm_prm"]
    # dust-only albedo must not depend on where the ALGAE are assumed to sit
    a1, _, _ = r.run_multi(split, 55, [(dust, 5e5)])
    a2, _, _ = r.run_multi(uni, 55, [(dust, 5e5)])
    assert np.max(np.abs(a1 - a2)) < 1e-12


def test_vacuole_joint_absorption_bounded_and_consistent():
    import pigment_packaging as PP
    for shape, r, L in (("sphere", 6.0, 0.0), ("cylinder", 5.0, 20.0)):
        g = PP.CellGeometry(shape, r, L)
        A = g.projected_area_um2 * 1e-12
        V = g.volume_um3 * 1e-18
        a = np.array([1e2, 1e4, 1e5, 1e7])                       # dilute ... opaque (m^-1)
        # f = 1: joint calculation with all absorption in the "vacuole" equals Q* x a V
        sv, sc = PP.joint_absorption(g, 1.0, a, 0 * a, n=200_000)
        assert np.allclose(sv, PP.q_star(a, g) * a * V, rtol=0.01)
        # dilute limit: no shading, sigma = a_vac V_vac + a_cell V
        sv, sc = PP.joint_absorption(g, 0.3, np.array([1.0]), np.array([2.0]), n=200_000)
        assert np.isclose(sv[0], 0.3 * V, rtol=0.02) and np.isclose(sc[0], 2.0 * V, rtol=0.02)
        # opaque vacuole and opaque cytoplasm: total bounded by the geometric cross-section
        sv, sc = PP.joint_absorption(g, 0.5, np.full(1, 1e7), np.full(1, 1e7), n=200_000)
        assert sv[0] + sc[0] <= A * 1.01
        # the former independent treatment exceeded that bound (regression of the defect)
        old = PP.q_star(np.full(1, 1e7), g.scaled(0.5)) * 1e7 * 0.5 * V + PP.q_star(np.full(1, 1e7), g) * 1e7 * V
        assert old[0] > 1.3 * A
