#!/usr/bin/env python
"""
Optics audit: is the low-biomass forcing (2-7x below the published site-differencing estimates) caused by a
bug or an assumption in the molecule -> cell -> BioSNICAR chain? Four numerical checks:

  1. Assay factor. The calibrated strength factor f (= 8.7 x the 1:1 glucoside->phenol value) multiplies the
     TD-DFT MAC per kg glucoside to give MAC per kg PHENOL EQUIVALENT; the cell's pigment load is also in
     phenol equivalents (4-AAP, Williamson et al. 2020). Check: (a) the calibrated tier-D MAC reproduces the
     measured extract MAC (same units) 350-700 nm; (b) the intracellular absorption coefficient
     a_i = MAC x c_i is the same as with the measured extract MAC (the assay factor cancels).
  2. Packaging. Q* (Duysens / Morel-Bricaud, chord-averaged cylinder) is a per-cell property; it must not
     depend on cell density, must be <= 1 and must make sigma_abs <= the geometric projected area.
  3. Units into BioSNICAR. ext_xsc [m^2 cell^-1], unit=1, conc x 0.917 -> column cells m^-2 = count x rho x dz
     (BioSNICAR column_OPs.mix_in_impurities). Check the optical depth BioSNICAR builds against a hand
     calculation, and our ext_xsc magnitude against BioSNICAR's own glacier-algae entries.
  4. Absolute per-cell absorption cross-sections, linear scale, 350-700 nm, against the lab-derived
     entries in BioSNICAR (Chevrollier et al. 2023 in vivo cell absorption = our tier A; Cook et al. 2020
     optics) and the unpackaged "extract MAC x mass per cell" cross-section (Williamson et al. 2020 style).
  5. Consequence: broadband albedo reduction and noon forcing at the published abundance classes for each
     optics, same ice column and dust.

    python phase2/audit_optics.py --outdir records/optics_audit
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path[:0] = [HERE, os.path.join(HERE, "..", "phase4")]

import biosnicar_bridge as bb  # noqa: E402
import cell_optics as co  # noqa: E402
import empirical_data as ED  # noqa: E402

BAND = (350.0, 700.0)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--phase1-l2", default=os.path.join(HERE, "..", "phase1", "results", "level2"))
    p.add_argument("--outdir", default=os.path.join(HERE, "..", "records", "optics_audit"))
    a = p.parse_args(argv)
    os.makedirs(a.outdir, exist_ok=True)
    import emulator as E
    import tddft_calibration as TC
    root = bb.locate_biosnicar(None)
    wl = co.WVL_480_NM
    band = (wl >= BAND[0]) & (wl <= BAND[1])
    out, curves = {}, {"wavelength_nm": wl}

    # ------------------------------------------------------------ 1. assay factor
    cal = TC.cached_calibration(co.load_phase1(a.phase1_l2, "level2"), verbose=False)
    macD = co.to_480(cal.mac_D)
    macC = co.to_480(cal.mac_C)
    macW = ED.pigment_macs_480()["phenolics_williamson2020"]
    c_i = co.empirical_phenolic_concentration()
    r = macD[band] / macW[band]
    out["1_assay"] = dict(
        f=cal.mean("f"), f_stoichiometric=float(cal.spec.molar_mass / TC.PHENOL_MOLAR_MASS),
        c_internal_phenol_eq_kg_m3=float(c_i),
        macD_over_measured_350_700=dict(median=float(np.median(r)), min=float(r.min()), max=float(r.max())),
        macC_over_measured_350_700_median=float(np.median(macC[band] / macW[band])),
        a_i_tierD_over_a_i_measured_median=float(np.median(macD[band] * c_i / (macW[band] * c_i))),
        a_i_at_400_550_650_per_m=[float(np.interp(x, wl, macD * c_i)) for x in (400, 550, 650)])
    curves.update(mac_tierD=macD, mac_tierC=macC, mac_measured_extract=macW)

    # ------------------------------------------------------------ build cells exactly as the emulator does
    xs = {}
    for tier, phenol in (("D", "tddft"), ("C", "tddft"), ("C", "williamson2020")):
        cfg = E.EmulatorConfig(model="ours", tier=tier, phenol=phenol)
        b = E._Builder(cfg, a.phase1_l2, False, None)
        for name, imp in b.imps.items():
            xs[(tier, phenol, name)] = imp
    species = E.empirical_species()
    from pigment_packaging import CellGeometry
    rows = []
    for name, sp in species.items():
        g = CellGeometry("cylinder", sp.diameter_um / 2.0, sp.length_um)
        cell = co.CellModel(g, co.empirical_phenolic_concentration(name), vd_diagnostic=False,
                            extra_pigments=co.empirical_pigments(species=name))
        kw = co.water_k_480(root)
        oD = cell.optics(macD, kw, packaged=True)
        oB = cell.optics(macD, kw, packaged=False, scatter_from=oD)
        A = g.projected_area_um2 * 1e-12
        imp = xs[("D", "tddft", name)]
        same = float(np.max(np.abs(imp.mac - oD["ext_xsc"]) / oD["ext_xsc"]))
        q = oD["q_star"]
        rows.append(dict(species=name, volume_um3=g.volume_um3, projected_area_m2=A,
                         q_star_350_700_mean=float(q[band].mean()), q_star_min=float(q.min()), q_star_max=float(q.max()),
                         sigma_abs_packaged_over_area_350_700=float((oD["abs_xsc"][band] / A).mean()),
                         sigma_abs_unpackaged_over_area_350_700=float((oB["abs_xsc"][band] / A).mean()),
                         emulator_matches_rebuild_max_rel_diff=same))
        curves[f"abs_xsc_tierD_{name}"] = oD["abs_xsc"]
        curves[f"abs_xsc_unpackaged_{name}"] = oB["abs_xsc"]
        curves[f"abs_xsc_measuredMAC_{name}"] = xs[("C", "williamson2020", name)].mac * (
            1 - xs[("C", "williamson2020", name)].ssa)
        curves[f"area_{name}"] = np.full(wl.size, A)
    out["2_packaging"] = rows
    # density independence of Q*: Q* is computed per cell; check by construction (no conc argument) and
    # numerically: the same CustomImpurity is used at every abundance (emulator.node passes B to conc only)
    out["2_packaging_density_dependence"] = "none: CellModel.optics has no concentration argument; " \
        "abundance enters only BioSNICAR's conc (emulator.node)"

    # ------------------------------------------------------------ 3. units into BioSNICAR
    run = bb.BioSNICARRunner(root, incoming=3)
    spec = bb.IceSpec(600.0, 450.0, rho_bottom=690.0, mode="bubbly")
    imp = xs[("D", "tddft", "nordenskioeldii")]
    conc = 1e4
    cells_m2_hand = conc * 0.917 / 917 * 1e6 * 450.0 * 0.02           # = conc[1/mL meltwater]*rho*dz/1000
    tau_hand_500 = cells_m2_hand * float(np.interp(500, wl, imp.mac))
    out["3_units"] = dict(cells_per_m2_at_1e4_per_mL=cells_m2_hand,
                          expected_from_count_x_rho_x_dz=conc * 450.0e3 * 0.02,   # cells/g x g/m^3 x m
                          tau_algae_500nm_at_1e4=tau_hand_500,
                          ext_xsc_500nm_ours_m2=float(np.interp(500, wl, imp.mac)),
                          area_fraction_covered_at_1e4=float(cells_m2_hand * rows[0]["projected_area_m2"]))
    lap = np.load(os.path.join(root, "data", "OP_data", "480band", "lap.npz"))
    refs = {}
    for stem in ("ice_algae_empirical_Chevrollier2023", "Cook2020_glacier_algae_4_40", "Glacier_Algae_IS"):
        if f"{stem}__ext_xsc" in lap.files:
            e, s = lap[f"{stem}__ext_xsc"], lap[f"{stem}__ss_alb"]
            refs[stem] = (e, s)
            curves[f"abs_xsc_{stem}"] = e * (1 - s)
    out["3_units_reference_ext_xsc_500nm"] = {k: float(np.interp(500, wl, v[0])) for k, v in refs.items()}

    # ------------------------------------------------------------ 4. magnitude comparison (linear, 350-700)
    comm = {k: 0.6 * curves[f"{k}_nordenskioeldii"] + 0.4 * curves[f"{k}_alaskanum"]
            for k in ("abs_xsc_tierD", "abs_xsc_unpackaged", "abs_xsc_measuredMAC")}
    mag = {}
    for k, v in comm.items():
        mag[k] = float(v[band].mean())
    for stem, (e, s) in refs.items():
        mag[f"abs_xsc_{stem}"] = float((e * (1 - s))[band].mean())
    ref = mag.get("abs_xsc_ice_algae_empirical_Chevrollier2023")
    out["4_band_mean_abs_xsc_m2_350_700"] = mag
    out["4_ratio_to_Chevrollier2023_measured"] = {k: v / ref for k, v in mag.items()} if ref else None

    # ------------------------------------------------------------ 5. consequence at published abundances
    dust = [(b.dust, float(np.exp(ED.dust_prior()[0])))]
    a0, flx, _ = run.run_multi(spec, 47.0, dust)
    sw = bb.sw_down_clear_sky(47.0, ED.clear_sky_transmissivity()[0], 202)
    cons = []
    optics = {"tier D (ours)": [(xs[("D", "tddft", "nordenskioeldii")], 0.6), (xs[("D", "tddft", "alaskanum")], 0.4)],
              "measured extract MAC, packaged": [(xs[("C", "williamson2020", "nordenskioeldii")], 0.6),
                                                 (xs[("C", "williamson2020", "alaskanum")], 0.4)]}
    for stem, (e, s) in refs.items():
        optics[stem] = [(bb.CustomImpurity(stem, e, s, lap[f"{stem}__asm_prm"]), 1.0)]
    for label, conc_ in (("Cook Hbio", 2.9e4), ("Cook Lbio", 4.73e3), ("Williamson high", 8989.0)):
        for oname, parts in optics.items():
            alb, _, _ = run.run_multi(spec, 47.0, [(im, conc_ * w) for im, w in parts] + dust)
            cons.append(dict(case=label, cells_per_mL=conc_, optics=oname,
                             d_bba=float(run.broadband(a0, flx) - run.broadband(alb, flx)),
                             rf_noon_W_m2=float(run.forcing(a0, alb, flx, sw))))
    out["5_consequence"] = cons
    out["5_published_site_albedo_difference"] = dict(
        note="Cook et al. 2020 Table 2 (UAV): Hbio 0.25, Lbio 0.44, clean ice 0.53 -> d_bba 0.28 (Hbio), 0.09 (Lbio)")

    with open(os.path.join(a.outdir, "optics_audit.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    pd.DataFrame(curves).to_csv(os.path.join(a.outdir, "optics_audit_curves.csv"), index=False, float_format="%.5g")
    pd.DataFrame(cons).to_csv(os.path.join(a.outdir, "optics_audit_consequence.csv"), index=False, float_format="%.4g")
    figure(curves, comm, refs, wl, a.outdir)
    print(json.dumps({k: v for k, v in out.items() if k != "5_consequence"}, indent=1, default=float))
    print(pd.DataFrame(cons).round(4).to_string(index=False))


def figure(curves, comm, refs, wl, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.6))
    m = (wl >= 330) & (wl <= 750)
    sc = 1e10
    ax[0].plot(wl[m], comm["abs_xsc_tierD"][m] * sc, color="#1b6ca8", lw=1.6, label="ours, tier D (packaged)")
    ax[0].plot(wl[m], comm["abs_xsc_measuredMAC"][m] * sc, color="#e07b39", lw=1.2, label="measured extract MAC, packaged")
    unp = (comm["abs_xsc_unpackaged"][m] / (0.6 * curves["area_nordenskioeldii"][m] + 0.4 * curves["area_alaskanum"][m]))
    ax[0].text(0.03, 0.97, f"unpackaged (dissolved) pigment: {unp.min():.0f}-{unp.max():.0f}x the\nprojected area (off scale, unphysical)",
               transform=ax[0].transAxes, fontsize=5.5, va="top")
    cols = {"ice_algae_empirical_Chevrollier2023": "k", "Cook2020_glacier_algae_4_40": "#7a5195", "Glacier_Algae_IS": "#999999"}
    for stem, (e, s) in refs.items():
        ax[0].plot(wl[m], (e * (1 - s))[m] * sc, color=cols.get(stem, "0.5"), ls="--", lw=1.2, label=stem.replace("_", " "))
    A = 0.6 * curves["area_nordenskioeldii"] + 0.4 * curves["area_alaskanum"]
    ax[0].plot(wl[m], A[m] * sc, color="0.6", lw=0.8, label="projected area (black cell)")
    ax[0].set_ylabel(r"$\sigma_{abs}$ ($10^{-10}$ m$^2$ cell$^{-1}$)")
    ax[0].set_xlabel("wavelength (nm)")
    ax[0].set_ylim(0, None)
    ax[0].legend(fontsize=5.5, frameon=False, loc="lower left")
    ax[0].set_title("a  per-cell absorption, linear", loc="left", fontsize=8)
    ref = "ice_algae_empirical_Chevrollier2023"
    if ref in refs:
        e, s_ = refs[ref]
        ax[1].plot(wl[m], (comm["abs_xsc_tierD"] / (e * (1 - s_)))[m], color="#1b6ca8", label="ours tier D")
        ax[1].plot(wl[m], (comm["abs_xsc_measuredMAC"] / (e * (1 - s_)))[m], color="#e07b39", label="measured MAC")
        ax[1].axhline(1, color="0.5", lw=0.6)
        ax[1].set_ylim(0, 2)
        ax[1].set_ylabel("ratio to measured cells")
        ax[1].set_xlabel("wavelength (nm)")
        ax[1].legend(fontsize=6, frameon=False)
        ax[1].set_title("b  vs Chevrollier 2023", loc="left", fontsize=8)
    ax[2].semilogy(wl[m], curves["mac_measured_extract"][m], color="#e07b39", label="measured extract MAC")
    ax[2].semilogy(wl[m], curves["mac_tierD"][m], color="#1b6ca8", label="tier D (TD-DFT + Fe)")
    ax[2].semilogy(wl[m], curves["mac_tierC"][m], color="#1b6ca8", ls="--", label="tier C")
    ax[2].set_ylabel(r"MAC (m$^2$ kg$^{-1}$ phenol eq.)")
    ax[2].set_xlabel("wavelength (nm)")
    ax[2].legend(fontsize=6, frameon=False)
    ax[2].set_title("c  pigment MAC", loc="left", fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "FigS8_optics_audit.png"), dpi=600)
    fig.savefig(os.path.join(outdir, "FigS8_optics_audit.pdf"))


if __name__ == "__main__":
    main()
