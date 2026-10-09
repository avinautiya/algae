#!/usr/bin/env python
"""
Molecular-contribution study (docs/molecular_contribution_protocol.md). Outputs: records/molecular_contribution/.

Stages:
  raw         Q1: uncalibrated TD-B3LYP and TD-CAM-B3LYP spectra vs measured spectra. Band maxima, unit-area
              shape residual vs the isolated HPLC chromophore (dE = 0, FWHM 0.3 eV), visible share of the
              absorption vs the whole extract. The purpurogallin (aglycone) absorbance of Prochazkova 2025
              is used as a related-chromophore position check only.
  calibrate   identical calibration (AR(1), same bounds and data) for both functionals; cached on disk by
              the calibration module's own content fingerprint.
  compare     posterior parameters, correlations and calibrated-MAC differences between functionals,
              tier C (uncomplexed) and D (with the measured Fe increment).
  expA        Experiment A: identical cell model, packaging, photosynthetic background, column, dust,
              illumination and observation quantity; ONLY the phenolic absorption spectrum changes
              (T-B3-C/D, T-CAM-C/D, T-MEAS). Hemispherical spectral albedo, S2 bands, broadband 300-2500 nm
              and absorbed SW at SZA 50 deg clear sky, on a state grid plus zero algae. Joint posterior
              draws propagate (same draw through MAC -> cell -> albedo).

    python common/run_budgeted.py --background --mem-mb 900 --name molecular -- \
        python3 phase4/molecular_contribution.py --stages raw calibrate compare expA
"""
from __future__ import annotations

import argparse
import copy
import json
import os
import pickle
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, "..")
sys.path[:0] = [HERE, os.path.join(ROOT, "phase2"), os.path.join(ROOT, "phase1")]

import empirical_data as ED  # noqa: E402

OUT = os.path.join(ROOT, "records", "molecular_contribution")
SPECS = {"T-B3": (os.path.join(ROOT, "phase1", "results", "level2"), "B3LYP"),
         "T-CAM": (os.path.join(ROOT, "phase1", "results", "legacy_2026-10-09", "level2_cam"), "CAM-B3LYP")}
CAL_CACHE = os.path.join(HERE, "results", "cache", "calibrations")
WL = np.arange(250.0, 750.5, 1.0)
VIS = (WL >= 400) & (WL <= 700)
N_DRAWS = 24
SZA = 50.0
STATES = [dict(log_b=lb, r_um=r, dust_ppb=d) for lb in (3.0, 4.0, 4.5, 5.0) for r in (1000.0, 3000.0, 8000.0)
          for d in (3e4, 3e5)]


def load_specs():
    import cell_optics as co
    return {k: co.load_phase1(d, "level2", fn) for k, (d, fn) in SPECS.items()}


def local_maxima(x, y, lo=260, hi=700):
    k = (x >= lo) & (x <= hi)
    xs, ys = x[k], y[k]
    i = np.flatnonzero((ys[1:-1] > ys[:-2]) & (ys[1:-1] >= ys[2:])) + 1
    return [round(float(xs[j]), 1) for j in i if ys[j] > 0.05 * ys.max()]


def stage_raw(specs):
    import tddft_calibration as TC
    wl_h, S, sS = ED.chromophore_hplc_shape()
    wl_m, E, sE = ED.phenolic_extract_mac()
    wl_pg, A_pg, A_pgfe = ED.fe_purpurogallin_absorbance()
    rows = []
    hplc_max = local_maxima(wl_h, S, 265, 600)
    for k, sp in specs.items():
        M = TC.perturbed_mac(sp, 0.0, 1.0, 0.3)
        mh = M(wl_h)
        shape = mh / np.trapezoid(mh, wl_h)
        m = M(WL)
        sticks = sorted(zip(sp.energies_ev, sp.osc))
        rows.append(dict(treatment=k, n_roots=int(sp.energies_ev.size), molar_mass=sp.molar_mass,
                         S1_nm=round(1239.84 / sticks[0][0], 1), S1_f=round(float(sticks[0][1]), 4),
                         raw_maxima_nm=local_maxima(WL, m), hplc_maxima_nm=hplc_max,
                         shape_rmse_vs_hplc=float(np.sqrt(np.mean((shape - S) ** 2))),
                         shape_rmse_rel=float(np.sqrt(np.mean((shape - S) ** 2)) / S.mean()),
                         raw_visible_fraction=float(np.trapezoid(m[VIS], WL[VIS]) / np.trapezoid(m[WL >= 265], WL[WL >= 265])),
                         raw_mac_peak_m2_kg=float(m.max())))
    vis_e = (wl_m >= 400) & (wl_m <= 700)
    meas = dict(extract_visible_fraction=float(np.trapezoid(E[vis_e], wl_m[vis_e]) / np.trapezoid(E[wl_m >= 265], wl_m[wl_m >= 265])),
                hplc_maxima_nm=hplc_max, hplc_visible_fraction=float(np.trapezoid(S[(wl_h >= 400)], wl_h[wl_h >= 400])),
                purpurogallin_maxima_nm=local_maxima(wl_pg, A_pg, 260, 700),
                purpurogallin_Fe_maxima_nm=local_maxima(wl_pg, A_pgfe, 260, 700),
                note="HPLC shapes and extract MAC were used to calibrate (not independent); purpurogallin is the "
                     "aglycone without the carboxylic acid/glucoside: related chromophore, position check only")
    return pd.DataFrame(rows), meas


def calibration(spec):
    import tddft_calibration as TC
    os.makedirs(CAL_CACHE, exist_ok=True)
    for f in os.listdir(CAL_CACHE):
        if f.endswith(".pkl") and f[:-4] not in TC._CACHE:
            try:
                TC._CACHE[f[:-4]] = pickle.load(open(os.path.join(CAL_CACHE, f), "rb"))
            except Exception:  # noqa: BLE001
                pass
    cal = TC.cached_calibration(spec, verbose=False, residuals="ar1")
    for k, v in TC._CACHE.items():
        dst = os.path.join(CAL_CACHE, f"{k}.pkl")
        if not os.path.exists(dst):
            pickle.dump(v, open(dst + ".tmp", "wb"))
            os.replace(dst + ".tmp", dst)
    return cal


def draw_macs(cal, n, seed=0):
    import tddft_calibration as TC
    rng = np.random.default_rng(seed)
    idx = rng.choice(cal.samples.shape[0], n, replace=False)
    C, D = [], []
    for i in idx:
        dE, w, f, phi = cal.samples[i, :4]
        C.append(TC.perturbed_mac(cal.spec, dE, f, w)(WL))
        D.append(TC.complexed_mac(cal.spec, dE, f, w, phi)(WL))
    return np.array(C), np.array(D), idx


def stage_compare(cals):
    rows, macs = [], {}
    wl_m, E, sE = ED.phenolic_extract_mac()
    meas = np.interp(WL, wl_m, E, left=np.nan, right=np.nan)
    import biosnicar_bridge as bb
    run = bb.BioSNICARRunner(bb.locate_biosnicar(None))
    flx = np.interp(WL, run.wvl_um * 1000, run.illumination(SZA).flx_slr)
    for k, cal in cals.items():
        s = cal.summary()
        C, D, _ = draw_macs(cal, N_DRAWS)
        macs[k] = dict(C=C, D=D, C_mean=cal.mac_C(WL), D_mean=cal.mac_D(WL))
        rows.append(dict(treatment=k, **{f"{n}_mean": s[n]["mean"] for n in ("dE", "w", "phi")},
                         **{f"{n}_sd": s[n]["sd"] for n in ("dE", "w", "phi")}, log_f_mean=s["log_f"]["mean"],
                         log_f_sd=s["log_f"]["sd"], corr=json.dumps(s["corr_dE_w_f_phi"]),
                         D_rel_spread_vis=float(np.nanmean(D[:, VIS].std(0) / D[:, VIS].mean(0))),
                         D_vs_extract_logrmse=float(np.sqrt(np.nanmean((np.log(macs[k]["D_mean"]) - np.log(meas)) ** 2))),
                         C_vs_extract_logrmse=float(np.sqrt(np.nanmean((np.log(np.maximum(macs[k]["C_mean"], 1e-9)) - np.log(meas)) ** 2))),
                         fe_term_share_vis=float(1 - np.trapezoid(macs[k]["C_mean"][VIS], WL[VIS]) / np.trapezoid(macs[k]["D_mean"][VIS], WL[VIS]))))
    a, b = macs["T-B3"], macs["T-CAM"]
    diff = {}
    for t in ("C", "D"):
        d = (a[f"{t}_mean"] - b[f"{t}_mean"]) / np.maximum(0.5 * (a[f"{t}_mean"] + b[f"{t}_mean"]), 1e-12)
        sw = flx[WL >= 300]
        diff[t] = dict(max_abs_rel=float(np.nanmax(np.abs(d))), vis_mean_abs_rel=float(np.nanmean(np.abs(d[VIS]))),
                       solar_weighted_rel=float(np.sum((a[f"{t}_mean"] - b[f"{t}_mean"])[WL >= 300] * sw) /
                                                np.sum(0.5 * (a[f"{t}_mean"] + b[f"{t}_mean"])[WL >= 300] * sw)),
                       posterior_spread_vis_B3=float(np.mean(a[t][:, VIS].std(0) / a[t][:, VIS].mean(0))))
    return pd.DataFrame(rows), diff, macs, meas


def cell_impurities(mac_wl_fn, runner_root):
    """Identical cell model for every treatment (as emulator._Builder 'ours'), only the phenolic MAC differs."""
    import cell_optics as co
    import biosnicar_bridge as bb
    from pigment_packaging import CellGeometry
    from emulator import empirical_species
    kw = co.water_k_480(runner_root)
    mac480 = co.to_480(mac_wl_fn)
    imps, pg = {}, {}
    for name, sp in empirical_species().items():
        cell = co.CellModel(CellGeometry("cylinder", sp.diameter_um / 2.0, sp.length_um),
                            co.empirical_phenolic_concentration(name), vd_diagnostic=False,
                            extra_pigments=co.empirical_pigments(species=name))
        o = cell.optics(mac480, kw, packaged=True)
        imps[name] = bb.CustomImpurity(name, o["ext_xsc"], o["ss_alb"], o["asm_prm"])
    return imps


def stage_expA(cals, macs):
    import biosnicar_bridge as bb
    import cell_optics as co
    root = bb.locate_biosnicar(None)
    run = bb.BioSNICARRunner(root)
    dust = copy.deepcopy(run.default_impurity("glacier_algae"))
    d = np.load(os.path.join(root, "data", "OP_data", "480band", "lap.npz"))
    st = "dust_greenland_Cook_CENTRAL_20190911"
    dust.name, dust.unit = "dust", 0
    dust.mac, dust.ssa, dust.g = d[st + "__ext_cff_mss"], d[st + "__ss_alb"], d[st + "__asm_prm"]
    srf = ED.s2_srf_480("S2A", ("B2", "B3", "B4", "B8"))
    sw = bb.sw_down_clear_sky(SZA)
    fn = 0.5
    meas_vivo = ED.pigment_macs_480()["phenolics_williamson2020"]          # what measured_mac_C uses
    treatments = {"T-MEAS": [lambda wl: np.interp(wl, co.WVL_480_NM, meas_vivo)]}
    for k in cals:
        for t in ("C", "D"):
            treatments[f"{k}-{t}"] = [(lambda wl, m=m: np.interp(wl, WL, m)) for m in macs[k][t]]
            treatments[f"{k}-{t}-mean"] = [lambda wl, m=macs[k][f"{t}_mean"]: np.interp(wl, WL, m)]
    rows = []
    clean = {}
    for s in STATES:
        spec = bb.IceSpec(s["r_um"], 450.0, 690.0, mode="bubbly")
        key = (s["r_um"], s["dust_ppb"])
        if key not in clean:
            a0, flx, _ = run.run_multi(spec, SZA, [(copy.copy(dust), s["dust_ppb"])])
            clean[key] = (a0, flx)
    for name, fns in treatments.items():
        for di, mfn in enumerate(fns):
            imps = cell_impurities(mfn, root)
            for s in STATES:
                spec = bb.IceSpec(s["r_um"], 450.0, 690.0, mode="bubbly")
                a0, flx = clean[(s["r_um"], s["dust_ppb"])]
                B = 10 ** s["log_b"]
                alb, flx, _ = run.run_multi(spec, SZA, [(copy.copy(imps["nordenskioeldii"]), B * fn),
                                                        (copy.copy(imps["alaskanum"]), B * (1 - fn)),
                                                        (copy.copy(dust), s["dust_ppb"])])
                bands = {b: float(np.sum(alb * srf[b] * flx) / np.sum(srf[b] * flx)) for b in srf}
                rows.append(dict(treatment=name.replace("-mean", ""), draw=("mean" if name.endswith("-mean") or name == "T-MEAS" else di),
                                 **s, bba=run.broadband(alb, flx), bba_clean=run.broadband(a0, flx),
                                 absorbed_algal_W=run.forcing(a0, alb, flx, sw), **bands))
        print(f"[expA] {name}: {len(fns)} draw(s)", flush=True)
    return pd.DataFrame(rows)


def summarise_expA(R):
    num = ["bba", "bba_clean", "absorbed_algal_W", "B2", "B3", "B4", "B8"]
    m = R[R.draw.astype(str) == "mean"].set_index(["treatment", "log_b", "r_um", "dust_ppb"])[num]
    pairs = [("T-B3-D", "T-CAM-D"), ("T-B3-C", "T-CAM-C"), ("T-B3-D", "T-MEAS"), ("T-CAM-D", "T-MEAS"), ("T-B3-D", "T-B3-C")]
    out = []
    for a, b in pairs:
        d = (m.loc[a] - m.loc[b])
        out.append(dict(pair=f"{a} - {b}", max_abs_dbba=float(d.bba.abs().max()), mean_abs_dbba=float(d.bba.abs().mean()),
                        max_abs_dB2=float(d.B2.abs().max()), max_abs_dabsorbed_W=float(d.absorbed_algal_W.abs().max()),
                        state_of_max=str(d.bba.abs().idxmax())))
    dr = R[R.draw.astype(str) != "mean"]
    spread = (dr.groupby(["treatment", "log_b", "r_um", "dust_ppb"]).agg(bba_sd=("bba", "std"), W_sd=("absorbed_algal_W", "std"))
              .reset_index().groupby("treatment").agg(max_bba_sd=("bba_sd", "max"), max_W_sd=("W_sd", "max")).reset_index())
    return pd.DataFrame(out), spread


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--stages", nargs="+", default=["raw", "calibrate", "compare", "expA"])
    a = p.parse_args(argv)
    os.makedirs(OUT, exist_ok=True)
    specs = load_specs()
    import provenance as PV
    man = dict(spec_files={k: {f: PV.file_sha256(os.path.join(d, f"level2_{fn}_{f}.csv")) for f in ("states", "spectrum")}
                           for k, (d, fn) in SPECS.items()},
               data=PV.empirical_data_fingerprint(), code=PV.file_sha256(os.path.abspath(__file__)),
               calibration_code=PV.file_sha256(os.path.join(ROOT, "phase2", "tddft_calibration.py")), environment=PV.environment())
    json.dump(man, open(os.path.join(OUT, "manifest.json"), "w"), indent=1)
    if "raw" in a.stages:
        R, meas = stage_raw(specs)
        R.to_csv(os.path.join(OUT, "raw_vs_measured.csv"), index=False, float_format="%.5g")
        json.dump(meas, open(os.path.join(OUT, "measured_features.json"), "w"), indent=1)
        print(R.to_string(index=False)); print(json.dumps(meas, indent=1))
    if a.stages == ["summarise"]:
        R = pd.read_csv(os.path.join(OUT, "expA_states.csv"))
        S, spread = summarise_expA(R)
        S.to_csv(os.path.join(OUT, "expA_pair_differences.csv"), index=False, float_format="%.5g")
        spread.to_csv(os.path.join(OUT, "expA_posterior_spread.csv"), index=False, float_format="%.5g")
        print(S.to_string(index=False)); print(spread.to_string(index=False))
        return
    if not ({"calibrate", "compare", "expA"} & set(a.stages)):
        return
    # CAM first: the B3LYP calibration is being written to the shared disk cache by the forward diagnostics
    cals = {k: calibration(specs[k]) for k in ("T-CAM", "T-B3")}
    cals = {k: cals[k] for k in specs}
    print("[calibrate] done", flush=True)
    if {"compare", "expA"} & set(a.stages):
        C, diff, macs, meas = stage_compare(cals)
        C.to_csv(os.path.join(OUT, "calibration_comparison.csv"), index=False, float_format="%.5g")
        json.dump(diff, open(os.path.join(OUT, "calibrated_mac_differences.json"), "w"), indent=1)
        pd.DataFrame({"wl_nm": WL, "extract_mac": meas, **{f"{k}_{t}_mean": macs[k][f"{t}_mean"] for k in macs for t in "CD"}}) \
            .to_csv(os.path.join(OUT, "calibrated_macs.csv"), index=False, float_format="%.5g")
        print(C.to_string(index=False)); print(json.dumps(diff, indent=1))
    if "expA" in a.stages or "summarise" in a.stages:
        if "expA" in a.stages:
            R = stage_expA(cals, macs)
            R.to_csv(os.path.join(OUT, "expA_states.csv"), index=False, float_format="%.6g")
        else:
            R = pd.read_csv(os.path.join(OUT, "expA_states.csv"))
        S, spread = summarise_expA(R)
        S.to_csv(os.path.join(OUT, "expA_pair_differences.csv"), index=False, float_format="%.5g")
        spread.to_csv(os.path.join(OUT, "expA_posterior_spread.csv"), index=False, float_format="%.5g")
        print(S.to_string(index=False)); print(spread.to_string(index=False))


if __name__ == "__main__":
    main()
