#!/usr/bin/env python
"""
Uncertainty-aware contrasts for the controlled molecular-substitution experiment (Experiment A).

The original summary (`molecular_contribution.summarise_expA`, files expA_pair_differences.csv /
expA_posterior_spread.csv) is RETAINED as plug-in point-estimate contrasts (posterior-mean spectra). This
module adds posterior contrasts:

  * Within one functional (tier C vs tier D): PAIRED by posterior sample (both tiers are evaluated from the
    same calibration draw: same dE, w, f; D adds phi x the measured Fe increment).
  * Between separately calibrated functionals (B3LYP vs CAM-B3LYP): the two calibrations are independent
    posteriors with no shared parameter, so equal draw indices carry no meaning. Contrasts use INDEPENDENT
    draws (all n_a x n_b combinations = product of the two marginal posteriors). A common-cause coupling
    (both calibrations fit the same measured spectra) is not represented; this is stated as a limitation.
  * TD-DFT vs measured MAC: TD draws against measured-MAC draws. Measured-MAC draws (stage `extra`) perturb
    the Williamson (2020) extract MAC by its regression SE with ONE standard normal per draw applied at all
    wavelengths (fully correlated amplitude error). Per-wavelength correlation of the SE is unknown;
    full correlation maximises the broadband effect. Without the `extra` stage the measured MAC is a point
    value and its uncertainty is reported as omitted.

Outputs (records/molecular_contribution/uncertainty/): per-state contrast quantiles, P(|delta| >=
practical threshold), P(delta > 0), and a Monte Carlo stability table (24 vs all draws).

Practical thresholds (from the frozen protocol, not tuned): broadband albedo 0.01, absorbed SW 10 W m-2.
These are practical-importance thresholds; they are NOT measurement uncertainties or detection limits.

    python3 phase4/molecular_uncertainty.py --stages contrasts            # from saved expA_states.csv only
    python common/run_budgeted.py --background --mem-mb 900 --name mol_extra -- \
        python3 phase4/molecular_uncertainty.py --stages extra contrasts --n-extra 72
"""
from __future__ import annotations

import argparse
import itertools
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import molecular_contribution as MC  # noqa: E402

OUT = os.path.join(MC.OUT, "uncertainty")
THR = {"bba": 0.01, "absorbed_algal_W": 10.0}
KEYS = ["log_b", "r_um", "dust_ppb"]
QTY = ["bba", "absorbed_algal_W", "B2", "B4", "B8"]


# --------------------------------------------------------------------------- extra draws (resumable)
def extra_draws(n_extra, seed=1):
    import copy
    import pickle
    import biosnicar_bridge as bb
    import cell_optics as co
    import empirical_data as ED
    import tddft_calibration as TC
    os.makedirs(OUT, exist_ok=True)
    specs = MC.load_specs()
    cals = {k: MC.calibration(specs[k]) for k in ("T-B3", "T-CAM")}
    root = bb.locate_biosnicar(None)
    run = bb.BioSNICARRunner(root)
    dust = copy.deepcopy(run.default_impurity("glacier_algae"))
    d = np.load(os.path.join(root, "data", "OP_data", "480band", "lap.npz"))
    st = "dust_greenland_Cook_CENTRAL_20190911"
    dust.name, dust.unit = "dust", 0
    dust.mac, dust.ssa, dust.g = d[st + "__ext_cff_mss"], d[st + "__ss_alb"], d[st + "__asm_prm"]
    srf = ED.s2_srf_480("S2A", ("B2", "B3", "B4", "B8"))
    sw = bb.sw_down_clear_sky(MC.SZA)
    clean = {}

    def evaluate(mfn, treatment, draw):
        imps = MC.cell_impurities(mfn, root)
        rows = []
        for s in MC.STATES:
            spec = bb.IceSpec(s["r_um"], 450.0, 690.0, mode="bubbly")
            key = (s["r_um"], s["dust_ppb"])
            if key not in clean:
                clean[key] = run.run_multi(spec, MC.SZA, [(copy.copy(dust), s["dust_ppb"])])[:2]
            a0, _ = clean[key]
            B = 10 ** s["log_b"]
            alb, flx, _ = run.run_multi(spec, MC.SZA, [(copy.copy(imps["nordenskioeldii"]), B * 0.5),
                                                       (copy.copy(imps["alaskanum"]), B * 0.5),
                                                       (copy.copy(dust), s["dust_ppb"])])
            bands = {b: float(np.sum(alb * srf[b] * flx) / np.sum(srf[b] * flx)) for b in srf}
            rows.append(dict(treatment=treatment, draw=draw, **s, bba=run.broadband(alb, flx),
                             bba_clean=run.broadband(a0, flx), absorbed_algal_W=run.forcing(a0, alb, flx, sw), **bands))
        return rows

    def done(path):
        return set(pd.read_csv(path).draw.astype(int)) if os.path.isfile(path) else set()

    # TD treatments: draws disjoint from the original 24 (same rng/seed as molecular_contribution.draw_macs)
    for k, cal in cals.items():
        n = cal.samples.shape[0]
        orig = np.random.default_rng(0).choice(n, MC.N_DRAWS, replace=False)
        pool = np.setdiff1d(np.arange(n), orig)
        extra = np.random.default_rng(seed).choice(pool, n_extra, replace=False)
        path = os.path.join(OUT, f"extra_{k}.csv")
        have = done(path)
        for j, i in enumerate(extra):
            draw = MC.N_DRAWS + j
            if draw in have:
                continue
            dE, w, f, phi = cal.samples[i, :4]
            mC = TC.perturbed_mac(cal.spec, dE, f, w)(MC.WL)
            mD = TC.complexed_mac(cal.spec, dE, f, w, phi)(MC.WL)
            rows = []
            for t, m in (("C", mC), ("D", mD)):
                rows += evaluate(lambda wl, m=m: np.interp(wl, MC.WL, m), f"{k}-{t}", draw)
            pd.DataFrame(rows).assign(posterior_index=int(i)).to_csv(path, mode="a", header=not os.path.isfile(path),
                                                                      index=False, float_format="%.6g")
        print(f"[extra] {k}: {n_extra} additional draws", flush=True)
    # measured MAC draws: amplitude perturbed by the regression SE (one z per draw, all wavelengths)
    wl_m, E, sE = ED.phenolic_extract_mac()
    meas = ED.pigment_macs_480()["phenolics_williamson2020"]
    rel_se = np.interp(co.WVL_480_NM, wl_m, sE / E, left=0.0, right=0.0)
    path = os.path.join(OUT, "extra_T-MEAS.csv")
    have = done(path)
    z = np.random.default_rng(seed + 1).standard_normal(MC.N_DRAWS + n_extra)
    for draw, zz in enumerate(z):
        if draw in have:
            continue
        m = np.clip(meas * (1.0 + zz * rel_se), 0.0, None)
        pd.DataFrame(evaluate(lambda wl, m=m: np.interp(wl, co.WVL_480_NM, m), "T-MEAS", draw)).assign(z=zz).to_csv(
            path, mode="a", header=not os.path.isfile(path), index=False, float_format="%.6g")
    print("[extra] T-MEAS draws done", flush=True)


# --------------------------------------------------------------------------- contrasts
def load_states():
    R = pd.read_csv(os.path.join(MC.OUT, "expA_states.csv"))
    R["draw"] = R.draw.astype(str)
    parts = [R]
    for f in ("extra_T-B3.csv", "extra_T-CAM.csv", "extra_T-MEAS.csv"):
        p = os.path.join(OUT, f)
        if os.path.isfile(p):
            x = pd.read_csv(p)
            x["draw"] = x.draw.astype(int).astype(str)
            parts.append(x[R.columns])
    A = pd.concat(parts, ignore_index=True)
    return A.drop_duplicates(["treatment", "draw"] + KEYS, keep="first")


def draws_of(A, treatment, max_draws=None):
    x = A[(A.treatment == treatment) & (A.draw != "mean")].copy()
    x["d"] = x.draw.astype(int)
    if max_draws is not None:
        x = x[x.d < max_draws]
    return x


def summarise(delta, thr):
    delta = np.asarray(delta, float)
    return dict(n=int(delta.size), median=float(np.median(delta)), q025=float(np.quantile(delta, 0.025)),
                q975=float(np.quantile(delta, 0.975)), p_pos=float(np.mean(delta > 0)),
                p_exceeds_threshold=float(np.mean(np.abs(delta) >= thr)) if thr is not None else np.nan)


def contrasts(A, max_draws=None):
    plug = A[A.draw == "mean"].set_index(["treatment"] + KEYS)
    pairs = [("T-B3-D", "T-B3-C", "paired"), ("T-CAM-D", "T-CAM-C", "paired"),
             ("T-B3-D", "T-CAM-D", "independent"), ("T-B3-C", "T-CAM-C", "independent"),
             ("T-B3-D", "T-MEAS", "vs_measured"), ("T-CAM-D", "T-MEAS", "vs_measured")]
    meas_draws = draws_of(A, "T-MEAS", max_draws)
    rows = []
    for a, b, kind in pairs:
        da, db = draws_of(A, a, max_draws), draws_of(A, b, max_draws)
        for s in MC.STATES:
            sel = lambda x: x[(x.log_b == s["log_b"]) & (x.r_um == s["r_um"]) & (x.dust_ppb == s["dust_ppb"])]
            xa, xb = sel(da).set_index("d"), sel(db).set_index("d")
            for q in QTY:
                plug_d = float(plug.loc[(a, *s.values())][q] - plug.loc[(b, *s.values())][q])
                if kind == "paired":
                    common = xa.index.intersection(xb.index)
                    delta = (xa.loc[common, q] - xb.loc[common, q]).to_numpy()
                    meas_unc = "n/a"
                elif kind == "independent":
                    delta = (xa[q].to_numpy()[:, None] - xb[q].to_numpy()[None, :]).ravel()
                    meas_unc = "n/a"
                else:
                    xm = sel(meas_draws)
                    if len(xm):
                        delta = (xa[q].to_numpy()[:, None] - xm[q].to_numpy()[None, :]).ravel()
                        meas_unc = "included (regression SE, fully correlated amplitude)"
                    else:
                        delta = xa[q].to_numpy() - float(plug.loc[(b, *s.values())][q])
                        meas_unc = "OMITTED (point value)"
                rows.append(dict(pair=f"{a} - {b}", coupling=kind, **s, quantity=q, plug_in=plug_d,
                                 measured_mac_uncertainty=meas_unc, **summarise(delta, THR.get(q))))
    return pd.DataFrame(rows)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--stages", nargs="+", default=["contrasts"])
    p.add_argument("--n-extra", type=int, default=72)
    a = p.parse_args(argv)
    os.makedirs(OUT, exist_ok=True)
    if "extra" in a.stages:
        extra_draws(a.n_extra)
    if "contrasts" in a.stages:
        A = load_states()
        C = contrasts(A)
        C.to_csv(os.path.join(OUT, "posterior_contrasts.csv"), index=False, float_format="%.5g")
        n_td = draws_of(A, "T-B3-D").d.nunique()
        stab = None
        if n_td > MC.N_DRAWS:
            C24 = contrasts(A, max_draws=MC.N_DRAWS)
            m = C.merge(C24, on=["pair"] + KEYS + ["quantity"], suffixes=("_all", "_24"))
            stab = (m.assign(d_q025=(m.q025_all - m.q025_24).abs(), d_q975=(m.q975_all - m.q975_24).abs(),
                             d_p=(m.p_exceeds_threshold_all - m.p_exceeds_threshold_24).abs())
                    .groupby(["pair", "quantity"]).agg(max_d_q025=("d_q025", "max"), max_d_q975=("d_q975", "max"),
                                                       max_d_p_exceed=("d_p", "max")).reset_index())
            stab.to_csv(os.path.join(OUT, "mc_stability_24_vs_all.csv"), index=False, float_format="%.5g")
        S = (C[C.quantity.isin(THR)].groupby(["pair", "coupling", "quantity"])
             .agg(max_abs_plug_in=("plug_in", lambda x: float(np.max(np.abs(x)))),
                  min_q025=("q025", "min"), max_q975=("q975", "max"),
                  max_p_exceeds_threshold=("p_exceeds_threshold", "max"), n_draw_pairs=("n", "max"),
                  measured_mac_uncertainty=("measured_mac_uncertainty", "first")).reset_index())
        S.to_csv(os.path.join(OUT, "posterior_contrast_summary.csv"), index=False, float_format="%.5g")
        json.dump(dict(n_draws_per_TD_treatment=int(n_td), thresholds=THR,
                       n_measured_draws=int(draws_of(A, "T-MEAS").d.nunique()),
                       coupling_note=__doc__.split("Outputs")[0]), open(os.path.join(OUT, "settings.json"), "w"), indent=1)
        pd.set_option("display.width", 250)
        print(S.to_string(index=False))
        if stab is not None:
            print(stab.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
