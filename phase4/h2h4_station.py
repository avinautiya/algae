#!/usr/bin/env python
"""
H2 (absorbed SW) and H4 (albedo beyond empirical satellite albedo) at PROMICE KAN_L / KAN_M, per the
frozen protocol (docs/glacier_model_validation_protocol.md §3, §5). Input: phase4/station_pixels.py.

Predictions per valid station-day (identical pixel, mask, hours and SW-down for every model):

  * M3, M1, M1b, M3-alt: posterior-mean broadband albedo of the pixel from the exact grid posterior.
    Frozen held-out emulators and the primary-fold field calibration (sigma, radius prior); nothing
    is fitted to station data.
  * M0: the same retrieval with algae fixed at the lowest node (10 cells/mL), i.e. the albedo the
    pixel's ice + dust state implies without algae.
  * M2: "simple empirical" (the published Naegeli et al. 2017 coefficients could not be obtained
    verbatim: MDPI 403, repository mirror behind bot protection; protocol fallback):
        alpha_M2 = sum_b w_b R_b,  b in {B2, B3, B4, B8}
    with w_b the share of the BioSNICAR clear-sky irradiance (300-2500 nm, at the scene SZA) closest
    in wavelength to band b (centres 492, 560, 665, 833 nm). Fixed, no fitting.
  * Absorbed SW prediction = mean SW-down(3 h) x (1 - alpha).

Information used (reported, protocol §6):
  * M3/M1/M1b/M0 use a k prior on HCRF/albedo (field ARF, deviation D3), field-calibrated sigma and the
    radius prior;
  * M2 uses none of these.

SZA: the scene SZA is rounded to the nearest frozen emulator node; nodes exist at 44-50, 52, 53 and 56.
The nearest node within 1 degree is used and its offset recorded; rows further away are excluded. A
sensitivity restricted to exact nodes is reported.

Blocks: station-year. Primary: 2016-2023 except 2019 (amendment A1); 2019 is reported separately.
With fewer than 5 blocks no interval is given (protocol §6).

    python phase4/h2h4_station.py --pixels phase4/results/station_pixels/station_pixels.csv \
        --outdir phase4/results/h2h4
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "phase2"))

import scene_product as SP  # noqa: E402

PHYS_MODELS = {"M3": "tddft_D", "M1": "tierA_empirical", "M1b": "measured_mac_C", "M3alt_C": "tddft_C",
               "M3alt_iid": "tddft_D_iid"}
BAND_CENTRES = {"B2": 492.0, "B3": 560.0, "B4": 665.0, "B8": 833.0}
PRIMARY_YEARS = (2016, 2017, 2018, 2020, 2021, 2022, 2023)   # protocol amendment A1
MIN_H2, MIN_H4 = 10.0, 0.01
MAX_NODE_OFFSET = 1.0


def available_nodes(opt):
    return sorted(int(f.split("sza")[1][:-4]) for f in os.listdir(SP.EMU_DIR)
                  if f.startswith(f"heldout_{opt}_sza"))


def m2_weights(runner, sza):
    wl = runner.wvl_um * 1000.0
    flx = runner.illumination(sza).flx_slr * runner.band
    c = np.array(list(BAND_CENTRES.values()))
    nearest = np.argmin(np.abs(wl[:, None] - c[None, :]), axis=1)
    w = np.array([flx[nearest == i].sum() for i in range(len(c))])
    return w / w.sum()


def predict(df, biosnicar=None):
    import emulator as E
    import biosnicar_bridge as bb
    runner = bb.BioSNICARRunner(bb.locate_biosnicar(biosnicar))
    s = json.load(open(SP.SETTINGS))["settings"]
    nodes = available_nodes("tddft_D")
    df = df.copy()
    df["node"] = [min(nodes, key=lambda n: abs(n - z)) for z in df.sza]
    df["node_offset"] = df.node - df.sza
    df["excluded_node"] = df.node_offset.abs() > MAX_NODE_OFFSET
    R = df[["B2", "B3", "B4", "B8"]].to_numpy(float)
    df["alpha_M2"] = [float(m2_weights(runner, z) @ r) for z, r in zip(df.sza, R)]
    for key, opt in PHYS_MODELS.items():
        df[f"alpha_{key}"] = np.nan
        if key == "M3":
            df["alpha_M0"] = np.nan
        for n in sorted(df.node[~df.excluded_node].unique()):
            idx = np.flatnonzero((df.node == n) & ~df.excluded_node)
            em = SP.add_indicator(E.Emulator.load(os.path.join(SP.EMU_DIR, f"heldout_{opt}_sza{n}.npz")).refine_log_b(0.05))
            cal = s[f"{SP.CAL_FOLD}/{opt}"]
            r, _ = SP.posterior(em, R[idx], cal)
            df.loc[df.index[idx], f"alpha_{key}"] = r["bba_mean"]
            if key == "M3":
                r0, _ = SP.posterior(em, R[idx], cal, restrict_no_algae=True)
                df.loc[df.index[idx], "alpha_M0"] = r0["bba_mean"]
                df.loc[df.index[idx], "M3_log_b_mean"] = r["log_b_mean"]
                df.loc[df.index[idx], "M3_ppp"] = r["ppp"]
        print(f"[H2/H4] {key} done", flush=True)
    for m in ["M0", "M2"] + list(PHYS_MODELS):
        df[f"absorbed_{m}"] = df.sw_down * (1 - df[f"alpha_{m}"])
        df[f"abserr_absorbed_{m}"] = (df[f"absorbed_{m}"] - df.absorbed_obs).abs()
        df[f"abserr_albedo_{m}"] = (df[f"alpha_{m}"] - df.albedo_obs).abs()
        df[f"err_albedo_{m}"] = df[f"alpha_{m}"] - df.albedo_obs
    return df


def block_contrast(d, a, b, col, min_blocks=5, n=4000, seed=0):
    """mean(err_a - err_b) (positive = b better) with a station-year cluster bootstrap."""
    diff = (d[f"{col}_{a}"] - d[f"{col}_{b}"]).to_numpy()
    blocks = (d.station + "_" + d.year.astype(str)).to_numpy()
    ub = np.unique(blocks)
    out = dict(contrast=f"{col}: {a} - {b}", mean=float(np.mean(diff)), n=int(diff.size), n_blocks=int(ub.size),
               per_block=json.dumps({u: round(float(diff[blocks == u].mean()), 4) for u in ub}))
    if ub.size >= min_blocks:
        rng = np.random.default_rng(seed)
        bs = [np.mean(np.concatenate([diff[blocks == u] for u in rng.choice(ub, ub.size)])) for _ in range(n)]
        out.update(lo=float(np.quantile(bs, 0.025)), hi=float(np.quantile(bs, 0.975)))
    else:
        out.update(lo=np.nan, hi=np.nan, note="insufficient evidence for a generalisation interval (< 5 blocks)")
    return out


def evaluate(df):
    ok = df[~df.excluded_node & df.alpha_M3.notna()]
    sets = {"primary": ok[ok.year.isin(PRIMARY_YEARS)], "2019_separate": ok[ok.year == 2019],
            "primary_exact_nodes": ok[ok.year.isin(PRIMARY_YEARS) & (ok.node_offset.abs() < 0.5)]}
    summ, contr = [], []
    for name, d in sets.items():
        if d.empty:
            continue
        for m in ["M0", "M1", "M1b", "M2", "M3", "M3alt_C", "M3alt_iid"]:
            summ.append(dict(set=name, model=m, n=len(d), n_blocks=int((d.station + d.year.astype(str)).nunique()),
                             mae_absorbed=d[f"abserr_absorbed_{m}"].mean(), mae_albedo=d[f"abserr_albedo_{m}"].mean(),
                             bias_albedo=d[f"err_albedo_{m}"].mean()))
        for a in ("M1", "M2", "M0", "M1b"):
            contr.append(dict(set=name, hypothesis="H2", **block_contrast(d, a, "M3", "abserr_absorbed")))
        contr.append(dict(set=name, hypothesis="H4", **block_contrast(d, "M2", "M3", "abserr_albedo")))
    S, C = pd.DataFrame(summ), pd.DataFrame(contr)
    verdict = {}
    p = C[C.set == "primary"]
    for hyp, a, mn in (("H2_vs_M1", "M1", MIN_H2), ("H2_vs_M2", "M2", MIN_H2), ("H4", "M2", MIN_H4)):
        col = "abserr_albedo" if hyp == "H4" else "abserr_absorbed"
        row = p[p.contrast == f"{col}: {a} - M3"]
        if row.empty:
            verdict[hyp] = "UNTESTED (no valid station-days)"
            continue
        r = row.iloc[0]
        if not np.isfinite(r.lo):
            verdict[hyp] = f"INSUFFICIENT EVIDENCE ({r.n_blocks} blocks < 5); point {r['mean']:.4g}"
        elif r.lo > 0 and r["mean"] >= mn:
            verdict[hyp] = f"SUPPORTED (point {r['mean']:.4g} >= {mn}, interval {r.lo:.4g}..{r.hi:.4g})" + \
                ("" if hyp == "H4" else "; requires H1 not contradicted")
        else:
            verdict[hyp] = f"NOT SUPPORTED (point {r['mean']:.4g}, interval {r.lo:.4g}..{r.hi:.4g}, minimum {mn})"
    return S, C, verdict


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pixels", default=os.path.join(HERE, "results", "station_pixels", "station_pixels.csv"))
    p.add_argument("--outdir", default=os.path.join(HERE, "results", "h2h4"))
    p.add_argument("--biosnicar", default=None)
    a = p.parse_args(argv)
    px = pd.read_csv(a.pixels)
    df = predict(px[px.status == "valid"].reset_index(drop=True), a.biosnicar)
    S, C, verdict = evaluate(df)
    import provenance as PV
    os.makedirs(a.outdir, exist_ok=True)
    PV.atomic_to_csv(df, os.path.join(a.outdir, "h2h4_predictions.csv"), index=False, float_format="%.6g")
    PV.atomic_to_csv(S, os.path.join(a.outdir, "h2h4_summary.csv"), index=False, float_format="%.5g")
    PV.atomic_to_csv(C, os.path.join(a.outdir, "h2h4_contrasts.csv"), index=False, float_format="%.5g")
    PV.atomic_write_text(os.path.join(a.outdir, "h2h4_verdict.json"), json.dumps(dict(
        verdict=verdict, pixels_sha256=PV.file_sha256(a.pixels), m2="simple empirical (protocol fallback)",
        environment=PV.environment()), indent=1))
    pd.set_option("display.width", 220)
    print(S.round(4).to_string(index=False))
    print(C.drop(columns=["per_block"]).round(4).to_string(index=False))
    print(json.dumps(verdict, indent=1))


if __name__ == "__main__":
    main()
