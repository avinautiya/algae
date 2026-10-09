#!/usr/bin/env python
"""
H1 (docs/glacier_model_validation_protocol.md): does a pigment/cell representation predict independent
plot-scale broadband ALBEDO better than established algae optics, with abundance GIVEN (measured counts)?

Target: S6 2017 hemispherical spectral albedo (Cook et al. 2020 archive, biosnicar-py Albedo_master.csv,
not used in any ALBEDO fit; NB it entered the HCRF/albedo k prior via ARF_master, which H1 does not use - deviation D3) integrated over 300-2500 nm with the BioSNICAR clear-sky
irradiance at the plot's solar zenith (300-350 nm held at the 350 nm value) - the same weighting as the
model's broadband albedo (BBA).

Models: the frozen held-out emulators (records/heldout_v2/FROZEN_PROVENANCE.json: posterior-mean optics):
  M0 no algae (log B at the lowest node is NOT zero algae -> the no-algae albedo is BBA - rf_algae/SW),
  M1 tierA_empirical, M1b measured_mac_C, M3 tddft_D, M3-alt tddft_C / tddft_D_iid.
Nuisance state (ice radius population, f_n, dust) is marginalised under priors; the radius population
prior (log-normal mu, sd) and an albedo discrepancy SD are chosen by empirical Bayes on the TRAINING
days' albedos only (leave-one-day-out). Nothing is fitted to a scored plot.

Outputs (outdir): per-plot predictions, per-fold settings, summary with day-block bootstrap contrasts
(n_blocks reported; < 5 blocks -> no interval), coverage AND sharpness, CRPS.

    python phase4/h1_albedo.py --outdir phase4/results/h1_albedo
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.special import logsumexp
from scipy.stats import norm

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "phase2"))

import empirical_data as ED  # noqa: E402

OPTICS = ("no_algae", "tierA_empirical", "measured_mac_C", "tddft_D", "tddft_C", "tddft_D_iid")
PRIMARY, CONTRASTS = "tddft_D", ("tierA_empirical", "measured_mac_C")
EMU_DIR = os.path.join(HERE, "results", "heldout_v2", "cache")
MU_LNR = np.log([400.0, 600.0, 900.0, 1350.0, 2000.0, 3000.0, 4500.0, 7000.0])
SD_LNR = np.array([0.35, 0.6, 0.9, 1.2, 1.6])
DISC_SD = np.array([0.005, 0.01, 0.015, 0.02, 0.03, 0.045, 0.07, 0.1])     # albedo discrepancy + obs SD
MIN_IMPROVEMENT = 0.01


def albedo_data(biosnicar_root):
    """Plots with a finite albedo spectrum and a count; exclusion reasons recorded."""
    a = pd.read_csv(os.path.join(biosnicar_root, "data", "additional_data", "Albedo_master.csv"))
    c = pd.read_csv(ED.path("cook2020_archive_cell_counts.csv"))
    wl = a.Wavelength.to_numpy(float)
    rows, excluded = [], []
    for r in c.itertuples():
        if r.sample not in a:
            excluded.append(dict(sample=r.sample, reason="no albedo spectrum"))
            continue
        sp = a[r.sample].to_numpy(float)
        vis = (wl >= 400) & (wl <= 1300)
        if np.mean(~np.isfinite(sp)) > 0.05:
            excluded.append(dict(sample=r.sample, reason=">5% missing values"))
            continue
        if np.nanmax(sp[vis]) > 1.05 or np.nanmin(sp[vis]) < 0:
            excluded.append(dict(sample=r.sample, reason="albedo outside [0, 1.05] in 400-1300 nm"))
            continue
        day, month = (int(x) for x in r.sample.split("_")[:2])
        doy = pd.Timestamp(year=2017, month=month, day=day).dayofyear
        rows.append(dict(sample=r.sample, day=f"{day}_{month}", cells=float(r.cells_per_ml),
                         cells_counted=float(r.cells_counted), sza=ED.solar_zenith_noon(doy)))
    return pd.DataFrame(rows), pd.DataFrame(excluded), a, wl


def broadband(spec, wl, sza, runner):
    grid = np.arange(300.0, 2500.0 + 0.5, 1.0)
    ok = np.isfinite(spec)
    h = np.interp(grid, wl[ok], spec[ok])
    flx = np.interp(grid, runner.wvl_um * 1000.0, runner.illumination(sza).flx_slr)
    return float(np.sum(h * flx) / np.sum(flx))


def load_emulators(optics, szas):
    import emulator as E
    src = "tierA_empirical" if optics == "no_algae" else optics
    out = {}
    for z in szas:
        e = E.Emulator.load(os.path.join(EMU_DIR, f"heldout_{src}_sza{int(z)}.npz"))
        out[int(z)] = e.refine_log_b(0.05)
    return out


def node_bba(em, log_b, no_algae=False):
    """BBA on the (f_n, r, dust) grid at abundance log_b (interpolated between log B nodes); for the
    no-algae model the algae-free albedo BBA + rf_algae / SW (rf_algae is relative to an algae-free run)."""
    lb = em.axes["log_b"]
    j = int(np.clip(np.searchsorted(lb, log_b) - 1, 0, len(lb) - 2))
    t = float(np.clip((log_b - lb[j]) / (lb[j + 1] - lb[j]), 0, 1))
    b = (1 - t) * em.data["bba"][j] + t * em.data["bba"][j + 1]
    if no_algae:
        rf = (1 - t) * em.data["rf_algae"][j] + t * em.data["rf_algae"][j + 1]
        b = b + rf / em.meta["sw_down"]
    return b                                               # (f_n, r, dust)


def prior_weights(em, mu, sd):
    from priors import PriorConfig, prior_logpdfs
    import dataclasses
    pc = dataclasses.replace(PriorConfig.for_density(690.0), mu_lnr=mu, sd_lnr=sd)
    lp, _, _, _ = prior_logpdfs(em.axes, 1, pc)
    w = sum(lp[a][0].reshape([-1 if b == a else 1 for b in em.names[1:]]) for a in em.names[1:])
    w = np.exp(w - w.max())
    return w / w.sum()


def plot_logp(bba_grid, w, y, s):
    """log p(y | mixture over nuisance nodes, Gaussian discrepancy SD s)."""
    return float(logsumexp(np.log(np.maximum(w, 1e-300)) + norm.logpdf(y, bba_grid, s)))


def predictive(bba_grid, w, s, q=(0.05, 0.5, 0.95)):
    """Mixture mean, quantiles and CRPS-ready samples of the predictive broadband albedo."""
    x = np.linspace(0.0, 1.0, 2001)
    pdf = (w.reshape(-1)[:, None] * norm.pdf(x[None, :], bba_grid.reshape(-1)[:, None], s)).sum(axis=0)
    cdf = np.cumsum(pdf)
    cdf /= cdf[-1]
    return x, cdf, float(np.sum(w * bba_grid)), [float(x[np.argmax(cdf >= qq)]) for qq in q]


def crps(x, cdf, y):
    return float(np.sum((cdf - (x >= y)) ** 2) * (x[1] - x[0]))


def run(biosnicar=None, outdir=None):
    import biosnicar_bridge as bb
    root = bb.locate_biosnicar(biosnicar)
    runner = bb.BioSNICARRunner(root)
    df, excl, a, wl = albedo_data(root)
    df["bba_meas"] = [broadband(a[s].to_numpy(float), wl, z, runner) for s, z in zip(df["sample"], df.sza)]
    df["log_b"] = np.log10(np.clip(df.cells, 10.0, None))      # zero counts -> lowest node (10 cells/mL)
    szas = sorted(set(np.round(df.sza).astype(int)))
    days = sorted(df.day.unique(), key=lambda d: (int(d.split("_")[1]), int(d.split("_")[0])))
    preds, settings = [], {}
    for opt in OPTICS:
        ems = load_emulators(opt, szas)
        grids = [node_bba(ems[int(round(z))], lb, opt == "no_algae") for z, lb in zip(df.sza, df.log_b)]
        wcache = {}

        def W(z, mu, sd):
            k = (z, mu, sd)
            if k not in wcache:
                wcache[k] = prior_weights(ems[z], mu, sd)
            return wcache[k]
        # log likelihood of every plot under every (mu, sd, s)
        hyp = [(m, d, s) for m in MU_LNR for d in SD_LNR for s in DISC_SD]
        ll = np.array([[plot_logp(g, W(int(round(z)), m, d), y, s) for (m, d, s) in hyp]
                       for g, z, y in zip(grids, df.sza, df.bba_meas)])
        for day in days:
            tr = (df.day != day).to_numpy()
            h = int(np.argmax(ll[tr].sum(axis=0)))
            m, d, s = hyp[h]
            settings[f"{opt}/{day}"] = dict(radius_median_um=float(np.exp(m)), radius_ln_sd=float(d), disc_sd=float(s))
            for i in np.flatnonzero(~tr):
                z = int(round(df.sza.iloc[i]))
                x, cdf, mean, (q05, q50, q95) = predictive(grids[i], W(z, m, d), s)
                y = df.bba_meas.iloc[i]
                preds.append(dict(model=opt, sample=df["sample"].iloc[i], day=day, cells=df.cells.iloc[i],
                                  bba_meas=y, pred_mean=mean, pred_median=q50, q05=q05, q95=q95,
                                  abs_err=abs(mean - y), err=mean - y, crps=crps(x, cdf, y),
                                  covered90=float(q05 <= y <= q95), width90=q95 - q05,
                                  logp=plot_logp(grids[i], W(z, m, d), y, s)))
        print(f"[H1] {opt} done", flush=True)
    P = pd.DataFrame(preds)
    summ = (P.groupby("model").agg(n=("sample", "size"), mae=("abs_err", "mean"), bias=("err", "mean"),
                                   crps=("crps", "mean"), coverage90=("covered90", "mean"),
                                   width90=("width90", "mean"), mean_logp=("logp", "mean")).reset_index())
    contr = []
    for c in CONTRASTS + ("no_algae",):
        a_ = P[P.model == c].set_index("sample")
        b_ = P[P.model == PRIMARY].set_index("sample").reindex(a_.index)
        dmae = (a_.abs_err - b_.abs_err).to_numpy()          # positive = primary better
        blocks = a_.day.to_numpy()
        ub = np.unique(blocks)
        rng = np.random.default_rng(0)
        bs = [np.mean(np.concatenate([dmae[blocks == u] for u in rng.choice(ub, ub.size)])) for _ in range(4000)]
        contr.append(dict(contrast=f"MAE({c}) - MAE({PRIMARY})", mean=float(dmae.mean()),
                          lo=float(np.quantile(bs, 0.025)) if ub.size >= 5 else np.nan,
                          hi=float(np.quantile(bs, 0.975)) if ub.size >= 5 else np.nan,
                          n=int(dmae.size), n_blocks=int(ub.size),
                          per_block=json.dumps({u: round(float(dmae[blocks == u].mean()), 4) for u in ub})))
    C = pd.DataFrame(contr)
    if outdir:
        os.makedirs(outdir, exist_ok=True)
        import provenance as PV
        PV.atomic_to_csv(P, os.path.join(outdir, "h1_predictions.csv"), index=False, float_format="%.5g")
        PV.atomic_to_csv(summ, os.path.join(outdir, "h1_summary.csv"), index=False, float_format="%.5g")
        PV.atomic_to_csv(C, os.path.join(outdir, "h1_contrasts.csv"), index=False, float_format="%.5g")
        PV.atomic_to_csv(excl, os.path.join(outdir, "h1_excluded.csv"), index=False)
        PV.atomic_write_text(os.path.join(outdir, "h1_settings.json"), json.dumps(dict(
            settings=settings, n_plots=int(len(df)), n_days=len(days), min_improvement=MIN_IMPROVEMENT,
            emulators=EMU_DIR, environment=PV.environment()), indent=1, default=float))
    return df, P, summ, C


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--biosnicar", default=None)
    p.add_argument("--outdir", default=os.path.join(HERE, "results", "h1_albedo"))
    a = p.parse_args(argv)
    df, P, summ, C = run(a.biosnicar, a.outdir)
    pd.set_option("display.width", 200)
    print(summ.round(4).to_string(index=False))
    print(C.round(4).to_string(index=False))


if __name__ == "__main__":
    main()
