"""
Field validation against co-located measurements (Cook et al. 2020, S6, July 2017).

Data: 31 surface-ice samples, each with a nadir HCRF spectrum (ASD FieldSpec Pro 3, 10 deg lens,
clear sky, solar noon +/- 2 h) and a haemocytometer count of the algal cells in the ice that was
inside the spectrometer's field of view (cells mL^-1 of melted ice). Spectra are band-averaged with
the official ESA Sentinel-2A spectral response functions, so the test exercises exactly the
retrieval applied to Sentinel-2 pixels - minus the atmosphere and the 10-20 m footprint.

Methods compared (all evaluated out of sample):
  physics-informed Bayesian (ours)   emulator + grid posterior; reflectance noise sigma chosen by
                                     maximising the marginal likelihood of the OTHER 30 samples
  Tier A Bayesian                    same machinery with BioSNICAR's empirical algal optics
  empirical band-ratio regression    log10 B = c0 + c1 I + c2 I^2, I = (B4-B2)/(B4+B2), fitted to the
                                     other samples (leave-one-out)
  Cook et al. (2020) inversion       the published BioSNICAR-GO retrieval for the same spectra
                                     (column algae_cells_inv_model_S2 of their metadata), as reported

Zero-count samples (below the counting resolution of ~62.5 cells mL^-1) cannot be scored on a log
scale; they are reported separately.
"""

from __future__ import annotations

import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "phase2"))
sys.path.insert(0, HERE)

import empirical_data as ED  # noqa: E402

SIGMA_GRID = np.array([0.01, 0.015, 0.02, 0.03, 0.04, 0.06, 0.08])


def field_band_reflectance(spacecraft="S2A", bands=("B2", "B3", "B4", "B8")):
    m, h = ED.field_samples()
    srf = pd.read_csv(ED.path("esa_s2_srf_TN-15-0007_v4.0.csv")).set_index("wavelength_nm")
    wl = h.wavelength_nm.to_numpy()
    rows = []
    for _, r in m.iterrows():
        sp = h[r.filename].to_numpy(float)
        vals = []
        for b in bands:
            s = srf[f"{spacecraft}_{b}"].reindex(wl).fillna(0).to_numpy()
            ok = np.isfinite(sp)
            vals.append(np.sum(sp[ok] * s[ok]) / s[ok].sum())
        day, month = (int(x) for x in r.filename.split("_")[:2])
        doy = pd.Timestamp(year=ED.FIELD_SITE["year"], month=month, day=day).dayofyear
        rows.append(dict(sample=r.filename, cells=r.measured_cells, cook2020_inversion=r.algae_cells_inv_model_S2,
                         sza=ED.solar_zenith_noon(doy), **dict(zip(bands, vals))))
    return pd.DataFrame(rows)


def _metrics(t, e, lo=None, hi=None):
    from scipy.stats import spearmanr
    t, e = np.asarray(t, float), np.asarray(e, float)
    ok = np.isfinite(t) & np.isfinite(e)
    t, e = t[ok], e[ok]
    err = e - t
    out = dict(n=int(ok.sum()), bias=err.mean(), rmse=np.sqrt(np.mean(err ** 2)), mae=np.abs(err).mean(),
               r2=1 - np.sum(err ** 2) / np.sum((t - t.mean()) ** 2), spearman=spearmanr(t, e).statistic)
    if lo is not None:
        lo, hi = np.asarray(lo, float)[ok], np.asarray(hi, float)[ok]
        out.update(coverage95=float(np.mean((t >= lo) & (t <= hi))), ci_width=float(np.mean(hi - lo)))
    return out


def run(phase1_l2=None, demo=False, biosnicar=None, workers=1, cache_dir=".", tier="C",
        phenol="tddft", photosynthetic=True, spacecraft="S2A", verbose=True):
    import emulator as E
    import inversion as INV
    from priors import PriorConfig, prior_logpdfs

    df = field_band_reflectance(spacecraft)
    R = df[["B2", "B3", "B4", "B8"]].to_numpy()
    pc = PriorConfig()
    szas = sorted(set(int(round(s)) for s in df.sza))
    res = {}
    for model in ("ours", "tierA"):
        ems = {}
        for z in szas:
            cfg = E.EmulatorConfig(model=model, tier=tier, phenol=phenol, photosynthetic=photosynthetic, sza=z,
                                   spacecraft=spacecraft)
            tag = f"{model}_{tier}_{phenol}_{int(photosynthetic)}_{spacecraft}_sza{z}"
            ems[z] = E.build_emulator(cfg, phase1_l2=phase1_l2, demo=demo, biosnicar=biosnicar, workers=workers,
                                      cache=os.path.join(cache_dir, f"field_emulator_{tag}.npz"),
                                      verbose=verbose).refine_log_b(0.05)
        # evidence of every sample under every sigma (each sample uses the emulator at its own SZA)
        logz = np.zeros((len(df), len(SIGMA_GRID)))
        post = {}
        for j, s in enumerate(SIGMA_GRID):
            for z, em in ems.items():
                idx = np.flatnonzero(np.round(df.sza).astype(int).to_numpy() == z)
                lp, sk, mk, _ = prior_logpdfs(em.axes, len(idx), pc)
                r = INV.GridPosterior(em, np.full(4, s)).run(R[idx], lp, sk, mk=mk)
                logz[idx, j] = r["log_evidence"]
                for q in ("log_b_mean", "log_b_q025", "log_b_q975", "log_b_sd", "f_n_mean", "r_um_mean", "chi2",
                          "k_map", "rf_algae_mean", "bba_mean"):
                    post.setdefault((j, q), np.full(len(df), np.nan))[idx] = r[q]
        # leave-one-out choice of sigma: maximise the evidence of the other samples
        loo = {q: np.full(len(df), np.nan) for q in ("log_b_mean", "log_b_q025", "log_b_q975", "log_b_sd",
                                                     "f_n_mean", "r_um_mean", "chi2", "k_map", "rf_algae_mean",
                                                     "bba_mean")}
        sig_loo = np.empty(len(df))
        for i in range(len(df)):
            others = np.delete(np.arange(len(df)), i)
            j = int(np.argmax(logz[others].sum(axis=0)))
            sig_loo[i] = SIGMA_GRID[j]
            for q in loo:
                loo[q][i] = post[(j, q)][i]
        j_all = int(np.argmax(logz.sum(axis=0)))
        res[model] = dict(loo=loo, sigma_loo=sig_loo, sigma_all=SIGMA_GRID[j_all],
                          total_log_evidence=logz[:, j_all].sum(), logz=logz)
        if verbose:
            print(f"[field] {model}: sigma (max marginal likelihood, all samples) = {SIGMA_GRID[j_all]:.3f}; "
                  f"total log evidence {logz[:, j_all].sum():.1f}")

    # empirical band-ratio regression, leave-one-out over samples with cells > 0
    I = (R[:, 2] - R[:, 0]) / (R[:, 2] + R[:, 0])
    pos = df.cells.to_numpy() > 0
    y = np.where(pos, np.log10(np.where(pos, df.cells, 1.0)), np.nan)
    emp = np.full(len(df), np.nan)
    for i in range(len(df)):
        tr = pos.copy()
        tr[i] = False
        X = np.column_stack([np.ones(tr.sum()), I[tr], I[tr] ** 2])
        c, *_ = np.linalg.lstsq(X, y[tr], rcond=None)
        emp[i] = c[0] + c[1] * I[i] + c[2] * I[i] ** 2

    for model in ("ours", "tierA"):
        for q, v in res[model]["loo"].items():
            df[f"{model}_{q}"] = v
        df[f"{model}_sigma"] = res[model]["sigma_loo"]
    df["empirical_log_b"] = emp
    cook = df.cook2020_inversion.to_numpy(float)
    df["cook2020_log_b"] = np.where(cook > 0, np.log10(np.where(cook > 0, cook, 1.0)), np.nan)
    df["log_b_obs"] = y

    rows = []
    sel = pos
    for name, est, lo, hi in (
            ("Physics-informed Bayesian (ours)", df.ours_log_b_mean, df.ours_log_b_q025, df.ours_log_b_q975),
            ("Tier A Bayesian", df.tierA_log_b_mean, df.tierA_log_b_q025, df.tierA_log_b_q975),
            ("Empirical band-ratio (LOO)", df.empirical_log_b, None, None),
            ("Cook et al. 2020 inversion (published)", df.cook2020_log_b, None, None)):
        rows.append(dict(method=name, **_metrics(y[sel], np.asarray(est)[sel],
                                                 None if lo is None else np.asarray(lo)[sel],
                                                 None if hi is None else np.asarray(hi)[sel])))
    metrics = pd.DataFrame(rows)
    n_cook_zero = int(np.sum(sel & ~np.isfinite(df.cook2020_log_b)))
    metrics.loc[metrics.method.str.startswith("Cook"), "note"] = (
        f"{n_cook_zero} of {sel.sum()} counted samples were retrieved as 0 cells and are excluded")
    return df, metrics, res
