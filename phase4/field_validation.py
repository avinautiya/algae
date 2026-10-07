"""
Field validation against co-located measurements: every sample has a nadir HCRF spectrum and a
haemocytometer count of the algae in the ice under the spectrometer (phase2/empirical_data.field_samples):

  s6_2017     Cook et al. (2020) archive, S6 (SW Greenland), July 2017, 46 samples with a spectrum
              (5 with 0 cells; 20_7_SB4's archived spectrum is empty)
  sgris_2021  Chevrollier et al. (2023), southern Greenland ice sheet, Aug 2021, 18 counted samples -
              an independent site, year, team and spectrometer (ASD FieldSpec4)

Spectra are band-averaged with the official ESA Sentinel-2 spectral response functions, so the test
exercises exactly the retrieval applied to Sentinel-2 pixels - minus the atmosphere and the 10-20 m
footprint. Each sample uses the emulator at its own solar zenith angle.

Methods compared (all out of sample):
  physics-informed Bayesian (ours)   emulator + grid posterior; reflectance noise sigma chosen by
                                     maximising the marginal likelihood of the OTHER samples
  Tier A Bayesian                    same machinery with BioSNICAR's empirical algal optics
  empirical band-ratio regression    log10 B = c0 + c1 I + c2 I^2, I = (B4-B2)/(B4+B2), leave-one-out
  Cook et al. (2020) inversion       the published BioSNICAR-GO retrieval (S6 samples only), as reported

Ice-radius prior: the measured specific surface area of bubbly ice (Cooper et al. 2021; Dadic et al.
2013) characterises ice 0.1-1 m deep, not the weathering crust at the surface. The population
distribution of the surface radius is therefore estimated from the field spectra themselves by
empirical Bayes (log-normal hyper-parameters maximising the marginal likelihood of the other samples,
jointly with sigma), and compared with the measured-SSA prior by their marginal likelihoods.

Observation error: counts have Poisson error 1/sqrt(N counted) (S6; N is in the count workbook);
'coverage95_obs' scores the 95 % interval against the observation widened by that error.
Zero counts (below one counted cell) cannot be scored on a log scale and are reported separately.
"""

from __future__ import annotations

import dataclasses
import hashlib
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "phase2"))
sys.path.insert(0, HERE)

import empirical_data as ED  # noqa: E402

SIGMA_GRID = np.array([0.01, 0.015, 0.02, 0.03, 0.04, 0.06, 0.08])
# empirical-Bayes grid for the population distribution of the (surface) ice radius: ln r ~ N(mu, sd)
EB_MU_LNR = np.log([600.0, 900.0, 1350.0, 2000.0, 3000.0, 4500.0])
EB_SD_LNR = np.array([0.35, 0.6, 0.9, 1.2, 1.6, 2.0, 2.5])
DATASET_LABEL = {"s6_2017": "S6 2017 (Cook et al. 2020)", "sgris_2021": "S Greenland 2021 (Chevrollier et al. 2023)",
                 "s6_2014": "S6 2014 (Stibal et al. 2017; manual ingestion)"}
CACHE_VERSION = "v3"


def field_band_reflectance(spacecraft="S2A", bands=("B2", "B3", "B4", "B8")):
    tab, spectra = ED.field_samples()
    srf = pd.read_csv(ED.path("esa_s2_srf_TN-15-0007_v4.0.csv")).set_index("wavelength_nm")
    rows = []
    for r in tab.itertuples():
        h = spectra[r.dataset]
        wl = h.wavelength_nm.to_numpy()
        sp = h[r.sample].to_numpy(float)
        vals = []
        ok = np.isfinite(sp)
        if not ok.any():               # 20_7_SB4: the archived HCRF column is empty
            continue
        for b in bands:
            s = srf[f"{spacecraft}_{b}"].reindex(wl).fillna(0).to_numpy()
            vals.append(np.sum(sp[ok] * s[ok]) / s[ok].sum())
        rows.append(dict(dataset=r.dataset, sample=r.sample, cells=r.cells, cells_counted=r.cells_counted,
                         sza=r.sza, quantity=r.quantity, **dict(zip(bands, vals))))
    out = pd.DataFrame(rows)
    out.attrs["albedo_k_sd"] = tab.attrs.get("albedo_k_sd", np.nan)
    return out


def _metrics(t, e, lo=None, hi=None, sd=None, obs_sd=None):
    from scipy.stats import spearmanr
    t, e = np.asarray(t, float), np.asarray(e, float)
    ok = np.isfinite(t) & np.isfinite(e)
    err = e[ok] - t[ok]
    tt = t[ok]
    out = dict(n=int(ok.sum()), bias=err.mean(), rmse=np.sqrt(np.mean(err ** 2)), mae=np.abs(err).mean(),
               r2=1 - np.sum(err ** 2) / np.sum((tt - tt.mean()) ** 2), spearman=spearmanr(tt, e[ok]).statistic)
    if lo is not None:
        lo, hi = np.asarray(lo, float)[ok], np.asarray(hi, float)[ok]
        out.update(coverage95=float(np.mean((tt >= lo) & (tt <= hi))), ci_width=float(np.mean(hi - lo)))
    if sd is not None:
        o = np.nan_to_num(np.asarray(obs_sd, float)[ok], nan=0.0)
        z = err / np.sqrt(np.asarray(sd, float)[ok] ** 2 + o ** 2)
        out["coverage95_obs"] = float(np.mean(np.abs(z) <= 1.96))
    return out


def run(phase1_l2=None, demo=False, biosnicar=None, workers=1, cache_dir=".", tier="C",
        phenol="tddft", photosynthetic=True, spacecraft="S2A", rho_bottom=690.0, verbose=True,
        emu_overrides: dict | None = None, models=("ours", "tierA"), dust: bool = True):
    """emu_overrides: extra EmulatorConfig fields for both models (e.g. dust_ppb=E.DUST_NODES_PPB,
    film_dz=0.002, f_n=(0, 1, 0.2)) - used by bias_study.py."""
    import emulator as E
    import inversion as INV
    from priors import PriorConfig, prior_logpdfs

    df = field_band_reflectance(spacecraft)
    R = df[["B2", "B3", "B4", "B8"]].to_numpy()
    pc = PriorConfig.for_density(rho_bottom)
    zs = np.round(df.sza).astype(int).to_numpy()
    res = {}
    ov = dict(dict(dust_ppb=E.DUST_NODES_PPB) if dust else {}, **(emu_overrides or {}))
    ov_tag = "_".join(f"{k}{v}" for k, v in sorted(ov.items())).replace(" ", "")
    for model in models:
        ems = {}
        for z in sorted(set(zs)):
            cfg = E.EmulatorConfig(model=model, tier=tier, phenol=phenol, photosynthetic=photosynthetic, sza=z,
                                   spacecraft=spacecraft, rho_bottom=rho_bottom, **ov)
            tag = (f"{CACHE_VERSION}_{model}_{tier}_{phenol}_{int(photosynthetic)}_{spacecraft}_rb{rho_bottom:.0f}"
                   f"_sza{z}_{hashlib.md5(ov_tag.encode()).hexdigest()[:8]}")
            ems[z] = E.build_emulator(cfg, phase1_l2=phase1_l2, demo=demo, biosnicar=biosnicar, workers=workers,
                                      cache=os.path.join(cache_dir, f"field_emulator_{tag}.npz"),
                                      verbose=verbose).refine_log_b(0.05)
        qs = ("log_b_mean", "log_b_q025", "log_b_q975", "log_b_sd", "f_n_mean", "r_um_mean", "chi2", "k_map",
              "rf_algae_mean", "bba_mean")
        # hyper-parameter candidates: reflectance noise sigma x ice-radius prior (measured-SSA prior, or a
        # log-normal population distribution estimated from the field spectra by empirical Bayes)
        hyp = [(s_, "measured_ssa", pc.mu_lnr, pc.sd_lnr) for s_ in SIGMA_GRID]
        hyp += [(s_, "empirical_bayes", m_, d_) for s_ in SIGMA_GRID for m_ in EB_MU_LNR for d_ in EB_SD_LNR]
        logz = np.zeros((len(df), len(hyp)))
        post = {}
        for j, (s_, _, m_, d_) in enumerate(hyp):
            cfg_j = dataclasses.replace(pc, mu_lnr=m_, sd_lnr=d_)
            for z, em in ems.items():
                for qty in sorted(set(df.quantity)):
                    idx = np.flatnonzero((zs == z) & (df.quantity == qty).to_numpy())
                    if idx.size == 0:
                        continue
                    lp, sk, mk, _ = prior_logpdfs(em.axes, len(idx), cfg_j)
                    if qty == "albedo":       # albedo spectra: no HCRF/albedo anisotropy, k ~ N(1, sd)
                        sk, mk = df.attrs["albedo_k_sd"], 1.0
                    r = INV.GridPosterior(em, np.full(4, s_)).run(R[idx], lp, sk, mk=mk)
                    logz[idx, j] = r["log_evidence"]
                    for q in qs:
                        post.setdefault((j, q), np.full(len(df), np.nan))[idx] = r[q]
        # leave-one-out choice of all hyper-parameters: maximise the evidence of the other samples
        loo = {q: np.full(len(df), np.nan) for q in qs}
        sig_loo = np.empty(len(df))
        for i in range(len(df)):
            j = int(np.argmax(np.delete(logz, i, axis=0).sum(axis=0)))
            sig_loo[i] = hyp[j][0]
            for q in loo:
                loo[q][i] = post[(j, q)][i]
        tot = logz.sum(axis=0)
        j_all = int(np.argmax(tot))
        j_meas = int(np.argmax(np.where([h[1] == "measured_ssa" for h in hyp], tot, -np.inf)))
        res[model] = dict(loo=loo, sigma_loo=sig_loo, sigma_all=hyp[j_all][0], r_prior=hyp[j_all][1],
                          mu_lnr=hyp[j_all][2], sd_lnr=hyp[j_all][3], total_log_evidence=tot[j_all],
                          log_evidence_measured_ssa_prior=tot[j_meas], sigma_measured_ssa_prior=hyp[j_meas][0],
                          logz=logz, hyper=hyp)
        if verbose:
            print(f"[field] {model}: sigma {hyp[j_all][0]:.3f}, radius prior {hyp[j_all][1]} (median "
                  f"{np.exp(hyp[j_all][2]):.0f} um, ln-SD {hyp[j_all][3]:.2f}); log evidence {tot[j_all]:.1f} vs "
                  f"{tot[j_meas]:.1f} with the measured-SSA prior")

    # empirical band-ratio regression, leave-one-out over counted samples
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

    for model in models:
        for q, v in res[model]["loo"].items():
            df[f"{model}_{q}"] = v
        df[f"{model}_sigma"] = res[model]["sigma_loo"]
    df["empirical_log_b"] = emp
    pub = ED.cook2020_published_inversion()
    cook = df["sample"].map(pub).astype(float).where(df.dataset == "s6_2017").to_numpy()
    df["cook2020_log_b"] = np.where(cook > 0, np.log10(np.where(cook > 0, cook, 1.0)), np.nan)
    df["log_b_obs"] = y
    df["log_b_obs_sd"] = np.where(pos, 1.0 / np.log(10) / np.sqrt(df.cells_counted.to_numpy()), np.nan)

    rows = []
    for ds, sel in [("all", pos)] + [(d, pos & (df.dataset == d).to_numpy()) for d in DATASET_LABEL]:
        bayes = [(lab, df[f"{m}_log_b_mean"], df[f"{m}_log_b_q025"], df[f"{m}_log_b_q975"], df[f"{m}_log_b_sd"])
                 for m, lab in (("ours", "Physics-informed Bayesian (ours)"), ("tierA", "Tier A Bayesian"))
                 if m in models]
        for name, est, lo, hi, sd in bayes + [
                ("Empirical band-ratio (LOO)", df.empirical_log_b, None, None, None),
                ("Cook et al. 2020 inversion (published)", df.cook2020_log_b, None, None, None)]:
            e = np.asarray(est)[sel]
            if np.isfinite(e).sum() < 3:
                continue
            rows.append(dict(dataset=ds, method=name, **_metrics(
                y[sel], e, None if lo is None else np.asarray(lo)[sel], None if hi is None else np.asarray(hi)[sel],
                None if sd is None else np.asarray(sd)[sel], df.log_b_obs_sd.to_numpy()[sel])))
    metrics = pd.DataFrame(rows)
    s6 = pos & (df.dataset == "s6_2017").to_numpy()
    n_cook_zero = int(np.sum(s6 & ~np.isfinite(df.cook2020_log_b)))
    metrics["note"] = ""
    metrics.loc[metrics.method.str.startswith("Cook"), "note"] = (
        f"published values exist for 31 S6 samples; {n_cook_zero} of {s6.sum()} counted S6 samples have no "
        "published value or were retrieved as 0 cells and are excluded")
    return df, metrics, res
