#!/usr/bin/env python
"""
Held-out (leave-site-out) predictive evaluation: does pigment-aware optics predict an UNSEEN site better
than empirical algal optics and statistical reflectance baselines?

Data: co-located nadir HCRF spectra + algal counts (phase2/empirical_data.field_samples)
  s6_2017    S6, SW Greenland, Cook et al. (2020): 46 samples (5 zero counts), counted cells N and
             counted volume V = N / (cells mL^-1) known -> Poisson count likelihood
  sgris_2021 southern Greenland, Chevrollier et al. (2023): 18 samples, number counted not published ->
             log-normal observation error with assumed SD (OBS_SD_SGRIS, sensitivity runs)

Folds (everything fitted on the training site only, nothing tuned on the test site):
  primary    train S6 2017      -> test S Greenland 2021
  secondary  train S Greenland  -> test S6 2017

Inverse test (abundance from reflectance), per frozen optical model:
  * the forward model (emulator) is fixed before any field data are seen (BioSNICAR + optics; the TD-DFT
    spectrum is calibrated on laboratory spectra only);
  * hyper-parameters (reflectance noise sigma, ice-radius population prior) by empirical Bayes on the
    TRAINING site (sum of training log evidences);
  * structural discrepancy tau (dex) by maximum likelihood of the TRAINING observations, using inner
    leave-one-out posteriors (hyper-parameters refitted without the sample) so tau is not fitted on
    in-sample posteriors;
  * test predictive density of log10 B = posterior marginal (grid) convolved with N(0, tau^2);
  * scored with the count likelihood (zeros kept: P(N = 0 | B) = exp(-B V)).
Baselines (same training data, same scoring):
  climatology   N(mu, s) of log10 B fitted to the training counts
  band ratio    log10 B ~ N(c0 + c1 I + c2 I^2, s), I = (B4 - B2)/(B4 + B2)
  ridge         log10 B ~ N(linear in standardised reflectance features, s), lambda by training LOO
  literature prior  N(3.56, 0.78) (S6 2016, Williamson et al. 2020; not fitted here)
  Cook et al. 2020 published inversion (S6 only; point values, not trained here)

Forward test (nadir HCRF from the MEASURED abundance; positive counts; field spectra are HCRFs, not albedos):
  log p(4-band HCRF | log10 B_obs) under each optical model with the training hyper-parameters (k, ice
  radius, community and dust priors marginalised) vs. a training-site Gaussian climatology of the bands and
  a training linear regression of the bands on log10 B. Broadband: prior-predictive k x BBA (BBA = model
  broadband ALBEDO, 300-2500 nm; k = HCRF/albedo anisotropy factor) vs the measured broadband HCRF on the
  same 300-2500 nm support (300-350 nm extrapolated from 350 nm; its irradiance share is reported).

Metrics (summary_stats): every metric refers to the observed log10 count through the predictive
distribution of the OBSERVATION; the log score is reported on all samples (zeros included, exact count
likelihood) and on the positive-count subset on which RMSE, bias, coverage and CRPS are computed. The
observation SD for S Greenland is the requested --obs-sd-sgris everywhere.
Uncertainty of paired contrasts: cluster bootstrap over sampling days (few blocks are flagged).
Pre-registration: docs/preregistration_heldout.md (primary model, contrast, decision rule).

Outputs: per-sample scores, per-fold summary tables and paired differences with bootstrap 95 % CIs.

    python phase4/heldout.py --outdir phase4/results/heldout
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import os
import sys
import time

import numpy as np
import pandas as pd
from scipy.special import logsumexp
from scipy.stats import norm

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "phase2"))

import empirical_data as ED  # noqa: E402
import field_validation as FV  # noqa: E402

FOLDS = (("primary", "s6_2017", "sgris_2021"), ("secondary", "sgris_2021", "s6_2017"))
OPTICS = {
    "tierA_empirical": dict(model="tierA"),
    "measured_mac_C": dict(model="ours", tier="C", phenol="williamson2020"),
    "tddft_C": dict(model="ours", tier="C", phenol="tddft"),
    "tddft_D": dict(model="ours", tier="D", phenol="tddft"),
    # calibration structural scenario: independent-residual likelihood (phi 0.49) instead of AR(1) (phi 0.11)
    "tddft_D_iid": dict(model="ours", tier="D", phenol="tddft", cal_residuals="iid"),
}
PRIMARY = "tddft_D"             # pre-registered (docs/preregistration_heldout.md)
PRIMARY_CONTRAST = "tierA_empirical"
OBS_SD_SGRIS = 0.10            # dex, assumed (counts in quintuplicate, number counted not published)
V_ZERO_ML = 0.016              # counted volume assumed for zero counts (largest volume in the S6 workbook)
XGRID = np.arange(-1.0, 7.0 + 1e-9, 0.01)
DX = XGRID[1] - XGRID[0]
BANDS = ["B2", "B3", "B4", "B8"]


# --------------------------------------------------------------------------- observations
def observation_loglik(df, obs_sd_sgris=OBS_SD_SGRIS, v_zero=V_ZERO_ML):
    """L[i, g] = log p(observation_i | log10 B = XGRID[g])."""
    L = np.zeros((len(df), XGRID.size))
    B = 10.0 ** XGRID
    for i, r in enumerate(df.itertuples()):
        if np.isfinite(r.cells_counted):
            N = float(r.cells_counted)
            V = N / r.cells if r.cells > 0 else v_zero
            lam = B * V
            L[i] = N * np.log(np.maximum(lam, 1e-300)) - lam - (np.sum(np.log(np.arange(1, N + 1))) if N > 0 else 0.0)
        else:
            if not r.cells > 0:
                raise ValueError(f"{r.sample}: zero count without counted volume")
            L[i] = norm.logpdf(np.log10(r.cells), XGRID, obs_sd_sgris)
    return L


def log_score(logp_x, L):
    """log of the predictive probability of each observation: log sum_g p(x_g) dx exp(L_g)."""
    return logsumexp(logp_x + L, axis=1) + np.log(DX)


def gaussian_logp(mu, s):
    mu, s = np.atleast_1d(mu)[:, None], np.atleast_1d(s)[:, None]
    return norm.logpdf(XGRID[None, :], mu, s)


def observation_sd(df, obs_sd_sgris=OBS_SD_SGRIS):
    """SD (dex) of the observed log10 count: Poisson 1/(ln10 sqrt(N)) where the counted number is known,
    otherwise the assumed (requested) SD. Zero counts have no log value (NaN)."""
    return np.where(np.isfinite(df.cells_counted), 1 / np.log(10) / np.sqrt(np.maximum(df.cells_counted, 1)),
                    obs_sd_sgris)


def summary_stats(logp_x, df, L, obs_sd_sgris=OBS_SD_SGRIS):
    """Metrics that all refer to the SAME target - the observed log10 count - through the predictive
    distribution of the OBSERVATION (latent predictive convolved with the observation error, Gaussian in
    log10 with observation_sd; the log score uses the exact count likelihood). Point prediction, 95 %
    interval coverage, CRPS and the positive-sample log score are computed on the same positive-count set.
    Returns dict of per-sample arrays."""
    from scipy.ndimage import gaussian_filter1d
    p = np.exp(logp_x - logsumexp(logp_x, axis=1, keepdims=True))
    y = np.where(df.cells > 0, np.log10(np.where(df.cells > 0, df.cells, 1.0)), np.nan)
    osd = observation_sd(df, obs_sd_sgris)
    n = len(df)
    med, cover, crps = (np.full(n, np.nan) for _ in range(3))
    for i in range(n):
        if not np.isfinite(y[i]):
            continue
        q = gaussian_filter1d(p[i], osd[i] / DX, mode="constant", cval=0.0, truncate=6.0)
        q = q / q.sum()
        c = np.cumsum(q)
        med[i] = XGRID[np.argmax(c >= 0.5)]
        lo, hi = XGRID[np.argmax(c >= 0.025)], XGRID[np.argmax(c >= 0.975)]
        cover[i] = float(lo <= y[i] <= hi)
        crps[i] = np.sum((c - (XGRID >= y[i])) ** 2) * DX
    return dict(median_obs=med, y=y, covered=cover, crps=crps, log_score=log_score(logp_x, L))


def sampling_block(df):
    """Sampling-day block of each sample: S6 '13_7_SB5' -> '13_7'; S Greenland '210805-S1' -> '210805'."""
    s = df["sample"].astype(str)
    return np.where(df.dataset == "s6_2017", s.str.rsplit("_", n=1).str[0], s.str.split("-").str[0])


# --------------------------------------------------------------------------- physics models
def physics_tables(df, optics_key, phase1_l2, biosnicar, workers, cache_dir, rho_bottom=690.0, verbose=True):
    """For every sample and hyper-parameter candidate h: log evidence, posterior marginal of log10 B on
    the emulator nodes, and the forward log p(R | log10 B node) for the radius/sigma candidate.
    Returns dict(logz (n, H), marg (n, H, nb), nodes (nb,), hyp list, fwd (n, H), bb_pred (n, H, 2))."""
    import emulator as E
    import inversion as INV
    from priors import PriorConfig, prior_logpdfs
    R = df[BANDS].to_numpy()
    zs = np.round(df.sza).astype(int).to_numpy()
    pc0 = PriorConfig.for_density(rho_bottom)
    hyp = [(s_, "measured_ssa", pc0.mu_lnr, pc0.sd_lnr) for s_ in FV.SIGMA_GRID]
    hyp += [(s_, "empirical_bayes", m_, d_) for s_ in FV.SIGMA_GRID for m_ in FV.EB_MU_LNR for d_ in FV.EB_SD_LNR]
    ems = {}
    for z in sorted(set(zs)):
        cfg = E.EmulatorConfig(sza=int(z), spacecraft="S2A", rho_bottom=rho_bottom, dust_ppb=E.DUST_NODES_PPB,
                               photosynthetic=True, **OPTICS[optics_key])
        ems[z] = E.build_emulator(cfg, phase1_l2=phase1_l2, biosnicar=biosnicar, workers=workers,
                                  cache=os.path.join(cache_dir, f"heldout_{optics_key}_sza{z}.npz"),
                                  verbose=verbose).refine_log_b(0.05)
    nodes = next(iter(ems.values())).axes["log_b"]
    n, H = len(df), len(hyp)
    logz = np.full((n, H), np.nan)
    marg = np.full((n, H, nodes.size), np.nan)
    fwd = np.full((n, H), np.nan)
    bb = np.full((n, H, 2), np.nan)
    y = np.where(df.cells > 0, np.log10(np.where(df.cells > 0, df.cells, 1.0)), np.nan)
    sigmas = sorted({h[0] for h in hyp})
    for z, em in ems.items():
        idx = np.flatnonzero(zs == z)
        names = em.names
        ir = names.index("r_um")
        lp0, sk, mk, _ = prior_logpdfs(em.axes, len(idx), pc0)
        # prior weights over (f_n, dust) for the broadband prior-predictive (do not depend on the hyper-params)
        other = [a_ for a_ in names if a_ not in ("log_b", "r_um")]
        w_fd = np.exp(sum(lp0[a_][0].reshape([-1 if b_ == a_ else 1 for b_ in names]) for a_ in other))
        bba = em.data["bba"]
        sum_fd = tuple(i for i, a_ in enumerate(names) if a_ in other)
        e1 = np.sum(w_fd * bba, axis=sum_fd) / np.sum(w_fd)                    # (nb, nr)
        e2 = np.sum(w_fd * bba ** 2, axis=sum_fd) / np.sum(w_fd)
        pos_local = np.flatnonzero(np.isfinite(y[idx]))
        inode = np.abs(nodes[None, :] - y[idx][pos_local][:, None]).argmin(axis=1)
        for s_ in sigmas:
            gp = INV.GridPosterior(em, np.full(4, s_), derived=())
            gp.half_logden = 0.5 * np.log(1.0 + sk ** 2 * gp.a)
            ll, _, _ = gp._loglik(R[idx] / s_, sk, mk)
            ll = ll.reshape(len(idx), *em.shape)
            for a_ in other:                                                   # add f_n and dust priors
                shp = [1] + [1] * len(names)
                shp[names.index(a_) + 1] = -1
                ll = ll + lp0[a_][0].reshape(shp)
            B0 = logsumexp(ll, axis=tuple(i + 1 for i, a_ in enumerate(names) if a_ in other))   # (p, nb, nr)
            lb = lp0["log_b"][0][None, :, None]
            for j, (s_h, _, m_, d_) in enumerate(hyp):
                if s_h != s_:
                    continue
                lr = prior_logpdfs(em.axes, 1, dataclasses.replace(pc0, mu_lnr=m_, sd_lnr=d_))[0]["r_um"][0]
                post = B0 + lb + lr[None, None, :]
                lz = logsumexp(post, axis=(1, 2))
                logz[idx, j] = lz + gp.norm_const
                mb = logsumexp(post, axis=2)
                marg[idx, j] = np.exp(mb - lz[:, None])
                if pos_local.size:
                    fwd[idx[pos_local], j] = logsumexp(B0[pos_local, inode] + lr[None, :], axis=1) + gp.norm_const
                    wr = np.exp(lr - logsumexp(lr))
                    mbb, m2 = e1[inode] @ wr, e2[inode] @ wr
                    vb = m2 - mbb ** 2
                    bb[idx[pos_local], j, 0] = mk * mbb
                    bb[idx[pos_local], j, 1] = np.sqrt(sk ** 2 * (vb + mbb ** 2) + mk ** 2 * vb)
        if verbose:
            print(f"[heldout] {optics_key}: SZA {z} done", flush=True)
    return dict(logz=logz, marg=marg, nodes=nodes, hyp=hyp, fwd=fwd, bb=bb)


def table_code_fingerprint():
    """Fingerprint of the code that builds the physics tables (physics_tables and its helpers only), so
    that edits to the scoring code do not invalidate expensive cached tables."""
    import hashlib
    import inspect
    src = "".join(inspect.getsource(f) for f in (physics_tables, marg_to_logp))
    return hashlib.sha256(src.encode()).hexdigest()


def marg_to_logp(marg_row, nodes, tau):
    """Posterior marginal (mass per node) -> log density on XGRID, convolved with N(0, tau^2)."""
    step = nodes[1] - nodes[0]
    dens = np.interp(XGRID, nodes, marg_row / step, left=0.0, right=0.0)
    if tau > 0:
        from scipy.ndimage import gaussian_filter1d
        dens = gaussian_filter1d(dens, tau / DX, mode="constant", cval=0.0, truncate=6.0)
    dens = dens / max(dens.sum() * DX, 1e-300)
    with np.errstate(divide="ignore"):
        return np.log(np.maximum(dens, 1e-300))


TAU_GRID = np.round(np.arange(0.0, 1.5001, 0.02), 3)


def physics_fold(tab, train, test, L):
    """Training-only hyper-parameters and tau; test log scores and densities."""
    tr, te = np.flatnonzero(train), np.flatnonzero(test)
    tot = np.nansum(tab["logz"][tr], axis=0)
    h = int(np.argmax(tot))
    # inner LOO: hyper-parameters refitted without each training sample, then tau by ML on those posteriors
    lp_loo = {}
    for i in tr:
        hi = int(np.argmax(tot - np.nan_to_num(tab["logz"][i])))
        lp_loo[i] = tab["marg"][i, hi]
    ll_tau = np.array([sum(log_score(marg_to_logp(lp_loo[i], tab["nodes"], t)[None, :], L[i:i + 1])[0] for i in tr)
                       for t in TAU_GRID])
    tau = float(TAU_GRID[int(np.argmax(ll_tau))])
    logp = np.array([marg_to_logp(tab["marg"][i, h], tab["nodes"], tau) for i in te])
    return dict(h=h, hyper=tab["hyp"][h], tau=tau, logp=logp, fwd=tab["fwd"][te, h], bb=tab["bb"][te, h])


# --------------------------------------------------------------------------- statistical baselines
def features(df):
    R = df[BANDS].to_numpy()
    I = (R[:, 2] - R[:, 0]) / (R[:, 2] + R[:, 0])
    return R, I


def fit_gaussian_model(X, L, l2=0.0):
    """log10 B ~ N(X beta, s); maximise sum_i log int N(x; X_i beta, s) exp(L_i(x)) dx - l2 |beta[1:]|^2."""
    from scipy.optimize import minimize
    y0 = XGRID[np.argmax(L, axis=1)]
    beta0, *_ = np.linalg.lstsq(X, y0, rcond=None)

    def nll(t):
        beta, s = t[:-1], np.exp(t[-1])
        return -np.sum(log_score(gaussian_logp(X @ beta, np.full(len(X), s)), L)) + l2 * np.sum(beta[1:] ** 2)
    r = minimize(nll, np.r_[beta0, np.log(0.5)], method="L-BFGS-B")
    return r.x[:-1], float(np.exp(r.x[-1])), bool(r.success)


def baselines(df, train, test, L):
    out = {}
    R, I = features(df)
    tr, te = np.flatnonzero(train), np.flatnonzero(test)
    # climatology
    b, s, ok = fit_gaussian_model(np.ones((tr.size, 1)), L[tr])
    out["climatology"] = dict(logp=gaussian_logp(np.full(te.size, b[0]), np.full(te.size, s)), fit_ok=ok)
    # literature prior (independent S6 2016 data)
    mu_b, sd_b, _, _ = ED.abundance_prior()
    out["literature_prior"] = dict(logp=gaussian_logp(np.full(te.size, mu_b), np.full(te.size, sd_b)), fit_ok=True)
    # band ratio
    Xq = np.column_stack([np.ones(len(df)), I, I ** 2])
    b, s, ok = fit_gaussian_model(Xq[tr], L[tr])
    out["band_ratio"] = dict(logp=gaussian_logp(Xq[te] @ b, np.full(te.size, s)), fit_ok=ok)
    # ridge on standardised features; lambda by training leave-one-out log score
    F = np.column_stack([R, np.log(R), R[:, 1] / R[:, 0], R[:, 2] / R[:, 1], R[:, 3] / R[:, 2], I])

    def design(rows_fit, rows_apply):
        """Standardisation fitted on rows_fit ONLY (inside each inner fold), applied to rows_apply."""
        mu, sd = F[rows_fit].mean(0), F[rows_fit].std(0) + 1e-12
        return (np.column_stack([np.ones(len(rows_fit)), (F[rows_fit] - mu) / sd]),
                np.column_stack([np.ones(len(rows_apply)), (F[rows_apply] - mu) / sd]))
    best = None
    for lam in (0.1, 1.0, 3.0, 10.0, 30.0, 100.0):
        sc = 0.0
        for k in range(tr.size):
            inner, left = np.delete(tr, k), tr[k:k + 1]
            Xf, Xa = design(inner, left)
            b, s, _ = fit_gaussian_model(Xf, L[inner], lam)
            sc += log_score(gaussian_logp(Xa @ b, np.array([s])), L[left])[0]
        if best is None or sc > best[0]:
            best = (sc, lam)
    Xf, Xa = design(tr, te)
    b, s, ok = fit_gaussian_model(Xf, L[tr], best[1])
    out["ridge"] = dict(logp=gaussian_logp(Xa @ b, np.full(te.size, s)), fit_ok=ok, lam=best[1])
    return out


def forward_baselines(df, train, test):
    """Gaussian log density of the 4-band HCRF: training climatology and training regression on log10 B."""
    from scipy.stats import multivariate_normal as mvn
    R = df[BANDS].to_numpy()
    y = np.where(df.cells > 0, np.log10(np.where(df.cells > 0, df.cells, 1.0)), np.nan)
    tr = np.flatnonzero(train & np.isfinite(y))
    te = np.flatnonzero(test)
    out = {}
    C = np.cov(R[tr].T) + 1e-8 * np.eye(4)
    out["climatology"] = np.array([mvn(R[tr].mean(0), C).logpdf(R[i]) if np.isfinite(y[i]) else np.nan for i in te])
    X = np.column_stack([np.ones(tr.size), y[tr]])
    B, *_ = np.linalg.lstsq(X, R[tr], rcond=None)
    res = R[tr] - X @ B
    Cr = np.cov(res.T) + 1e-8 * np.eye(4)
    out["regression_on_logB"] = np.array([mvn(np.array([1.0, y[i]]) @ B, Cr).logpdf(R[i]) if np.isfinite(y[i])
                                          else np.nan for i in te])
    return out


def measured_broadband(df):
    """Broadband nadir HCRF of each field spectrum, weighted by the clear-sky irradiance at the sample's SZA,
    on the SAME 300-2500 nm support as the emulator's broadband albedo (BBA). The spectra start at 350 nm;
    300-350 nm (a small share of the irradiance, returned per sample) takes the 350 nm value - a stated
    extrapolation. This is an HCRF, not an albedo: predictions are k x BBA, with k the HCRF/albedo
    anisotropy factor of the retrieval prior. Returns (bb_hcrf, irradiance share 300-350 nm)."""
    import biosnicar_bridge as bb
    tab, spectra = ED.field_samples()
    run = bb.BioSNICARRunner(bb.locate_biosnicar(None))
    wl480 = run.wvl_um * 1000.0
    out, share = np.full(len(df), np.nan), np.full(len(df), np.nan)
    grid = np.arange(300.0, 2500.0 + 0.5, 1.0)
    for i, r in enumerate(df.itertuples()):
        sp = spectra[r.dataset]
        h = sp[r.sample].to_numpy(float)
        wl = sp.wavelength_nm.to_numpy(float)
        ok = np.isfinite(h)
        hg = np.interp(grid, wl[ok], h[ok])                      # below 350 nm: held at the 350 nm value
        flx = np.interp(grid, wl480, run.illumination(r.sza).flx_slr)
        out[i] = np.sum(hg * flx) / np.sum(flx)
        share[i] = flx[grid < 350].sum() / flx.sum()
    return out, share


def block_bootstrap(a, b, blocks, n=4000, seed=0, min_blocks=5):
    """Paired difference mean(a - b) with a cluster bootstrap over sampling-day blocks (samples from one
    day are not independent). Flags when there are too few blocks for the interval to mean much."""
    d = np.asarray(a, float) - np.asarray(b, float)
    blocks = np.asarray(blocks)
    ok = np.isfinite(d)
    d, blocks = d[ok], blocks[ok]
    ub = np.unique(blocks)
    if d.size < 2:
        return dict(mean=np.nan, lo=np.nan, hi=np.nan, n=int(d.size), n_blocks=int(ub.size), few_blocks=True)
    rng = np.random.default_rng(seed)
    idx = {u: np.flatnonzero(blocks == u) for u in ub}
    bs = np.empty(n)
    for k in range(n):
        pick = rng.choice(ub, size=ub.size, replace=True)
        rows = np.concatenate([idx[u] for u in pick])
        bs[k] = d[rows].mean()
    return dict(mean=float(d.mean()), lo=float(np.quantile(bs, 0.025)), hi=float(np.quantile(bs, 0.975)),
                n=int(d.size), n_blocks=int(ub.size), few_blocks=bool(ub.size < min_blocks))


# --------------------------------------------------------------------------- evaluation
def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--phase1-l2", default=os.path.join(HERE, "..", "phase1", "results", "level2"))
    p.add_argument("--biosnicar", default=None)
    p.add_argument("--workers", type=int, default=2)
    p.add_argument("--optics", nargs="+", default=list(OPTICS))
    p.add_argument("--obs-sd-sgris", type=float, default=OBS_SD_SGRIS)
    p.add_argument("--v-zero-ml", type=float, default=V_ZERO_ML)
    p.add_argument("--outdir", default=os.path.join(HERE, "results", "heldout"))
    a = p.parse_args(argv)
    os.makedirs(a.outdir, exist_ok=True)
    cache = os.path.join(a.outdir, "cache")
    os.makedirs(cache, exist_ok=True)
    t0 = time.time()
    df = FV.field_band_reflectance("S2A")
    df = df[df.dataset.isin(["s6_2017", "sgris_2021"])].reset_index(drop=True)
    L = observation_loglik(df, a.obs_sd_sgris, a.v_zero_ml)
    bb_meas, bb_share = measured_broadband(df)
    blocks = sampling_block(df)
    tables = {}
    for k in a.optics:
        path = os.path.join(cache, f"tables_{k}.npz")
        import provenance as PV
        fp = PV.fingerprint(dict(k=k, df=df[BANDS + ["cells", "sza"]], code=PV.code_fingerprint(),
                                 phase1=PV.phase1_fingerprint(a.phase1_l2), data=PV.empirical_data_fingerprint(),
                                 inv=PV.code_fingerprint(("phase4/inversion.py", "phase4/priors.py")),
                                 tables_code=table_code_fingerprint()))
        ok, _ = PV.cached_ok(path, fp)
        if ok:
            z = np.load(path, allow_pickle=True)
            tables[k] = {f: z[f] for f in z.files}
            tables[k]["hyp"] = [tuple(h) for h in tables[k]["hyp"]]
        else:
            tables[k] = physics_tables(df, k, a.phase1_l2, a.biosnicar, a.workers, cache)
            tmp = path + ".tmp.npz"
            np.savez(tmp, **{f: (np.array(v, dtype=object) if f == "hyp" else v) for f, v in tables[k].items()})
            os.replace(tmp, path)
            PV.write_meta(path, fp)
    per_sample, summary, paired, fwd_rows, settings = [], [], [], [], {}
    from scipy.stats import spearmanr
    for fold, trn, tst in FOLDS:
        train = (df.dataset == trn).to_numpy()
        test = (df.dataset == tst).to_numpy()
        te = np.flatnonzero(test)
        methods = {}
        for k in a.optics:
            r = physics_fold(tables[k], train, test, L)
            methods[k] = r["logp"]
            settings[f"{fold}/{k}"] = dict(sigma=r["hyper"][0], radius_prior=r["hyper"][1],
                                           radius_median_um=float(np.exp(r["hyper"][2])), radius_ln_sd=r["hyper"][3],
                                           tau_dex=r["tau"])
            for jj, i in enumerate(te):
                fwd_rows.append(dict(fold=fold, sample=df["sample"][i], block=blocks[i], method=k,
                                     logp_hcrf_bands=r["fwd"][jj], bb_hcrf_pred=r["bb"][jj, 0],
                                     bb_hcrf_pred_sd=r["bb"][jj, 1], bb_hcrf_meas=bb_meas[i],
                                     irradiance_share_300_350=bb_share[i]))
        for k, v in baselines(df, train, test, L).items():
            methods[k] = v["logp"]
            settings[f"{fold}/{k}"] = {kk: vv for kk, vv in v.items() if kk != "logp"}
        for k, v in forward_baselines(df, train, test).items():
            for jj, i in enumerate(te):
                fwd_rows.append(dict(fold=fold, sample=df["sample"][i], block=blocks[i], method=k, logp_hcrf_bands=v[jj],
                                     bb_hcrf_meas=bb_meas[i], bb_hcrf_pred=np.nan, bb_hcrf_pred_sd=np.nan,
                                     irradiance_share_300_350=bb_share[i]))
        Lte = L[te]
        stats_by = {}
        for k, lp in methods.items():
            st_ = summary_stats(lp, df.iloc[te], Lte, a.obs_sd_sgris)
            stats_by[k] = st_
            pos = np.isfinite(st_["y"])
            err = st_["median_obs"][pos] - st_["y"][pos]
            summary.append(dict(fold=fold, train=trn, test=tst, method=k, primary=(k == PRIMARY),
                                n=int(te.size), n_zero=int((~pos).sum()), n_blocks=int(np.unique(blocks[te]).size),
                                mean_log_score_all=float(st_["log_score"].mean()),
                                mean_log_score_positive=float(st_["log_score"][pos].mean()),
                                rmse_dex=float(np.sqrt(np.mean(err ** 2))), bias_dex=float(err.mean()),
                                spearman=float(spearmanr(st_["y"][pos], st_["median_obs"][pos]).statistic),
                                coverage95=float(np.nanmean(st_["covered"])), mean_crps_dex=float(np.nanmean(st_["crps"])),
                                obs_sd_sgris=a.obs_sd_sgris))
            for jj, i in enumerate(te):
                per_sample.append(dict(fold=fold, sample=df["sample"][i], block=blocks[i], method=k,
                                       log_score=st_["log_score"][jj], median_obs_pred=st_["median_obs"][jj],
                                       log_b_obs=st_["y"][jj], covered=st_["covered"][jj], crps=st_["crps"][jj]))
        if tst == "s6_2017":                          # published retrieval as a reference point (no density)
            pub = ED.cook2020_published_inversion()
            cv = df["sample"].iloc[te].map(pub).astype(float).to_numpy()
            y = np.log10(df.cells.iloc[te].to_numpy().clip(min=1e-9))
            ok = (cv > 0) & (df.cells.iloc[te].to_numpy() > 0)
            summary.append(dict(fold=fold, train="(published)", test=tst, method="cook2020_published", n=int(ok.sum()),
                                rmse_dex=float(np.sqrt(np.mean((np.log10(cv[ok]) - y[ok]) ** 2))),
                                bias_dex=float(np.mean(np.log10(cv[ok]) - y[ok]))))
        # contrasts of the pre-registered primary against every other method (descriptive except the
        # pre-registered primary contrast), cluster bootstrap over sampling days
        if PRIMARY in stats_by:
            for k in methods:
                if k == PRIMARY:
                    continue
                for metric in ("log_score", "crps"):
                    a_, b_ = stats_by[PRIMARY][metric], stats_by[k][metric]
                    if metric == "crps":
                        a_, b_ = -a_, -b_                    # positive = primary better
                    paired.append(dict(fold=fold, comparison=f"{PRIMARY} - {k}", metric=metric,
                                       preregistered=(k == PRIMARY_CONTRAST and metric == "log_score"),
                                       **block_bootstrap(a_, b_, blocks[te])))
    fw = pd.DataFrame(fwd_rows)
    fsum = []
    for (fold, m), g in fw.groupby(["fold", "method"]):
        e = (g.bb_hcrf_pred - g.bb_hcrf_meas).to_numpy(float)
        fsum.append(dict(fold=fold, method=m, n=int(np.isfinite(g.logp_hcrf_bands).sum()),
                         mean_logp_hcrf_bands=float(np.nanmean(g.logp_hcrf_bands)),
                         bb_hcrf_mae=float(np.nanmean(np.abs(e))) if np.isfinite(e).any() else np.nan,
                         bb_hcrf_bias=float(np.nanmean(e)) if np.isfinite(e).any() else np.nan,
                         max_irradiance_share_300_350=float(g.irradiance_share_300_350.max())))
    for fold in fw.fold.unique():
        g = fw[fw.fold == fold].pivot(index="sample", columns="method", values="logp_hcrf_bands")
        bl = fw[fw.fold == fold].drop_duplicates("sample").set_index("sample").block.reindex(g.index).to_numpy()
        for k in g.columns:
            if k != PRIMARY and PRIMARY in g:
                paired.append(dict(fold=fold, comparison=f"{PRIMARY} - {k}", metric="forward_logp_hcrf_bands",
                                   preregistered=False, **block_bootstrap(g[PRIMARY], g[k], bl)))
    import provenance as PV
    for name, obj in (("heldout_summary.csv", pd.DataFrame(summary)), ("heldout_per_sample.csv", pd.DataFrame(per_sample)),
                      ("heldout_paired.csv", pd.DataFrame(paired)), ("heldout_forward_per_sample.csv", fw),
                      ("heldout_forward_summary.csv", pd.DataFrame(fsum))):
        PV.atomic_to_csv(obj, os.path.join(a.outdir, name), index=False, float_format="%.5g")
    PV.atomic_write_text(os.path.join(a.outdir, "heldout_settings.json"), json.dumps(dict(
        settings=settings, obs_sd_sgris=a.obs_sd_sgris, v_zero_ml=a.v_zero_ml, environment=PV.environment(),
        code=PV.code_fingerprint(), runtime_s=round(time.time() - t0)), indent=1, default=float))
    pd.set_option("display.width", 220)
    print(pd.DataFrame(summary).round(3).to_string(index=False))
    print(pd.DataFrame(fsum).round(3).to_string(index=False))
    print(pd.DataFrame(paired).round(3).to_string(index=False))


if __name__ == "__main__":
    main()
