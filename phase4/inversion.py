"""
Bayesian inversion of Sentinel-2 surface reflectance.

Observation model (per pixel, bands b = B2, B3, B4, B8):

    R_b = k * F_b(z) + eps_b,   eps_b ~ N(0, sigma_b^2),   k ~ N(1, s_k^2)

F(z) is the BioSNICAR emulator, z the state (log10 B, f_n, r[, dust]) and k a
multiplicative nuisance (illumination / anisotropy / residual atmosphere).

1) Grid quadrature (all pixels). The state space is 3-4 dimensional, so the
   posterior is evaluated EXACTLY on the (refined) emulator grid - no sampling
   noise, no convergence issues. k is integrated out analytically: with whitened
   vectors r = R/sigma, f = F/sigma, a = f.f and d = r - f,

       p(R | z) = N(d; 0, I + s_k^2 f f^T)
       -2 log p = d.d - s_k^2 (f.d)^2 / (1 + s_k^2 a) + log(1 + s_k^2 a) + const

   (Sherman-Morrison / matrix determinant lemma). For a chunk of pixels this is one
   matrix product r @ f^T plus element-wise operations.
   Outputs: posterior mean / SD / 2.5-50-97.5 % quantiles of each state variable,
   MAP state, posterior mean and SD of derived quantities (pigment mass, broadband
   albedo, radiative forcing), the log evidence (for model comparison / Bayes
   factors) and a chi-square goodness-of-fit at the MAP.

2) MCMC (selected pixels). emcee affine-invariant ensemble sampler on the
   continuous (interpolated) emulator with k explicit, to cross-check the grid
   posterior and to draw corner plots; diagnostics: integrated autocorrelation
   time, effective sample size and split-R-hat.
"""

from __future__ import annotations

import logging

import numpy as np
from scipy.special import logsumexp

logging.getLogger("emcee").setLevel(logging.ERROR)     # diagnostics are returned, not logged


class GridPosterior:
    def __init__(self, emulator, sigma, derived=("pigment_ug_l", "bba", "rf_algae", "rf_total")):
        self.em = emulator
        self.sigma = np.asarray(sigma, dtype=float)
        F = emulator.data["bands"].reshape(-1, 4)
        self.Fw = (F / self.sigma).astype(np.float32)
        self.a = (self.Fw.astype(np.float64) ** 2).sum(axis=1)
        self.a32 = self.a.astype(np.float32)
        self.half_logden = None
        self.shape = emulator.shape
        self.derived = {k: emulator.data[k].reshape(-1) for k in derived if k in emulator.data}
        self.norm_const = -0.5 * 4 * np.log(2 * np.pi) - np.log(self.sigma).sum()

    def _loglik(self, Rw, sk):
        G = Rw @ self.Fw.T                                   # (p, N), float32
        rr = (Rw ** 2).sum(axis=1, keepdims=True)
        a = self.a32[None, :]
        s2 = sk * sk
        c = (s2 / (1.0 + s2 * self.a32))[None, :]
        fd = G - a
        ll = rr - 2.0 * G
        ll += a
        ll -= c * fd * fd
        ll *= -0.5
        ll -= self.half_logden[None, :]
        return ll, G, rr

    def run(self, R, log_priors: dict, sk: float, chunk: int | None = None, keep_full: np.ndarray | None = None):
        """R: (P, 4) reflectance. log_priors: axis -> (P or 1, n_axis). Returns dict of arrays."""
        self.half_logden = (0.5 * np.log(1.0 + sk ** 2 * self.a)).astype(np.float32)
        R = np.asarray(R, dtype=float)
        P = R.shape[0]
        N = self.Fw.shape[0]
        chunk = chunk or max(8, int(1.5e7 // N))
        names = self.em.names
        res = {f"{n}_{s}": np.full(P, np.nan) for n in names for s in ("mean", "sd", "q025", "q50", "q975", "map")}
        for k in self.derived:
            res[f"{k}_mean"] = np.full(P, np.nan)
            res[f"{k}_sd"] = np.full(P, np.nan)
        res.update(log_evidence=np.full(P, np.nan), chi2=np.full(P, np.nan), k_map=np.full(P, np.nan))
        full = {}
        ok = np.all(np.isfinite(R), axis=1)
        idx_all = np.flatnonzero(ok)
        for c0 in range(0, idx_all.size, chunk):
            idx = idx_all[c0:c0 + chunk]
            Rw = (R[idx] / self.sigma).astype(np.float32)
            ll, G, rr = self._loglik(Rw, np.float32(sk))
            lp = ll.reshape(len(idx), *self.shape)
            for ax, n in enumerate(names):
                pri = log_priors[n]
                pri = pri[idx] if pri.shape[0] == P and P > 1 else np.broadcast_to(pri, (len(idx), pri.shape[1]))
                shp = [len(idx)] + [1] * len(names)
                shp[ax + 1] = pri.shape[1]
                lp = lp + pri.reshape(shp)
            flat = lp.reshape(len(idx), -1)
            logz = logsumexp(flat, axis=1).astype(np.float32)
            w = np.exp(flat - logz[:, None]).astype(np.float64)
            res["log_evidence"][idx] = logz + self.norm_const
            imap = flat.argmax(axis=1)
            for ax, n in enumerate(names):
                vals = self.em.axes[n]
                other = tuple(i + 1 for i in range(len(names)) if i != ax)
                pm = w.reshape(len(idx), *self.shape).sum(axis=other)        # (p, n_axis)
                m = pm @ vals
                res[f"{n}_mean"][idx] = m
                res[f"{n}_sd"][idx] = np.sqrt(np.clip(pm @ vals ** 2 - m ** 2, 0, None))
                cdf = np.cumsum(pm, axis=1)
                for q, key in ((0.025, "q025"), (0.5, "q50"), (0.975, "q975")):
                    res[f"{n}_{key}"][idx] = _quantile_from_cdf(vals, cdf, q)
                res[f"{n}_map"][idx] = vals[np.unravel_index(imap, self.shape)[ax]]
            for k, g in self.derived.items():
                ok_g = np.isfinite(g)
                if not ok_g.any():
                    continue
                gg = np.where(ok_g, g, 0.0)
                m = w @ gg
                res[f"{k}_mean"][idx] = m
                res[f"{k}_sd"][idx] = np.sqrt(np.clip(w @ gg ** 2 - m ** 2, 0, None))
            Gm = G[np.arange(len(idx)), imap]
            am = self.a[imap]
            kh = (1.0 / sk ** 2 + Gm) / (1.0 / sk ** 2 + am)
            res["k_map"][idx] = kh
            res["chi2"][idx] = rr[:, 0] - 2 * kh * Gm + kh ** 2 * am
            if keep_full is not None:
                for j, p in enumerate(idx):
                    if p in keep_full:
                        full[int(p)] = w[j].reshape(self.shape)
        res["full_posteriors"] = full
        return res


def _quantile_from_cdf(vals, cdf, q):
    out = np.empty(cdf.shape[0])
    for i in range(cdf.shape[0]):
        out[i] = np.interp(q, cdf[i], vals) if cdf[i, -1] > 0 else np.nan
    return out


# --------------------------------------------------------------------------- #
# MCMC cross-check                                                             #
# --------------------------------------------------------------------------- #
def mcmc_pixel(emulator, R, sigma, prior_params: dict, sk: float, n_walkers: int = 32,
               n_steps: int = 6000, burn: int = 2000, seed: int = 0):
    """emcee sampling of (state..., k) for one pixel on the continuous emulator.

    prior_params: dict with mu_b, sd_b, f_alpha, f_beta, mu_r, sd_r (as used by the grid).
    Returns dict(samples, names, tau, ess, rhat, acceptance).
    """
    import emcee
    from scipy import stats

    I = emulator.interpolator("bands")
    act = emulator.active
    lo = np.array([emulator.coord(n)[0] for n in act])
    hi = np.array([emulator.coord(n)[-1] for n in act])
    R = np.asarray(R, dtype=float)
    sig = np.asarray(sigma, dtype=float)
    pp = prior_params
    ib, ir = act.index("log_b"), act.index("r_um")
    i_f = act.index("f_n") if "f_n" in act else None

    def logp(th):
        z, k = th[:-1], th[-1]
        if np.any(z < lo) or np.any(z > hi):
            return -np.inf
        lpr = (stats.norm.logpdf(z[ib], pp["mu_b"], pp["sd_b"]) + stats.norm.logpdf(k, 1.0, sk)
               + stats.norm.logpdf(z[ir], pp["mu_r"], pp["sd_r"]))
        if i_f is not None:
            lpr += stats.beta.logpdf(np.clip(z[i_f], 1e-6, 1 - 1e-6), pp["f_alpha"], pp["f_beta"])
        F = I(z[None, :])[0]
        return lpr - 0.5 * np.sum(((R - k * F) / sig) ** 2)

    rng = np.random.default_rng(seed)
    ndim = len(act) + 1
    p0 = np.empty((n_walkers, ndim))
    p0[:, :-1] = lo + (hi - lo) * rng.uniform(0.25, 0.75, size=(n_walkers, len(act)))
    p0[:, ib] = np.clip(pp["mu_b"] + 0.1 * rng.normal(size=n_walkers), lo[ib] + 0.01, hi[ib] - 0.01)
    p0[:, -1] = 1.0 + 0.01 * rng.normal(size=n_walkers)
    sampler = emcee.EnsembleSampler(n_walkers, ndim, logp, moves=[(emcee.moves.DEMove(), 0.8),
                                                                    (emcee.moves.DESnookerMove(), 0.2)])
    sampler.run_mcmc(p0, n_steps, progress=False)
    chain = sampler.get_chain(discard=burn)                    # (steps, walkers, ndim)
    tau = sampler.get_autocorr_time(discard=burn, quiet=True)
    flat = chain.reshape(-1, ndim)
    return dict(samples=flat, names=act + ["k"], tau=tau, ess=flat.shape[0] / np.maximum(tau, 1.0),
                rhat=_split_rhat(chain), acceptance=float(np.mean(sampler.acceptance_fraction)))


def _split_rhat(chain):
    """Split-R-hat (Gelman et al. 2013) treating each walker half as a chain."""
    n = chain.shape[0] // 2
    c = np.concatenate([chain[:n], chain[n:2 * n]], axis=1)       # (n, 2W, ndim)
    m = c.shape[1]
    means = c.mean(axis=0)
    W = c.var(axis=0, ddof=1).mean(axis=0)
    B = n * means.var(axis=0, ddof=1)
    var = (n - 1) / n * W + B / n
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.sqrt(var / W)
