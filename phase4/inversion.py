"""
Bayesian inversion of Sentinel-2 surface reflectance.

Observation model (per pixel, bands b = B2, B3, B4, B8):

    R_b = k * F_b(z) + eps_b,   eps_b ~ N(0, sigma_b^2),   k ~ N(m_k, s_k^2)

F(z) is the BioSNICAR emulator, z the state (log10 B, f_n, r[, dust]) and k a
multiplicative nuisance (illumination / anisotropy / residual atmosphere).

1) Grid quadrature (all pixels). The state space is 3-4 dimensional, so the
   posterior is evaluated EXACTLY on the (refined) emulator grid - no sampling
   noise, no convergence issues. k is integrated out analytically: with whitened
   vectors r = R/sigma, f = F/sigma, a = f.f and d = r - f,

       p(R | z) = N(d; 0, I + s_k^2 f f^T),   d = r - m_k f
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
        self.Fw = (F / self.sigma).astype(np.float64)       # float64: rr - 2G + a cancels ~1e4 -> O(1)
        self.a = (self.Fw ** 2).sum(axis=1)
        self.half_logden = None
        self.shape = emulator.shape
        self.derived = {k: emulator.data[k].reshape(-1) for k in derived if k in emulator.data}
        self.norm_const = -0.5 * 4 * np.log(2 * np.pi) - np.log(self.sigma).sum()

    def _loglik(self, Rw, sk, mk=1.0):
        """k ~ N(mk, sk^2) integrated out: d = r - mk f ~ N(0, I + sk^2 f f^T)."""
        G = Rw @ self.Fw.T                                   # (p, N), float64
        rr = (Rw ** 2).sum(axis=1, keepdims=True)
        a = self.a[None, :]
        s2 = sk * sk
        c = (s2 / (1.0 + s2 * self.a))[None, :]
        fd = G - mk * a
        ll = rr - (2.0 * mk) * G
        ll += (mk * mk) * a
        ll -= c * fd * fd
        ll *= -0.5
        ll -= self.half_logden[None, :]
        return ll, G, rr

    def run(self, R, log_priors: dict, sk: float, chunk: int | None = None, keep_full: np.ndarray | None = None,
            mk: float = 1.0):
        """R: (P, 4) reflectance. log_priors: axis -> (P or 1, n_axis). Returns dict of arrays."""
        self.half_logden = 0.5 * np.log(1.0 + sk ** 2 * self.a)
        R = np.asarray(R, dtype=float)
        P = R.shape[0]
        N = self.Fw.shape[0]
        chunk = chunk or max(8, int(1e7 // N))
        names = self.em.names
        res = {f"{n}_{s}": np.full(P, np.nan) for n in names for s in ("mean", "sd", "q025", "q50", "q975", "map")}
        for k in self.derived:
            res[f"{k}_mean"] = np.full(P, np.nan)
            res[f"{k}_sd"] = np.full(P, np.nan)
            res[f"{k}_nonfinite_mass"] = np.full(P, np.nan)
        res.update(log_evidence=np.full(P, np.nan), chi2=np.full(P, np.nan), k_map=np.full(P, np.nan),
                   mahal_map=np.full(P, np.nan), ppp=np.full(P, np.nan))
        full = {}
        ok = np.all(np.isfinite(R), axis=1)
        idx_all = np.flatnonzero(ok)
        for c0 in range(0, idx_all.size, chunk):
            idx = idx_all[c0:c0 + chunk]
            Rw = R[idx] / self.sigma
            ll, G, rr = self._loglik(Rw, float(sk), float(mk))
            # k-marginal Mahalanobis distance of every node: D = -2 (ll + half_logden) ~ chi2(4) given z
            D = -2.0 * (ll + self.half_logden[None, :])
            lp = ll.reshape(len(idx), *self.shape)
            for ax, n in enumerate(names):
                pri = log_priors[n]
                pri = pri[idx] if pri.shape[0] == P and P > 1 else np.broadcast_to(pri, (len(idx), pri.shape[1]))
                shp = [len(idx)] + [1] * len(names)
                shp[ax + 1] = pri.shape[1]
                lp = lp + pri.reshape(shp)
            flat = lp.reshape(len(idx), -1)
            logz = logsumexp(flat, axis=1)
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
                # nodes where a derived quantity is undefined are NOT counted as zero: the posterior is
                # renormalised over finite nodes, the missing mass reported, and > 1e-3 gives NaN
                ok_g = np.isfinite(g)
                if not ok_g.any():
                    continue
                gg = np.where(ok_g, g, 0.0)
                wf = w @ ok_g.astype(float)
                m = (w @ gg) / np.maximum(wf, 1e-300)
                sd = np.sqrt(np.clip((w @ gg ** 2) / np.maximum(wf, 1e-300) - m ** 2, 0, None))
                bad = (1.0 - wf) > 1e-3
                res[f"{k}_mean"][idx] = np.where(bad, np.nan, m)
                res[f"{k}_sd"][idx] = np.where(bad, np.nan, sd)
                res[f"{k}_nonfinite_mass"][idx] = 1.0 - wf
            Gm = G[np.arange(len(idx)), imap]
            am = self.a[imap]
            kh = (mk / sk ** 2 + Gm) / (1.0 / sk ** 2 + am)
            res["k_map"][idx] = kh
            res["chi2"][idx] = rr[:, 0] - 2 * kh * Gm + kh ** 2 * am        # conditional on k at its mode
            res["mahal_map"][idx] = D[np.arange(len(idx)), imap]
            # posterior predictive p-value of the k-marginal discrepancy: sum_z w(z) P(chi2_4 >= D(z)).
            # With 4 bands and 3-4 states plus k this check has little power (it is conservative).
            from scipy.stats import chi2 as _chi2
            res["ppp"][idx] = np.sum(w * _chi2.sf(D, 4), axis=1)
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
               n_steps: int = 6000, burn: int = 2000, seed: int = 0, init=None, n_ensembles: int = 4):
    """emcee sampling of (state..., k) for one pixel on the continuous emulator, with the SAME priors as
    the grid (priors.mcmc_prior_params): Normal log B, Beta f_n, log-normal r (sampled as ln r),
    log-normal dust (sampled as ln dust, if the dust axis is active) and Normal k.

    Convergence: `n_ensembles` independent ensembles, each started from its own prior draws (seeded
    seed, seed+1, ...). Walkers of one ensemble are coupled by the moves, so R-hat across walkers is
    not a test of independent chains; `rhat` is the split-R-hat across ensembles (each ensemble's
    pooled samples as one chain, split in halves). `rhat_walkers` (across walkers of ensemble 0) is
    kept for reference only. `init` (grid-posterior draws) is used only when n_ensembles == 1.

    Returns dict(samples, names, tau, ess, rhat, rhat_walkers, acceptance, n_ensembles)."""
    import emcee
    from scipy import stats

    I = emulator.interpolator("bands")
    act = emulator.active
    ir = act.index("r_um")
    idust = act.index("dust_ppb") if "dust_ppb" in act else None
    logax = [ir] + ([idust] if idust is not None else [])          # sampled in ln
    lo = np.array([emulator.coord(n)[0] for n in act], dtype=float)
    hi = np.array([emulator.coord(n)[-1] for n in act], dtype=float)
    for j in logax:                                     # bounds of the SAMPLED variable: ln of the axis values
        ax = emulator.axes[act[j]]
        lo[j], hi[j] = np.log(max(ax[0], 1e-12)), np.log(ax[-1])

    def to_coord(z):
        """Sampled variables -> the emulator's interpolation coordinates (r in um; dust as log10(d + 100))."""
        zc = z.copy()
        zc[ir] = np.exp(z[ir])
        if idust is not None:
            zc[idust] = np.log10(np.exp(z[idust]) + 100.0)
        return zc
    R = np.asarray(R, dtype=float)
    sig = np.asarray(sigma, dtype=float)
    pp = prior_params
    if idust is not None and "mu_lndust" not in pp:
        raise ValueError("dust axis active but no dust prior (mu_lndust, sd_lndust) given")
    ib = act.index("log_b")
    i_f = act.index("f_n") if "f_n" in act else None

    def log_prior(z, k):
        lpr = stats.norm.logpdf(z[ib], pp["mu_b"], pp["sd_b"]) + stats.norm.logpdf(k, pp.get("mu_k", 1.0), sk)
        lr = z[ir]                                      # densities below are in ln r
        if pp.get("r_prior", "loguniform") == "lognormal":
            lpr += stats.norm.logpdf(lr, pp["mu_lnr"], pp["sd_lnr"])
        elif pp["r_prior"] == "normal":
            lpr += stats.norm.logpdf(np.exp(lr), pp["mu_r"], pp["sd_r"]) + lr
        if i_f is not None:
            lpr += stats.beta.logpdf(np.clip(z[i_f], 1e-3, 1 - 1e-3), pp["f_alpha"], pp["f_beta"])
        if idust is not None:
            lpr += stats.norm.logpdf(z[idust], pp["mu_lndust"], pp["sd_lndust"])
        return lpr

    def logp(th):
        z, k = th[:-1], th[-1]
        if np.any(z < lo) or np.any(z > hi):
            return -np.inf
        F = I(to_coord(z)[None, :])[0]
        return log_prior(z, k) - 0.5 * np.sum(((R - k * F) / sig) ** 2)

    ndim = len(act) + 1

    def prior_draws(rng):
        """Over-dispersed start: draws from the prior, truncated to the emulator support."""
        p0 = np.empty((n_walkers, ndim))
        for j, n in enumerate(act):
            if j == ib:
                v = rng.normal(pp["mu_b"], pp["sd_b"], n_walkers)
            elif j == i_f:
                v = rng.beta(pp["f_alpha"], pp["f_beta"], n_walkers)
            elif j == ir and pp.get("r_prior") == "lognormal":
                v = rng.normal(pp["mu_lnr"], pp["sd_lnr"], n_walkers)
            elif j == idust:
                v = rng.normal(pp["mu_lndust"], pp["sd_lndust"], n_walkers)
            else:
                v = rng.uniform(lo[j], hi[j], n_walkers)
            p0[:, j] = np.clip(v, lo[j] + 1e-6 * (hi[j] - lo[j]), hi[j] - 1e-6 * (hi[j] - lo[j]))
        p0[:, -1] = rng.normal(pp.get("mu_k", 1.0), sk, n_walkers)
        return p0

    chains, taus, acc = [], [], []
    for e in range(n_ensembles):
        rng = np.random.default_rng(seed + e)
        if init is not None and n_ensembles == 1:
            idx = rng.choice(len(init), size=n_walkers, replace=len(init) < n_walkers)
            p0 = np.asarray(init, float)[idx].copy()
            p0[:, logax] = np.log(p0[:, logax])
            p0[:, :-1] = np.clip(p0[:, :-1] + 1e-3 * (hi - lo) * rng.normal(size=(n_walkers, len(act))),
                                 lo + 1e-6, hi - 1e-6)
            p0[:, -1] += 0.25 * sk * rng.normal(size=n_walkers)
        else:
            p0 = prior_draws(rng)
        sampler = emcee.EnsembleSampler(n_walkers, ndim, logp, moves=[(emcee.moves.DEMove(), 0.8),
                                                                        (emcee.moves.DESnookerMove(), 0.2)])
        sampler.run_mcmc(p0, n_steps, progress=False, skip_initial_state_check=True)
        chains.append(sampler.get_chain(discard=burn))          # (steps, walkers, ndim)
        taus.append(sampler.get_autocorr_time(discard=burn, quiet=True))
        acc.append(float(np.mean(sampler.acceptance_fraction)))
    rhat_w = _split_rhat(chains[0])
    rhat = _split_rhat(np.stack([c.reshape(-1, ndim) for c in chains], axis=1)) if n_ensembles > 1 else rhat_w
    tau = np.max(taus, axis=0)
    out = []
    for c in chains:
        c = c.copy()
        c[..., logax] = np.exp(c[..., logax])                   # back to um / ppb
        out.append(c.reshape(-1, ndim))
    flat = np.concatenate(out)
    return dict(samples=flat, names=act + ["k"], tau=tau, ess=flat.shape[0] / np.maximum(tau, 1.0),
                rhat=rhat, rhat_walkers=rhat_w, acceptance=float(np.mean(acc)), n_ensembles=n_ensembles)


def _split_rhat(chain):
    """Split-R-hat (Gelman et al. 2013): chain (n, m, ndim) with m chains; each is split in halves."""
    n = chain.shape[0] // 2
    c = np.concatenate([chain[:n], chain[n:2 * n]], axis=1)       # (n, 2W, ndim)
    m = c.shape[1]
    means = c.mean(axis=0)
    W = c.var(axis=0, ddof=1).mean(axis=0)
    B = n * means.var(axis=0, ddof=1)
    var = (n - 1) / n * W + B / n
    with np.errstate(invalid="ignore", divide="ignore"):
        return np.sqrt(var / W)
