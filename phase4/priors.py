"""
Empirical priors for the per-pixel Bayesian inversion.

Every prior is derived from published field data (phase2/empirical_data.py;
citations in data/empirical/SOURCES.md):

  log10 B  ~ Normal(3.56, 0.78)  fitted to 180 natural surface-ice samples with cells > 0 at
             site S6, July-August 2016 (Williamson et al. 2020 data deposit). These samples are
             independent of the 2017 field-validation samples (Cook et al. 2020).
  f_n      ~ Beta(17.4, 11.4)    method-of-moments fit to the three published Greenland community
             surveys (A. nordenskioeldii fraction 0.65, 0.66, 0.50).
  r        ~ log-uniform on [300, 20000] um (bubbly-ice optical radius): the range spanned by the
             field NIR reflectances at the measured ice densities; Cooper et al. (2021) measured
             ~9.3-10.6 mm in this ice. No measured distribution exists, so the prior is flat in log r.
  k        ~ Normal(0.90, 0.175) the band-averaged anisotropic reflectance factor HCRF/albedo of 51
             field spectra (biosnicar-py ARF_master.csv), linking directional reflectance to albedo.
  dust     log-uniform over its nodes, only if the dust axis is enabled. Field studies at S6 found
             local mineral dust weakly absorbing (Cook et al. 2020; Tedstone et al. 2020), so the
             default model has no dust axis.

Elevation, slope and melt stage: no published calibration links them quantitatively to algal
abundance or ice structure at the pixel scale, so they are NOT used in the prior. Positive degree
days from PROMICE station data are still computed and mapped (empirical_data.pdd_at_elevation).
"""

from __future__ import annotations

import os
import sys
from dataclasses import dataclass, field

import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "phase2"))
import empirical_data as ED  # noqa: E402


def _defaults():
    mu_b, sd_b, n_b, _ = ED.abundance_prior()
    fa, fb, _, _ = ED.community_prior()
    mk, sk, _ = ED.anisotropy_prior()
    return dict(mu_b=mu_b, sd_b=sd_b, f_alpha=fa, f_beta=fb, mu_k=mk, sd_k=sk)


@dataclass
class PriorConfig:
    mu_b: float = field(default_factory=lambda: _defaults()["mu_b"])
    sd_b: float = field(default_factory=lambda: _defaults()["sd_b"])
    f_alpha: float = field(default_factory=lambda: _defaults()["f_alpha"])
    f_beta: float = field(default_factory=lambda: _defaults()["f_beta"])
    mu_k: float = field(default_factory=lambda: _defaults()["mu_k"])
    sd_k: float = field(default_factory=lambda: _defaults()["sd_k"])
    scale: float = 1.0                          # multiply prior SDs (sensitivity test only)

    def describe(self):
        return dict(log_b=f"Normal({self.mu_b:.3f}, {self.sd_b:.3f})",
                    f_n=f"Beta({self.f_alpha:.2f}, {self.f_beta:.2f})",
                    r_um=f"log-uniform {ED.ICE_RADIUS_BOUNDS_UM}", k=f"Normal({self.mu_k:.3f}, {self.sd_k:.3f})")


def prior_logpdfs(emu_axes: dict, n_pixels: int, cfg: PriorConfig):
    """Log-prior arrays (normalised over the grid nodes) for each axis; shape (1, n_axis) as the
    empirical priors are spatially uniform. Returns (log_priors, sd_k, mu_k, mu_b_per_pixel)."""
    sc = cfg.scale
    out = {}
    lb = emu_axes["log_b"]
    out["log_b"] = stats.norm.logpdf(lb, cfg.mu_b, cfg.sd_b * sc)[None, :]
    fn = emu_axes["f_n"]
    if len(fn) > 1:
        fc = np.clip(fn, 1e-3, 1 - 1e-3)
        # widen by sc: keep the mean, reduce the concentration
        conc = (cfg.f_alpha + cfg.f_beta) / sc ** 2
        m = cfg.f_alpha / (cfg.f_alpha + cfg.f_beta)
        out["f_n"] = stats.beta.logpdf(fc, m * conc, (1 - m) * conc)[None, :]
    else:
        out["f_n"] = np.zeros((1, 1))
    r = emu_axes["r_um"]
    # log-uniform density on the (possibly non-uniform) node set: weight each node by its share of log r
    lr = np.log(r)
    edges = np.concatenate([[lr[0]], 0.5 * (lr[1:] + lr[:-1]), [lr[-1]]]) if len(r) > 1 else np.array([0, 1])
    out["r_um"] = np.log(np.maximum(np.diff(edges), 1e-12))[None, :]
    if "dust_ppb" in emu_axes:
        out["dust_ppb"] = np.zeros((1, len(emu_axes["dust_ppb"])))
    for k, v in out.items():
        out[k] = v - np.logaddexp.reduce(v, axis=1, keepdims=True)
    return out, cfg.sd_k * sc, cfg.mu_k, np.full(n_pixels, cfg.mu_b)


def mcmc_prior_params(cfg: PriorConfig):
    return dict(mu_b=cfg.mu_b, sd_b=cfg.sd_b * cfg.scale, f_alpha=cfg.f_alpha, f_beta=cfg.f_beta,
                mu_k=cfg.mu_k, r_prior="loguniform")
