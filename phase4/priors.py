"""
Environmental priors for the per-pixel Bayesian inversion.

The priors encode glaciological knowledge as explicit, tunable hyperparameters
(they are assumptions - report them, and test sensitivity with --prior-scale):

  melt stage s in [0,1] from positive degree days (PDD), with air temperature
  lapsed from a reference station:  T(z, d) = T_ref(d) - Gamma (z - z_ref),
  PDD(z) = sum_d max(T(z, d), 0),   s = PDD / (PDD + PDD_half)

  log10 B   ~ Normal(mu_B, sd_B) truncated to the emulator range,
              mu_B = mu0 + a_melt (s - 0.5) - a_slope * slope / 10 deg
              (longer melt seasons favour blooms; steep, well-drained surfaces
               flush cells away)
  f_n       ~ Beta(alpha, beta)   (community ratio, from microscopy counts if available)
  r         ~ Normal(r0 + r_melt * s, sd_r) truncated (weathering crust coarsens with melt)
  dust      ~ log-uniform over the emulator's dust nodes (if present)
  k         ~ Normal(1, sd_k)     multiplicative reflectance factor (illumination/slope,
              HCRF-vs-albedo anisotropy, residual atmospheric error); marginalised
              analytically in the likelihood
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy import stats


@dataclass
class PriorConfig:
    mu0: float = 3.8            # log10 cells/mL at mid melt stage (Dark Zone summer ~10^3-10^4.5)
    sd_b: float = 0.8
    a_melt: float = 0.8
    a_slope: float = 0.3
    f_alpha: float = 2.0
    f_beta: float = 2.0
    r0: float = 1200.0
    r_melt: float = 1400.0
    sd_r: float = 600.0
    sd_k: float = 0.10
    # melt model
    z_ref: float = 1000.0       # m a.s.l. of the reference station (e.g. PROMICE S6 ~1010 m)
    lapse: float = 6.5e-3       # K m^-1
    pdd_half: float = 150.0     # degree-days at which s = 0.5
    scale: float = 1.0          # multiply all prior SDs (sensitivity test)


def climatological_tref(doy_end: int, t_peak: float = 3.0, doy_peak: int = 205, width: float = 40.0,
                        doy_start: int = 121):
    """Default daily mean air temperature at z_ref (deg C) if no station record is given:
    a smooth summer bump peaking at t_peak on doy_peak. Replace with AWS data (PROMICE) when
    available via --tref-csv (columns: date, t_air)."""
    d = np.arange(doy_start, doy_end + 1)
    return t_peak - 12.0 * ((d - doy_peak) / (2.5 * width)) ** 2


def melt_stage(elev_m: np.ndarray, t_ref_series: np.ndarray, cfg: PriorConfig) -> tuple[np.ndarray, np.ndarray]:
    """PDD (deg C d) and melt stage s for each pixel elevation."""
    z = np.asarray(elev_m, dtype=float)
    dz = (z - cfg.z_ref)[..., None]
    T = t_ref_series[None, ...] - cfg.lapse * dz if z.ndim == 1 else t_ref_series - cfg.lapse * dz
    pdd = np.clip(T, 0, None).sum(axis=-1)
    return pdd, pdd / (pdd + cfg.pdd_half)


def prior_logpdfs(emu_axes: dict, melt: np.ndarray, slope_deg: np.ndarray, cfg: PriorConfig):
    """Separable log-prior arrays for each pixel and axis.

    Returns dict axis -> array (P, n_axis) (or (1, n_axis) if pixel-independent),
    plus the k prior SD. Each row is normalised over the axis nodes.
    """
    s = np.asarray(melt, dtype=float).ravel()
    sl = np.asarray(slope_deg, dtype=float).ravel()
    sc = cfg.scale
    out = {}
    lb = emu_axes["log_b"]
    mu_b = cfg.mu0 + cfg.a_melt * (s - 0.5) - cfg.a_slope * sl / 10.0
    out["log_b"] = stats.norm.logpdf(lb[None, :], mu_b[:, None], cfg.sd_b * sc)
    fn = emu_axes["f_n"]
    if len(fn) > 1:
        fc = np.clip(fn, 1e-3, 1 - 1e-3)
        a = 1 + (cfg.f_alpha - 1) / sc
        b = 1 + (cfg.f_beta - 1) / sc
        out["f_n"] = stats.beta.logpdf(fc, a, b)[None, :]
    else:
        out["f_n"] = np.zeros((1, 1))
    r = emu_axes["r_um"]
    mu_r = cfg.r0 + cfg.r_melt * s
    out["r_um"] = stats.norm.logpdf(r[None, :], mu_r[:, None], cfg.sd_r * sc)
    if "dust_ppb" in emu_axes:
        out["dust_ppb"] = np.zeros((1, len(emu_axes["dust_ppb"])))
    for k, v in out.items():                                     # normalise over nodes
        out[k] = v - np.logaddexp.reduce(v, axis=1, keepdims=True)
    return out, cfg.sd_k * sc, mu_b
