"""
Standard empirical baseline: a band-ratio algal index regressed on biomass.

Index:  I = (B4 - B2) / (B4 + B2)      (red-blue normalised difference; algal phenolic
        pigments absorb blue far more than red, so I rises with biomass)
Model:  log10 B = c0 + c1 I + c2 I^2

Published Sentinel-2 glacier-algae indices use the red-edge band B5 (e.g. the 2-band
705/665 nm ratio); with the four 10 m bands used here we adopt the red-blue index.
Without field-calibrated coefficients, the regression is calibrated on simulations with
BioSNICAR's default empirical algal optics (tier A) over the full prior range of grain size
and illumination factor k - i.e. the "standard empirical" retrieval a user of the default
model would build. Pass --empirical-coefs c0 c1 c2 to use field-calibrated values instead.
"""

from __future__ import annotations

import numpy as np


def index(R):
    R = np.asarray(R, dtype=float)
    return (R[..., 2] - R[..., 0]) / (R[..., 2] + R[..., 0])


def calibrate(emulator_A, n: int = 20000, k_sd: float = 0.1, noise: float = 0.01, seed: int = 3,
              log_b_range=(2.0, 5.5)):
    """Least-squares fit of log10 B on the index using random tier-A emulator states."""
    rng = np.random.default_rng(seed)
    I = emulator_A.interpolator("bands")
    act = emulator_A.active
    pts = np.empty((n, len(act)))
    for j, name in enumerate(act):
        c = emulator_A.coord(name)
        pts[:, j] = rng.uniform(c[0], c[-1], n)
    pts[:, act.index("log_b")] = rng.uniform(*log_b_range, n)
    F = I(pts)
    R = F * rng.normal(1.0, k_sd, (n, 1)) + rng.normal(0, noise, F.shape)
    x = index(R)
    X = np.column_stack([np.ones(n), x, x ** 2])
    coef, *_ = np.linalg.lstsq(X, pts[:, act.index("log_b")], rcond=None)
    resid = pts[:, act.index("log_b")] - X @ coef
    return coef, float(np.sqrt(np.mean(resid ** 2)))


def predict(R, coef):
    x = index(R)
    return coef[0] + coef[1] * x + coef[2] * x ** 2
