"""
Re-implementation (NumPy, no TensorFlow) of the *published* Dark Zone retrieval, for use as the strongest applicable
baseline ("Baseline E": existing method, not ours).

Source (read-only, pinned): Zenodo 10.5281/zenodo.18826013, data_software.zip, md5 0955d9c3ddf93ff9b69229d4ac28aa17
  src/mcmc_pipeline.py  -> priors, likelihood, observation slice [0:851:5]
  data/emulator/mlp_icelaps_101025.h5 -> 12 dense layers, ELU, last layer linear (weights read with h5py)
  data/emulator/minmax_scaler_icelaps_101025.save -> NOT unpickled; constants below were read from its raw bytes
  data/spectra/HCRF_ASDFieldspec4_spectra_Chevrollier_et_al_2026.nc -> 105 HCRF spectra, 350-2500 nm at 1 nm

What is and is not reproduced
  * the emulator forward model is reproduced exactly as a function (no reference outputs are shipped, so equality to the
    TensorFlow outputs is NOT tested here);
  * the TARGET DENSITY (Student-t, nu=3, sigma=0.003, the stated Beta priors, the stated Jacobians) is re-implemented;
  * we compute the MAXIMUM A POSTERIORI point (multi-start L-BFGS-B), not the NUTS posterior; MAP in the unconstrained
    coordinates includes the sigmoid Jacobians exactly as the TFP target does. Multi-start spread printed below is
    OPTIMISER variability, not an uncertainty estimate.
  * the paper's "1.8 % average error" definition is not given; we report both MAE and the mean relative error.
Assumed (paper ordering, not checked against pickle metadata): scaler columns = [rho1, rho2, dz1, dz2, dust, cryoconite, algae].
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

import h5py
import numpy as np
from scipy.optimize import minimize
from scipy.special import expit, logit
from scipy.stats import beta as beta_dist

ARCHIVE = os.environ.get("DARKZONE_ARCHIVE", "")
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "..", "records", "inversion_m1", "darkzone_repro")

# MinMaxScaler constants read from raw bytes of the .save file (data_min_, data_max_)
DATA_MIN = np.array([250.0, 250.0, 1e-5, 1e-5, 0.0, 0.0, 0.0])
DATA_MAX = np.array([924.6, 924.7, 0.15, 2.5, 9.98e6, 4.9e6, 985000.0])
RANGE = DATA_MAX - DATA_MIN

SIGMA, NU = 0.003, 3.0
POLY = [0.02812912, -1.39078195, 5.12318554, -5.49135067, 2.72952866]
PRIOR_SCALE_DZ = (0.4666, 0.7999992)           # unit-space caps of the Beta(2, 3.5) depth priors (code l.64-76)
LAP_CAPS_PHYS = (5e6, 1.5e6, 5e5)               # physical caps used for Beta(0.95, 3) (code l.90)
LAP_CAPS_UNIT = tuple((c - DATA_MIN[4 + i]) / RANGE[4 + i] for i, c in enumerate(LAP_CAPS_PHYS))


def load_mlp(path):
    w = h5py.File(path, "r")["model_weights"]
    layers = []
    for i in range(12):
        n = "dense" if i == 0 else f"dense_{i}"
        g = w[n]["albedo_nn"][n]
        layers.append((np.asarray(g["kernel"]), np.asarray(g["bias"])))
    return layers


def mlp(layers, x):
    h = np.atleast_2d(x)
    for k, (W, b) in enumerate(layers):
        h = h @ W + b
        if k < len(layers) - 1:
            h = np.where(h > 0, h, np.expm1(np.minimum(h, 0)))      # ELU
    return h


def rho2_min(r1):
    return 0.07410701 if r1 < 0.479 else float(np.polyval(POLY[::-1], r1))


def z_to_theta(z):
    r1 = expit(z[0])
    lo = rho2_min(r1)
    r2 = lo + (1 - lo) * expit(z[1])
    return np.array([r1, r2, expit(z[2]), expit(z[3]), expit(z[4]), expit(z[5]), expit(z[6])]), lo


def log_beta_scaled(u, a, b, cap):
    """log density in z of u=sigmoid(z) when u/cap ~ Beta(a, b), including d(u)/dz and the 1/cap scale."""
    x = u / cap
    if not (0.0 < x < 1.0):
        return -np.inf
    return beta_dist.logpdf(x, a, b) - np.log(cap) + np.log(u * (1 - u))


def log_target(z, layers, y_obs):
    theta, lo = z_to_theta(z)
    pred = mlp(layers, theta)[0, 10:181]
    r = (pred - y_obs) / SIGMA
    ll = -0.5 * (NU + 1) * np.sum(np.log1p(r * r / NU))
    lp = log_beta_scaled(theta[2], 2.0, 3.5, PRIOR_SCALE_DZ[0]) + log_beta_scaled(theta[3], 2.0, 3.5, PRIOR_SCALE_DZ[1])
    for i in range(3):
        lp += log_beta_scaled(theta[4 + i], 0.95, 3.0, LAP_CAPS_UNIT[i])
    r1 = theta[0]
    ldj = np.log(r1 * (1 - r1) * (1 - lo) * expit(z[1]) * (1 - expit(z[1])) + 1e-12)
    return ll + lp + ldj, pred


def init_z(rng):
    r1 = rng.uniform(0.02, 0.98)
    lo = rho2_min(r1)
    r2 = rng.uniform(lo, 1.0)
    z0 = logit(r1)
    z1 = logit(np.clip((r2 - lo) / (1 - lo), 1e-3, 1 - 1e-3))
    z2 = logit(rng.uniform(0.02, PRIOR_SCALE_DZ[0] * 0.9))
    z3 = logit(rng.uniform(0.02, PRIOR_SCALE_DZ[1] * 0.9))
    laps = [logit(rng.uniform(1e-6, 0.25 * c)) for c in LAP_CAPS_UNIT]       # first-quartile-like start
    return np.array([z0, z1, z2, z3, *laps])


def fit_one(layers, y_obs, n_starts, rng):
    best, runs = None, []
    for _ in range(n_starts):
        res = minimize(lambda z: -(log_target(z, layers, y_obs)[0]) if np.isfinite(log_target(z, layers, y_obs)[0]) else 1e12,
                       init_z(rng), method="L-BFGS-B", options={"maxiter": 300})
        lt, pred = log_target(res.x, layers, y_obs)
        runs.append(float(lt))
        if best is None or lt > best[0]:
            best = (lt, res.x, pred)
    return best, runs


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--archive", default=ARCHIVE, help="directory holding the extracted Zenodo data/ folder")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--starts", type=int, default=6)
    ap.add_argument("--seed", type=int, default=20261031)
    a = ap.parse_args()
    if not a.archive:
        sys.exit("set --archive or DARKZONE_ARCHIVE to the extracted Zenodo folder (contains data/)")
    layers = load_mlp(os.path.join(a.archive, "data", "emulator", "mlp_icelaps_101025.h5"))
    f = h5py.File(os.path.join(a.archive, "data", "spectra", "HCRF_ASDFieldspec4_spectra_Chevrollier_et_al_2026.nc"), "r")
    wl = f["wavelength"][:]
    refl = f["reflectance"][:]                       # (wavelength, spectrum)
    tags = [t.decode() if isinstance(t, bytes) else t for t in f["measurement_tag"][:]]
    assert wl[0] == 350 and wl[850] == 1200 and refl.shape[0] == wl.size
    obs = refl[0:851:5, :]                           # same slice as the archived pipeline: 350..1200 nm at 5 nm
    assert obs.shape[0] == 171
    os.makedirs(OUT, exist_ok=True)
    rng = np.random.default_rng(a.seed)
    rows = []
    n = obs.shape[1] if not a.limit else a.limit
    t0 = time.time()
    for j in range(n):
        y = obs[:, j]
        if not np.all(np.isfinite(y)) or y.min() < -0.02 or y.max() > 1.2:
            rows.append(dict(tag=tags[j], status="skipped: non-finite or out-of-range reflectance"))
            continue
        (lt, z, pred), runs = fit_one(layers, y, a.starts, rng)
        th, _ = z_to_theta(z)
        phys = DATA_MIN + th * RANGE
        cap = np.array(LAP_CAPS_PHYS)
        res = pred - y
        rows.append(dict(
            tag=tags[j], status="ok", log_target=float(lt), mae=float(np.mean(np.abs(res))),
            mean_rel_err=float(np.mean(np.abs(res) / np.maximum(y, 1e-3))),
            rho1=phys[0], rho2=phys[1], dz1_m=phys[2], dz2_m=phys[3],
            dust_ppb=phys[4], cryoconite_ppb=phys[5], algae_cells_ml=phys[6],
            lap_at_prior_cap=[bool(phys[4 + i] > 0.98 * cap[i]) for i in range(3)],
            start_logtarget_spread=float(np.max(runs) - np.min(runs)),
        ))
        print(f"{j + 1}/{n} {tags[j]} mae={rows[-1]['mae']:.4f} rel={rows[-1]['mean_rel_err']:.3f} ({time.time() - t0:.0f}s)", flush=True)
    ok = [r for r in rows if r["status"] == "ok"]
    summ = dict(
        archive_md5="0955d9c3ddf93ff9b69229d4ac28aa17", n_spectra=len(rows), n_fitted=len(ok),
        mean_mae=float(np.mean([r["mae"] for r in ok])), mean_rel_err=float(np.mean([r["mean_rel_err"] for r in ok])),
        max_rel_err=float(np.max([r["mean_rel_err"] for r in ok])),
        published_reference="1.8 % mean, 11.8 % max (n=104); definition of 'error' not stated",
        n_cap_hits=[int(sum(r["lap_at_prior_cap"][i] for r in ok)) for i in range(3)],
        method="MAP (multi-start L-BFGS-B) of the archived target density; NOT NUTS", starts=a.starts, seed=a.seed,
    )
    json.dump(rows, open(os.path.join(OUT, "per_spectrum.json"), "w"), indent=1, default=float)
    json.dump(summ, open(os.path.join(OUT, "summary.json"), "w"), indent=1)
    print(json.dumps(summ, indent=1))


if __name__ == "__main__":
    main()
