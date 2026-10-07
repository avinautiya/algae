"""
Synthetic ground truth for testing the retrieval (no field data needed).

Truth fields live on a real georeferenced grid (the chosen AOI; real DEM if
available) and follow the same glaciological structure as the priors, plus
spatially correlated random variability:

    log10 B  = empirical prior mean + GRF(empirical SD, ~300 m)       clipped [1.5, 5.5]
    f_n      = logistic(logit(empirical community mean) + GRF(sd 0.8))
    r        = log-spaced field inside the empirical radius bounds
    k        = empirical anisotropy mean + smooth GRF(0.25 x empirical SD)

Observations are computed with DIRECT BioSNICAR runs at the exact (off-grid) states
(r snapped only to BioSNICAR's own 20 um LUT), multiplied by k, plus Gaussian noise.
This avoids the 'inverse crime' of generating and inverting with the same emulator.
`truth_model` selects the optics used to generate truth ('ours' or 'tierA') so the
retrieval can be stress-tested in both directions.
"""

from __future__ import annotations

import numpy as np
from scipy.ndimage import gaussian_filter


def grf(shape, sd, corr_px, rng):
    z = gaussian_filter(rng.normal(size=shape), corr_px, mode="reflect")
    return sd * z / z.std()


def make_truth(shape, prior_cfg, res_m, seed=7, with_dust=False):
    """Spatially correlated truth drawn around the EMPIRICAL priors (priors.PriorConfig):
    log10 B ~ N(mu_b, sd_b) field, f_n around the community mean, log r uniform-ish field, k ~ N(mu_k, 0.25 sd_k)."""
    import empirical_data as ED
    rng = np.random.default_rng(seed)
    corr = 300.0 / res_m
    m_f = prior_cfg.f_alpha / (prior_cfg.f_alpha + prior_cfg.f_beta)
    lo, hi = np.log(ED.ICE_RADIUS_BOUNDS_UM[0] * 1.5), np.log(ED.ICE_RADIUS_BOUNDS_UM[1] / 1.5)
    z = grf(shape, 1.0, corr, rng)
    t = dict(
        log_b=np.clip(prior_cfg.mu_b + grf(shape, prior_cfg.sd_b, corr, rng), 1.5, 5.5),
        f_n=1.0 / (1.0 + np.exp(-(np.log(m_f / (1 - m_f)) + grf(shape, 0.8, corr, rng)))),
        r_um=np.exp(lo + (hi - lo) * 0.5 * (1 + np.tanh(z))),
        k=prior_cfg.mu_k + grf(shape, 0.25 * prior_cfg.sd_k, 2 * corr, rng),
    )
    if with_dust:
        t["dust_ppb"] = np.clip(10 ** (3.5 + grf(shape, 0.4, corr, rng)), 0, 1e5)
    return t


def synthesise(truth, emu_cfg, phase1_l2=None, demo=False, biosnicar=None, noise=(0.01, 0.01, 0.01, 0.012),
               seed=11, workers=1):
    """Observed reflectance (H, W, 4) and true broadband albedo / forcing (H, W)."""
    import emulator as E

    H, W = truth["log_b"].shape
    keys = [k for k in ("log_b", "f_n", "r_um", "dust_ppb") if k in truth]
    states = [{k: float(truth[k].flat[i]) for k in keys} for i in range(H * W)]
    E._BUILDER = E._Builder(emu_cfg, phase1_l2, demo, biosnicar)      # warm in parent before fork
    E._BUILDER.runner.illumination(emu_cfg.sza)
    for r in np.unique([E._BUILDER.runner.snap_radius(s["r_um"], emu_cfg.ice_mode) for s in states]):
        E._BUILDER.runner.ice(E._BUILDER.spec(r))
    import multiprocessing as mp
    chunks = [states[i:i + 256] for i in range(0, len(states), 256)]
    if workers > 1 and "fork" in mp.get_all_start_methods():
        with mp.get_context("fork").Pool(workers) as pool:
            out = np.concatenate(pool.map(E.direct_forward_chunk, chunks))
    else:
        out = np.concatenate([E.direct_forward_chunk(c) for c in chunks])
    rng = np.random.default_rng(seed)
    F = out[:, :4].reshape(H, W, 4)
    R = truth["k"][..., None] * F + rng.normal(size=F.shape) * np.asarray(noise)
    extra = dict(bba=out[:, 4].reshape(H, W), rf_algae=out[:, 5].reshape(H, W),
                 rf_total=out[:, 6].reshape(H, W), pigment_ug_l=out[:, 7].reshape(H, W))
    return R.astype(np.float32), extra
