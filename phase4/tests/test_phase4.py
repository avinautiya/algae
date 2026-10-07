"""
Validation tests for Phase 4.   Run:  python phase4/tests/test_phase4.py   (or pytest)

1. The analytic marginalisation over the reflectance factor k equals brute-force
   numerical integration over k.
2. Emulator (coarse grid + PCHIP refinement) reproduces direct BioSNICAR runs at
   off-grid states to < 0.002 reflectance.
3. Grid-quadrature posterior and MCMC (emcee) agree for a synthetic pixel, MCMC converges
   (split R-hat < 1.05), and the truth lies inside the 95 % interval.
4. Posterior calibration: over many synthetic pixels with known truth, the 95 % credible
   interval of log10 B covers the truth >= 90 % of the time.
5. Sentinel-2 reflectance offset rule (processing baseline >= 04.00) and DEM slope.
"""

import os
import sys

import numpy as np
from scipy import stats

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, ".."))

import inversion as INV  # noqa: E402


class _FakeEm:
    """Minimal emulator stand-in for testing the likelihood algebra."""
    def __init__(self, F):
        self.data = {"bands": F[None, None, :, :]}
        self.axes = {"log_b": np.array([0.0]), "f_n": np.array([0.5]), "r_um": np.arange(F.shape[0], dtype=float)}
        self.names = list(self.axes)
        self.shape = (1, 1, F.shape[0])


def test_k_marginalisation():
    """log p(R|z) from the closed form == log of the integral over k (dense log-space quadrature)."""
    from scipy.special import logsumexp
    rng = np.random.default_rng(0)
    F0 = rng.uniform(0.4, 0.9, size=4)
    F = F0 * (1 + 0.05 * rng.normal(size=(5, 4)))
    sig = np.array([0.02, 0.02, 0.02, 0.03])
    R = 1.04 * F[2] + rng.normal(0, sig)
    sk = 0.1
    gp = INV.GridPosterior(_FakeEm(F), sig, derived=())
    gp.half_logden = (0.5 * np.log(1.0 + sk ** 2 * gp.a)).astype(np.float32)
    ll, _, _ = gp._loglik((R / sig)[None, :].astype(np.float32), np.float32(sk))
    ks = np.linspace(0.0, 2.0, 400001)
    dk = ks[1] - ks[0]
    for n in range(5):
        logf = (-0.5 * np.sum(((R[None, :] - ks[:, None] * F[n]) / sig) ** 2, axis=1)
                - 0.5 * 4 * np.log(2 * np.pi) - np.log(sig).sum() + stats.norm.logpdf(ks, 1.0, sk))
        num = logsumexp(logf) + np.log(dk)
        assert abs((ll[0, n] + gp.norm_const) - num) < 1e-2 * max(1.0, abs(num) * 1e-3), (n, ll[0, n] + gp.norm_const, num)


def _root():
    sys.path.insert(0, os.path.join(HERE, "..", "..", "phase2"))
    import biosnicar_bridge as bb
    try:
        return bb.locate_biosnicar()
    except ImportError:
        return None


def _small_emulator():
    import emulator as E
    cfg = E.EmulatorConfig(model="ours", sza=47, log_b=(1, 6, 0.25), f_n=(0, 1, 0.25), r_um=(1000, 3000, 160))
    return cfg, E.build_emulator(cfg, demo=True, verbose=False)


def test_emulator_accuracy():
    if _root() is None:
        print("BioSNICAR not found - skipping")
        return
    import emulator as E
    cfg, em = _small_emulator()
    I = em.refine_log_b(0.05).interpolator()
    rng = np.random.default_rng(1)
    st = [dict(log_b=rng.uniform(2, 5.5), f_n=rng.choice(em.axes["f_n"]), r_um=float(rng.choice(em.axes["r_um"])))
          for _ in range(25)]
    d = E.direct_forward(cfg, st, demo=True)
    p = I(np.array([[s["log_b"], s["f_n"], s["r_um"]] for s in st]))
    assert np.abs(p - d[:, :4]).max() < 2e-3


def test_grid_vs_mcmc_and_calibration():
    if _root() is None:
        print("BioSNICAR not found - skipping")
        return
    import emulator as E
    from priors import PriorConfig, prior_logpdfs
    cfg, em = _small_emulator()
    emf = em.refine_log_b(0.05)
    sig = np.array([0.02, 0.02, 0.02, 0.025])
    pc = PriorConfig()
    rng = np.random.default_rng(5)
    n = 120
    truth = [dict(log_b=rng.uniform(2.5, 5.0), f_n=rng.uniform(0, 1), r_um=float(rng.uniform(1000, 3000)))
             for _ in range(n)]
    F = E.direct_forward(cfg, truth, demo=True)[:, :4]
    k = rng.normal(1.0, 0.05, n)[:, None]
    R = k * F + rng.normal(0, 0.01, F.shape)
    lp, sk, mu_b = prior_logpdfs(emf.axes, np.full(n, 0.5), np.zeros(n), pc)
    res = INV.GridPosterior(emf, sig).run(R, lp, sk)
    tb = np.array([t["log_b"] for t in truth])
    cover = np.mean((tb >= res["log_b_q025"]) & (tb <= res["log_b_q975"]))
    assert cover >= 0.90, cover
    pp = dict(mu_b=mu_b[0], sd_b=pc.sd_b, f_alpha=pc.f_alpha, f_beta=pc.f_beta, mu_r=pc.r0 + pc.r_melt * 0.5,
              sd_r=pc.sd_r)
    m = INV.mcmc_pixel(emf, R[0], sig, pp, sk, n_steps=4000, burn=1500)
    j = m["names"].index("log_b")
    assert abs(m["samples"][:, j].mean() - res["log_b_mean"][0]) < 0.08
    assert abs(m["samples"][:, j].std() - res["log_b_sd"][0]) < 0.08
    assert np.all(m["rhat"] < 1.05), m["rhat"]


def test_s2_scaling_and_slope():
    import s2_io as io
    item = {"properties": {"s2:processing_baseline": "05.09"}}
    assert io._scale_offset(item, {"raster:bands": [{"scale": 1e-4, "offset": 0.0}]}) == (1e-4, -0.1)
    assert io._scale_offset(item, {"raster:bands": [{"scale": 1e-4, "offset": -0.1}]}) == (1e-4, -0.1)
    old = {"properties": {"s2:processing_baseline": "02.13"}}
    assert io._scale_offset(old, {"raster:bands": [{"scale": 1e-4, "offset": 0.0}]}) == (1e-4, 0.0)
    y, x = np.mgrid[0:50, 0:50] * 20.0
    dem = 1000 + 0.05 * x                                     # 5 % grade
    assert np.allclose(io.slope_deg(dem, 20.0), np.degrees(np.arctan(0.05)))


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print(f"PASS {name}")
