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
    gp = INV.GridPosterior(_FakeEm(F), sig, derived=())
    ks = np.linspace(-1.0, 3.0, 800001)
    dk = ks[1] - ks[0]
    for sk, mk in ((0.1, 1.0), (0.175, 0.898)):          # unit-mean and the empirical HCRF/albedo prior
        gp.half_logden = (0.5 * np.log(1.0 + sk ** 2 * gp.a)).astype(np.float32)
        ll, _, _ = gp._loglik((R / sig)[None, :].astype(np.float32), np.float32(sk), mk=np.float32(mk))
        for n in range(5):
            logf = (-0.5 * np.sum(((R[None, :] - ks[:, None] * F[n]) / sig) ** 2, axis=1)
                    - 0.5 * 4 * np.log(2 * np.pi) - np.log(sig).sum() + stats.norm.logpdf(ks, mk, sk))
            num = logsumexp(logf) + np.log(dk)
            assert abs((ll[0, n] + gp.norm_const) - num) < 1e-2 * max(1.0, abs(num) * 1e-3), \
                (sk, mk, n, ll[0, n] + gp.norm_const, num)


def _root():
    sys.path.insert(0, os.path.join(HERE, "..", "..", "phase2"))
    import biosnicar_bridge as bb
    try:
        return bb.locate_biosnicar()
    except ImportError:
        return None


def _small_emulator():
    import emulator as E
    # default (empirical) configuration: bubbly ice at the measured densities, ESA SRFs, Williamson
    # et al. (2020) in vivo pigment MACs; a reduced r grid keeps the test quick
    cfg = E.EmulatorConfig(model="ours", phenol="williamson2020", sza=47, log_b=(1, 6, 0.25), f_n=(0, 1, 0.25),
                           r_um=(2000, 16000, 6))
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
    truth = [dict(log_b=rng.normal(pc.mu_b, pc.sd_b), f_n=rng.beta(pc.f_alpha, pc.f_beta),
                  r_um=float(rng.choice(emf.axes["r_um"]))) for _ in range(n)]
    for t in truth:
        t["log_b"] = float(np.clip(t["log_b"], 1.2, 5.8))
    F = E.direct_forward(cfg, truth, demo=True)[:, :4]
    k = rng.normal(pc.mu_k, 0.05, n)[:, None]
    R = k * F + rng.normal(0, 0.01, F.shape)
    lp, sk, mk, mu_b = prior_logpdfs(emf.axes, n, pc)
    res = INV.GridPosterior(emf, sig).run(R, lp, sk, mk=mk)
    tb = np.array([t["log_b"] for t in truth])
    cover = np.mean((tb >= res["log_b_q025"]) & (tb <= res["log_b_q975"]))
    assert cover >= 0.90, cover
    from priors import mcmc_prior_params
    pp = mcmc_prior_params(pc)
    m = INV.mcmc_pixel(emf, R[0], sig, pp, sk, n_steps=4000, burn=1500)
    j = m["names"].index("log_b")
    assert abs(m["samples"][:, j].mean() - res["log_b_mean"][0]) < 0.08
    assert abs(m["samples"][:, j].std() - res["log_b_sd"][0]) < 0.08
    assert np.all(m["rhat"] < 1.05), m["rhat"]


def test_s2_scaling_and_slope():
    import pytest
    import s2_io as io
    rb = lambda off: {"raster:bands": [{"scale": 1e-4, **({} if off is None else {"offset": off})}]}  # noqa: E731
    new_flag = {"properties": {"s2:processing_baseline": "05.09", "earthsearch:boa_offset_applied": True}}
    new_old_meta = {"properties": {"s2:processing_baseline": "04.00", "earthsearch:boa_offset_applied": False}}
    old = {"properties": {"s2:processing_baseline": "02.13"}}
    assert io._scale_offset(new_flag, rb(-0.1))[:2] == (1e-4, -0.1)              # recorded, applied ONCE
    with pytest.raises(ValueError):
        io._scale_offset(new_flag, rb(0.0))                                       # contradictory metadata
    assert io._scale_offset(new_old_meta, rb(0.0))[:2] == (1e-4, -0.1)           # pre-flag items
    assert io._scale_offset(new_old_meta, rb(None))[:2] == (1e-4, -0.1)
    assert io._scale_offset(new_old_meta, rb(-0.1))[:2] == (1e-4, -0.1)          # never -0.2
    assert io._scale_offset(old, rb(0.0))[:2] == (1e-4, 0.0)
    with pytest.raises(ValueError):
        io._scale_offset({"properties": {}}, rb(None))                            # unknown baseline
    y, x = np.mgrid[0:50, 0:50] * 20.0
    dem = 1000 + 0.05 * x                                     # 5 % grade
    assert np.allclose(io.slope_deg(dem, 20.0), np.degrees(np.arctan(0.05)))


def _fake_tile(tmp_path, origin=(499980.0, 7500000.0), n10=60, epsg=32622):
    """10 m reflectance bands with a known pattern and a 20 m SCL, on a tile grid whose origin is not a
    multiple of 20 m (like real tiles)."""
    import rasterio
    from rasterio.transform import from_origin
    yy, xx = np.mgrid[0:n10, 0:n10]
    assets = {}
    for i, key in enumerate(("blue", "green", "red", "nir")):
        dn = (5000 + 10 * xx + 1000 * i).astype(np.uint16)
        dn[0, 0] = 0                                                   # one nodata pixel
        p = tmp_path / f"{key}.tif"
        with rasterio.open(p, "w", driver="GTiff", height=n10, width=n10, count=1, dtype="uint16", nodata=0,
                           crs=f"EPSG:{epsg}", transform=from_origin(origin[0], origin[1], 10, 10)) as d:
            d.write(dn, 1)
        assets[key] = {"href": str(p), "raster:bands": [{"scale": 1e-4, "offset": 0.0}],
                       "proj:transform": [10, 0, origin[0], 0, -10, origin[1]]}
    scl = np.full((n10 // 2, n10 // 2), 11, np.uint8)
    scl[:, :5] = 8                                                     # cloud in the 5 westernmost 20 m columns
    p = tmp_path / "scl.tif"
    with rasterio.open(p, "w", driver="GTiff", height=n10 // 2, width=n10 // 2, count=1, dtype="uint8", nodata=0,
                       crs=f"EPSG:{epsg}", transform=from_origin(origin[0], origin[1], 20, 20)) as d:
        d.write(scl, 1)
    assets["scl"] = {"href": str(p)}
    props = {"s2:processing_baseline": "02.13", "proj:epsg": epsg, "view:sun_elevation": 42.8861289918669,
             "datetime": "2019-07-23T15:14:03.547000Z"}
    return {"id": "fake", "properties": props, "assets": assets}


def test_s2_reader_grid_alignment_and_mask(tmp_path):
    import pytest
    import s2_io as io
    item = _fake_tile(tmp_path)
    x0, y0 = 499980.0, 7500000.0
    # requested 20 m window that is NOT on the 20 m grid; aligned_bounds snaps it outward
    b = io.aligned_bounds((x0 + 15, y0 - 405, x0 + 395, y0 - 25), (x0, y0), 20.0)
    assert b == (x0, y0 - 420, x0 + 400, y0 - 20)
    with pytest.raises(ValueError):
        io._dst_grid((0, 0, 30, 40), 20.0)                     # not whole pixels
    sc = io.read_scene(item, b, resolution=20.0, aoi_lonlat=(-49.86, 67.13), max_sun_mismatch_deg=1.0)
    assert sc.transform.c == b[0] and sc.transform.f == b[3] and sc.shape == (20, 20)
    # B2 at 20 m = mean of the 2x2 10 m block: DN 5000 + 10*(2j + 0.5)
    j = np.arange(20)
    expect = (5000 + 10 * (2 * j + 0.5)) * 1e-4
    good = sc.mask[5]
    assert np.allclose(sc.bands["B2"][5][good], expect[good], atol=1e-6)
    assert not sc.mask[:, :5].any() and sc.mask[:, 5:].all()   # SCL cloud columns line up exactly
    assert sc.item["_read"]["scaling"]["B2"]["offset"] == 0.0
    assert abs(sc.sza - sc.item["_read"]["sun"]["sza_tile_mean"]) < 1.0


def test_s2_reader_reprojects_and_checks_sun(tmp_path):
    import pytest
    import s2_io as io
    from rasterio.crs import CRS
    item = _fake_tile(tmp_path)
    # bounds in a different UTM zone (23N) are reprojected, not assumed to be in the tile CRS
    from pyproj import Transformer
    t = Transformer.from_crs(32622, 32623, always_xy=True)
    cx, cy = t.transform(499980 + 300, 7500000 - 300)
    cx, cy = round(cx, -1), round(cy, -1)
    sc = io.read_scene(item, (cx - 100, cy - 100, cx + 100, cy + 100), resolution=20.0,
                       dst_crs=CRS.from_epsg(32623), max_sun_mismatch_deg=5.0)
    assert sc.crs == CRS.from_epsg(32623) and sc.shape == (10, 10) and np.isfinite(sc.bands["B3"]).any()
    bad = dict(item, properties=dict(item["properties"], datetime="2019-07-23T03:14:03Z"))   # wrong time
    with pytest.raises(ValueError):
        io.read_scene(bad, io.aligned_bounds((499980, 7499600, 499980 + 400, 7500000), (499980, 7500000), 20.0))
    nosun = dict(item, properties={k: v for k, v in item["properties"].items() if k != "view:sun_elevation"})
    with pytest.raises(ValueError):
        io.read_scene(nosun, (499980, 7499600, 499980 + 400, 7500000))


def test_solar_zenith_matches_sentinel2_metadata():
    """Tile 22WEV centre (UTM 22N 554880 E, 7445100 N), S2A 2019-07-23 15:14:03 UTC: Element 84 records sun
    elevation 42.886 deg (tile mean)."""
    import s2_io as io
    from pyproj import Transformer
    lon, lat = Transformer.from_crs(32622, 4326, always_xy=True).transform(554880, 7445100)
    assert abs(io.solar_zenith_deg(lat, lon, "2019-07-23T15:14:03.547Z") - (90 - 42.886)) < 0.3


if __name__ == "__main__":
    for name, fn in list(globals().items()):
        if name.startswith("test_"):
            fn()
            print(f"PASS {name}")
