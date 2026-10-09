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
sys.path.insert(0, os.path.join(HERE, "..", "..", "phase2"))

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
        gp.half_logden = 0.5 * np.log(1.0 + sk ** 2 * gp.a)
        ll, _, _ = gp._loglik((R / sig)[None, :], sk, mk=mk)
        for n in range(5):
            logf = (-0.5 * np.sum(((R[None, :] - ks[:, None] * F[n]) / sig) ** 2, axis=1)
                    - 0.5 * 4 * np.log(2 * np.pi) - np.log(sig).sum() + stats.norm.logpdf(ks, mk, sk))
            num = logsumexp(logf) + np.log(dk)
            assert abs((ll[0, n] + gp.norm_const) - num) < 1e-2 * max(1.0, abs(num) * 1e-3), \
                (sk, mk, n, ll[0, n] + gp.norm_const, num)


def test_loglik_equals_multivariate_normal_and_is_float64_accurate():
    """Sigma = diag(sigma^2) + s_k^2 F F^T with mean m_k F, checked against scipy at bright-ice
    reflectance with small sigma (where float32 cancellation in rr - 2G + a was ~1e-3-1e-2)."""
    rng = np.random.default_rng(1)
    F = rng.uniform(0.55, 0.75, size=(50, 4))
    sig = np.array([0.004, 0.004, 0.005, 0.006])
    R = 0.9 * F[7] + rng.normal(0, sig)
    gp = INV.GridPosterior(_FakeEm(F), sig, derived=())
    sk, mk = 0.175, 0.9
    gp.half_logden = 0.5 * np.log(1.0 + sk ** 2 * gp.a)
    ll, _, _ = gp._loglik((R / sig)[None, :], sk, mk)
    ref = np.array([stats.multivariate_normal(mk * f, np.diag(sig ** 2) + sk ** 2 * np.outer(f, f)).logpdf(R)
                    for f in F])
    assert np.max(np.abs(ll[0] + gp.norm_const - ref)) < 1e-8


def test_derived_quantities_not_zero_filled():
    """A derived quantity undefined (NaN) on part of the grid must not be averaged as 0."""
    F = np.tile(np.array([0.6, 0.58, 0.55, 0.45]), (4, 1)) * np.array([1.0, 0.99, 0.98, 0.97])[:, None]
    em = _FakeEm(F)
    em.data["bba"] = np.array([0.5, 0.5, np.nan, np.nan])[None, None, :]
    sig = np.full(4, 0.05)
    lp = {"log_b": np.zeros((1, 1)), "f_n": np.zeros((1, 1)), "r_um": np.full((1, 4), -np.log(4))}
    res = INV.GridPosterior(em, sig, derived=("bba",)).run(F[0][None, :], lp, 0.1)
    assert np.isnan(res["bba_mean"][0]) and res["bba_nonfinite_mass"][0] > 0.3      # was 0.5 x (1 - mass)
    em.data["bba"] = np.array([0.5, 0.5, 0.5, np.nan])[None, None, :]
    lp["r_um"] = np.log(np.array([[0.5, 0.5, 1e-12, 1e-12]]))
    res = INV.GridPosterior(em, sig, derived=("bba",)).run(F[0][None, :], lp, 0.1)
    assert abs(res["bba_mean"][0] - 0.5) < 1e-9 and 0 <= res["ppp"][0] <= 1


def test_mcmc_prior_params_match_grid_priors():
    from priors import PriorConfig, mcmc_prior_params, prior_logpdfs
    for sc in (1.0, 1.5):
        pc = PriorConfig(scale=sc)
        pp = mcmc_prior_params(pc)
        f = np.linspace(0.0005, 0.9995, 2000)
        lp, _, _, _ = prior_logpdfs({"log_b": np.array([3.0]), "f_n": f, "r_um": np.array([3000.0]),
                                     "dust_ppb": np.geomspace(1e3, 1e7, 400)}, 1, pc)
        wf = np.exp(lp["f_n"][0])
        assert abs(wf @ f - pp["f_alpha"] / (pp["f_alpha"] + pp["f_beta"])) < 2e-3
        assert abs(np.sqrt(wf @ f ** 2 - (wf @ f) ** 2) - stats.beta(pp["f_alpha"], pp["f_beta"]).std()) < 2e-3
        ld = np.log(np.geomspace(1e3, 1e7, 400))
        wd = np.exp(lp["dust_ppb"][0])
        assert abs(wd @ ld - pp["mu_lndust"]) < 0.05 and abs(np.sqrt(wd @ ld ** 2 - (wd @ ld) ** 2) - pp["sd_lndust"]) < 0.05


class _LinEm:
    """Continuous fake emulator with a dust axis: reflectance falls linearly with log B and log dust."""
    def __init__(self):
        self.axes = {"log_b": np.linspace(1, 6, 11), "f_n": np.linspace(0, 1, 5), "r_um": np.geomspace(500, 20000, 6),
                     "dust_ppb": np.geomspace(3e4, 1.2e6, 6)}
        self.names = list(self.axes)
        self.active = self.names
        g = np.meshgrid(*[self.coord(n) for n in self.names], indexing="ij")
        base = 0.8 - 0.05 * g[0] - 0.01 * g[1] - 0.02 * np.log(g[2] / 1000) - 0.03 * (g[3] - 4.5)
        self.data = {"bands": np.stack([base * c for c in (1.0, 0.98, 0.95, 0.8)], axis=-1)}

    def coord(self, n):
        v = self.axes[n]
        return np.log10(v + 100.0) if n == "dust_ppb" else v

    def interpolator(self, key="bands"):
        from scipy.interpolate import RegularGridInterpolator
        return RegularGridInterpolator(tuple(self.coord(n) for n in self.active), self.data[key],
                                       bounds_error=False, fill_value=None)


def test_mcmc_uses_dust_prior_and_independent_ensembles():
    from priors import PriorConfig, mcmc_prior_params
    em = _LinEm()
    pp = mcmc_prior_params(PriorConfig())
    R = np.array([0.5, 0.49, 0.47, 0.4])
    m = INV.mcmc_pixel(em, R, np.full(4, 10.0), pp, 0.175, n_walkers=24, n_steps=3000, burn=1000, n_ensembles=3)
    j = m["names"].index("dust_ppb")
    ld = np.log(m["samples"][:, j])
    lo, hi = np.log(em.axes["dust_ppb"][[0, -1]])
    tn = stats.truncnorm((lo - pp["mu_lndust"]) / pp["sd_lndust"], (hi - pp["mu_lndust"]) / pp["sd_lndust"],
                         pp["mu_lndust"], pp["sd_lndust"])
    assert abs(ld.mean() - tn.mean()) < 0.1 and abs(ld.std() - tn.std()) < 0.1     # uninformative data -> prior
    assert m["n_ensembles"] == 3 and np.all(m["rhat"] < 1.05) and "rhat_walkers" in m


def test_split_rhat_detects_disagreeing_chains():
    rng = np.random.default_rng(0)
    same = rng.normal(size=(1000, 4, 1))
    apart = same + np.array([0.0, 0.0, 3.0, 3.0])[None, :, None]
    assert INV._split_rhat(same)[0] < 1.01 and INV._split_rhat(apart)[0] > 1.5


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


def _fake_run(tmp_path, rf=20.0, bba=0.5, sza=47.0, date="2019-07-23", sw_model=800.0, complete=True):
    import json
    import provenance as PV
    import s2_io as io
    from rasterio.transform import from_origin
    d = tmp_path / "run"
    (d / "geotiff").mkdir(parents=True)
    names = ("ours_log_b_mean", "ours_bba_mean", "ours_rf_algae_mean", "tierA_rf_algae_mean", "tierA_log_b_mean")
    vals = (3.5, bba, rf, 2 * rf, 3.4)
    io.write_geotiff(str(d / "geotiff" / "phase4_maps.tif"), {n: np.full((4, 4), v, np.float32) for n, v in zip(names, vals)},
                     from_origin(0, 0, 20, 20), "EPSG:32622")
    s = dict(scene="S2A_22WEV_20190723_0_L2A", date=date, sza=sza, overpass_datetime=f"{date}T15:14:03Z",
             sw_down_model_w_m2=sw_model)
    (d / "summary.json").write_text(json.dumps(s))
    if complete:
        (d / "COMPLETE.json").write_text(json.dumps(dict(artifacts={
            "summary.json": PV.file_sha256(str(d / "summary.json")),
            "geotiff/phase4_maps.tif": PV.file_sha256(str(d / "geotiff" / "phase4_maps.tif"))})))
    return d


def test_multi_scene_daily_forcing_uses_model_sw_and_overpass_time(tmp_path):
    import multi_scene as MS
    d = _fake_run(tmp_path)
    row = MS.summarise(str(d))
    sw_mean, sw_ov = MS.sw_at("2019-07-23", "2019-07-23T15:14:03Z")
    assert np.isclose(row["dalpha_rf_ours"], 20.0 / 800.0)
    assert np.isclose(row["rf_ours_daily_mean"], 20.0 / 800.0 * sw_mean)        # model SW, not measured
    assert np.isclose(row["rf_ours_daily_mean_legacy_formula"], 20.0 * sw_mean / sw_ov)
    # the overpass SW is interpolated at the actual time (hour-centre convention), not read at hour 15
    sw = MS._promice_day("2019-07-23")
    h = 15 + 14 / 60 + 3 / 3600
    assert np.isclose(sw_ov, np.interp(h, np.arange(24) + 0.5, sw)) and not np.isclose(sw_ov, sw[15])
    assert np.isclose(row["melt_potential_rf_ours_cm_we_d"], row["rf_ours_daily_mean"] * 86400 / 3.34e6)


def test_multi_scene_completion_requires_marker(tmp_path):
    import multi_scene as MS
    d = _fake_run(tmp_path, complete=False)
    assert MS.run_complete(str(d))[0] is False                  # summary.json alone is not completion
    d2 = _fake_run(tmp_path / "b")
    assert MS.run_complete(str(d2))[0] is True
    with open(d2 / "summary.json", "a") as fh:
        fh.write(" ")
    assert MS.run_complete(str(d2))[0] is False                 # changed after completion
