#!/usr/bin/env python
"""
Phase 4 driver: physics-informed Bayesian inversion of Sentinel-2 surface reflectance
(B2, B3, B4, B8) for glacier-algal abundance, community fraction and ice grain size,
with posterior-propagated albedo and radiative-forcing maps.

Examples
--------
# Synthetic-truth experiment on the real Dark Zone grid (quantifies error reduction):
python run_phase4.py --source synthetic --phase1-l2 ../phase1/results/level2

# Same, but truth generated with BioSNICAR's empirical algae (reverse stress test):
python run_phase4.py --source synthetic --truth-model tierA --phase1-l2 ...

# Real Sentinel-2 L2A scene over the SW Greenland Dark Zone (S6), 6 km x 6 km at 20 m:
python run_phase4.py --source s2 --scene-id S2A_22WEV_20190723_0_L2A --phase1-l2 ... \
       --ice-mode bubbly --dust

# Search instead of a fixed scene id (needs the earth-search STAC API):
python run_phase4.py --source s2 --search 2019-07-01 2019-08-31 --lat 67.08 --lon -49.35 ...

# Your own L2A GeoTIFFs (B02.tif, B03.tif, B04.tif, B08.tif, SCL.tif in one folder):
python run_phase4.py --source local --local-dir path/to/tifs --sza 47 ...

# Field validation (CSV with lon, lat, cells_per_ml):
python run_phase4.py --source s2 ... --field-csv cell_counts.csv
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

CHI2_99_DF4 = 13.277      # 99th percentile of chi-square with 4 dof


def parse_args(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    s = p.add_argument_group("scene")
    s.add_argument("--source", choices=["synthetic", "s2", "local"], default="synthetic")
    s.add_argument("--scene-id", default=None, help="e.g. S2A_22WEV_20190723_0_L2A (default for --source s2)")
    s.add_argument("--search", nargs=2, metavar=("START", "END"), default=None)
    s.add_argument("--max-cloud", type=float, default=10.0)
    s.add_argument("--lat", type=float, default=67.08)
    s.add_argument("--lon", type=float, default=-49.35)
    s.add_argument("--size-km", type=float, default=None, help="AOI size (default 6 km real, 3 km synthetic)")
    s.add_argument("--resolution", type=float, default=20.0, help="m (10 = native; 20 = 2x2 average)")
    s.add_argument("--local-dir", default=None)
    s.add_argument("--sza", type=float, default=None, help="solar zenith (deg) for --source local")
    s.add_argument("--no-dem", action="store_true", help="skip the Copernicus DEM download")

    m = p.add_argument_group("physics")
    m.add_argument("--phase1-l2", default=os.path.join(HERE, "..", "phase1", "results", "level2"))
    m.add_argument("--demo", action="store_true", help="BioSNICAR ppg.csv instead of Phase 1 (DEMO)")
    m.add_argument("--biosnicar", default=None)
    m.add_argument("--tier", choices=["C", "D"], default="C")
    m.add_argument("--ice-mode", choices=["grains", "bubbly"], default="grains")
    m.add_argument("--rho", type=float, default=650.0)
    m.add_argument("--r-range", type=float, nargs=3, default=[1000.0, 3000.0, 80.0], metavar=("MIN", "MAX", "STEP"))
    m.add_argument("--dust", action="store_true", help="add a mineral-dust nuisance dimension")
    m.add_argument("--log-b-step", type=float, default=0.1, help="inference grid step in log10 B (dex)")
    m.add_argument("--f-step", type=float, default=0.1, help="grid step of the community fraction")
    m.add_argument("--sigma", type=float, nargs=4, default=[0.02, 0.02, 0.02, 0.025],
                   help="per-band 1-sigma reflectance uncertainty (sensor + atmosphere + model)")

    q = p.add_argument_group("priors")
    q.add_argument("--mu0", type=float, default=3.8)
    q.add_argument("--sd-b", type=float, default=0.8)
    q.add_argument("--f-alpha", type=float, default=2.0)
    q.add_argument("--f-beta", type=float, default=2.0)
    q.add_argument("--sd-k", type=float, default=0.10)
    q.add_argument("--prior-scale", type=float, default=1.0, help="multiply all prior SDs")
    q.add_argument("--t-peak", type=float, default=3.0, help="peak daily mean air T at z_ref (deg C)")
    q.add_argument("--z-ref", type=float, default=1000.0)
    q.add_argument("--tref-csv", default=None, help="station record: columns date, t_air (deg C)")
    q.add_argument("--spatial-pooling", type=float, default=0.0,
                   help="empirical-Bayes second pass: Gaussian smoothing (pixels) of the log B prior mean")

    v = p.add_argument_group("validation / output")
    v.add_argument("--truth-model", choices=["ours", "tierA"], default="ours")
    v.add_argument("--noise", type=float, nargs=4, default=[0.01, 0.01, 0.01, 0.012])
    v.add_argument("--field-csv", default=None)
    v.add_argument("--empirical-coefs", type=float, nargs=3, default=None)
    v.add_argument("--mcmc-steps", type=int, default=6000, help="per walker; first third discarded as burn-in")
    v.add_argument("--workers", type=int, default=os.cpu_count() or 1)
    v.add_argument("--outdir", default=os.path.join(HERE, "results"))
    v.add_argument("--usetex", action="store_true")
    return p.parse_args(argv)


# --------------------------------------------------------------------------- #
def load_grid(a):
    """Scene reflectance (or None for synthetic), grid geometry, SZA, date."""
    import s2_io as io
    from pyproj import CRS
    if a.source == "local":
        import glob
        files = {b: glob.glob(os.path.join(a.local_dir, f"*{b}*.tif")) for b in ("B02", "B03", "B04", "B08", "SCL")}
        item = {"properties": {"view:sun_elevation": 90 - (a.sza or 45.0), "s2:processing_baseline": "0",
                               "datetime": "local"},
                "assets": {k: {"href": files[b][0], "raster:bands": [{"scale": 1e-4, "offset": 0.0}]}
                           for k, b in zip(("blue", "green", "red", "nir", "scl"),
                                           ("B02", "B03", "B04", "B08", "SCL"))}}
        import rasterio
        with rasterio.open(files["B02"][0]) as src:
            crs = src.crs
    else:
        if a.search:
            items = io.search_scenes(a.lat, a.lon, a.search[0], a.search[1], a.max_cloud)
            if not items:
                raise SystemExit("no scene found")
            item = items[0]
        else:
            item = io.item_from_s3(a.scene_id or io.DEFAULT_SCENE)
        crs = CRS.from_epsg(item["properties"]["proj:epsg"])
    size = a.size_km or (3.0 if a.source == "synthetic" else 6.0)
    bounds = io.aoi_bounds(a.lat, a.lon, size, crs)
    scene = io.read_scene(item, bounds, a.resolution)
    if a.source == "local" and a.sza:
        scene.sza = a.sza
    return scene, bounds


def environment(a, scene, bounds):
    import s2_io as io
    from priors import PriorConfig, melt_stage, climatological_tref
    if not a.no_dem:
        try:
            dem, slope = io.read_dem(bounds, scene.crs, scene.shape, scene.transform)
        except Exception as e:
            print(f"DEM unavailable ({e}); using flat 1050 m")
            dem = np.full(scene.shape, 1050.0)
            slope = np.zeros(scene.shape)
    else:
        dem, slope = np.full(scene.shape, 1050.0), np.zeros(scene.shape)
    dem = np.where(np.isfinite(dem), dem, np.nanmedian(dem))
    slope = np.where(np.isfinite(slope), slope, 0.0)
    pc = PriorConfig(mu0=a.mu0, sd_b=a.sd_b, f_alpha=a.f_alpha, f_beta=a.f_beta, sd_k=a.sd_k, z_ref=a.z_ref,
                     scale=a.prior_scale)
    if a.tref_csv:
        t = pd.read_csv(a.tref_csv, parse_dates=["date"])
        end = pd.Timestamp(scene.date) if scene.date != "local" else t.date.max()
        tser = t[(t.date >= f"{end.year}-05-01") & (t.date <= end)].t_air.to_numpy()
    else:
        doy = pd.Timestamp(scene.date).dayofyear if scene.date not in ("", "local") else 205
        tser = climatological_tref(doy, t_peak=a.t_peak)
    pdd, s = melt_stage(dem.ravel(), tser, pc)
    return dem, slope, pdd.reshape(dem.shape), s.reshape(dem.shape), pc


def main(argv=None):
    a = parse_args(argv)
    t0 = time.time()
    tab, figd, geo, cache = (os.path.join(a.outdir, d) for d in ("tables", "figures", "geotiff", "cache"))
    for d in (tab, figd, geo, cache):
        os.makedirs(d, exist_ok=True)

    import emulator as E
    import empirical as EMP
    import inversion as INV
    import maps4 as M
    import validation as VAL
    import s2_io as io
    from priors import prior_logpdfs

    # ------------------------------------------------------------- scene & priors
    scene, bounds = load_grid(a)
    dem, slope, pdd, melt, pc = environment(a, scene, bounds)
    H, W = scene.shape
    print(f"Grid {H}x{W} @ {a.resolution:g} m, SZA {scene.sza:.1f} deg, date {scene.date}, "
          f"elev {np.nanmin(dem):.0f}-{np.nanmax(dem):.0f} m, melt stage {np.nanmin(melt):.2f}-{np.nanmax(melt):.2f}")

    # ------------------------------------------------------------- emulators
    common = dict(ice_mode=a.ice_mode, rho=a.rho, sza=round(scene.sza), r_um=tuple(a.r_range),
                  dust_ppb=(0.0, 1e3, 3e3, 1e4, 3e4, 1e5) if a.dust else ())
    cfg_o = E.EmulatorConfig(model="ours", tier=a.tier, f_n=(0.0, 1.0, a.f_step), **common)
    cfg_a = E.EmulatorConfig(model="tierA", **common)
    kw = dict(phase1_l2=None if a.demo else a.phase1_l2, demo=a.demo, biosnicar=a.biosnicar, workers=a.workers)
    em_o = E.build_emulator(cfg_o, cache=os.path.join(cache, "emulator_ours.npz"), **kw).refine_log_b(a.log_b_step)
    em_a = E.build_emulator(cfg_a, cache=os.path.join(cache, "emulator_tierA.npz"), **kw).refine_log_b(a.log_b_step)

    # ------------------------------------------------------------- observations
    truth = None
    if a.source == "synthetic":
        import synthetic as SYN
        truth = SYN.make_truth(scene.shape, melt, slope, pc, a.resolution, with_dust=a.dust)
        tcfg = cfg_o if a.truth_model == "ours" else cfg_a
        print(f"Synthesising {H * W} pixels with direct BioSNICAR ({a.truth_model} optics) ...", flush=True)
        R, extra = SYN.synthesise(truth, tcfg, workers=a.workers, noise=tuple(a.noise), **{
            k: v for k, v in kw.items() if k != "workers"})
        truth.update(extra)
        mask = np.ones(scene.shape, bool)
    else:
        R = scene.stack()
        mask = scene.mask
    Rp = R.reshape(-1, 4)
    valid = mask.ravel() & np.all(np.isfinite(Rp), axis=1)
    print(f"{valid.sum()} valid pixels")

    def to_map(v):
        out = np.full(H * W, np.nan)
        out[valid] = v
        return out.reshape(H, W)

    # ------------------------------------------------------------- empirical baseline
    if a.empirical_coefs:
        coef, cal_rmse = np.array(a.empirical_coefs), np.nan
    else:
        coef, cal_rmse = EMP.calibrate(em_a)
    logb_emp = EMP.predict(R, coef)
    logb_emp[~mask] = np.nan
    print(f"Empirical index: log10 B = {coef[0]:.2f} + {coef[1]:.2f} I + {coef[2]:.2f} I^2 "
          f"(calibration RMSE {cal_rmse:.2f} dex)")

    # ------------------------------------------------------------- Bayesian inversions
    sigma = np.array(a.sigma)
    res = {}
    for name, em in (("ours", em_o), ("tierA", em_a)):
        lp, sk, mu_b = prior_logpdfs(em.axes, melt.ravel()[valid], slope.ravel()[valid], pc)
        t1 = time.time()
        r = INV.GridPosterior(em, sigma).run(Rp[valid], lp, sk)
        if a.spatial_pooling > 0 and name == "ours":
            from scipy.ndimage import gaussian_filter
            m1 = to_map(r["log_b_mean"])
            w = gaussian_filter(np.where(np.isfinite(m1), 1.0, 0.0), a.spatial_pooling)
            sm = gaussian_filter(np.nan_to_num(m1), a.spatial_pooling) / np.maximum(w, 1e-6)
            from scipy import stats as st
            lb = em.axes["log_b"]
            mu2 = 0.5 * sm.ravel()[valid] + 0.5 * mu_b
            lp["log_b"] = st.norm.logpdf(lb[None, :], mu2[:, None], 0.5 * pc.sd_b)
            lp["log_b"] -= np.logaddexp.reduce(lp["log_b"], axis=1, keepdims=True)
            r = INV.GridPosterior(em, sigma).run(Rp[valid], lp, sk)
        res[name] = r
        bad = np.mean(r["chi2"] > CHI2_99_DF4)
        print(f"Inversion [{name}]: {time.time() - t1:.0f} s; pixels failing the chi2 test (p<0.01): {100 * bad:.1f} %")
        if bad > 0.2:
            print(f"  WARNING: {100 * bad:.0f} % of pixels are not reproduced by the [{name}] forward model. "
                  "Try --ice-mode bubbly, --dust, a wider --r-range, or larger --sigma; check the Phase 1 spectrum.")
    ro, ra = res["ours"], res["tierA"]
    log_bf = to_map(ro["log_evidence"] - ra["log_evidence"])

    # identifiability of the community fraction: posterior / prior SD
    from scipy import stats as st
    prior_sd_f = st.beta(pc.f_alpha, pc.f_beta).std()
    shrink_f = np.nanmedian(ro["f_n_sd"]) / prior_sd_f
    print(f"Community fraction identifiability: median posterior/prior SD of f_n = {shrink_f:.2f} "
          f"(1 = data uninformative)")

    maps = {
        "ours_log_b_mean": to_map(ro["log_b_mean"]), "ours_log_b_sd": to_map(ro["log_b_sd"]),
        "ours_log_b_q025": to_map(ro["log_b_q025"]), "ours_log_b_q975": to_map(ro["log_b_q975"]),
        "ours_f_n_mean": to_map(ro["f_n_mean"]), "ours_f_n_sd": to_map(ro["f_n_sd"]),
        "ours_r_um_mean": to_map(ro["r_um_mean"]), "ours_pigment_ug_l_mean": to_map(ro["pigment_ug_l_mean"]),
        "ours_bba_mean": to_map(ro["bba_mean"]), "ours_bba_sd": to_map(ro["bba_sd"]),
        "ours_rf_algae_mean": to_map(ro["rf_algae_mean"]), "ours_rf_algae_sd": to_map(ro["rf_algae_sd"]),
        "ours_chi2": to_map(ro["chi2"]), "ours_k_map": to_map(ro["k_map"]),
        "tierA_log_b_mean": to_map(ra["log_b_mean"]), "tierA_rf_algae_mean": to_map(ra["rf_algae_mean"]),
        "tierA_bba_mean": to_map(ra["bba_mean"]), "tierA_chi2": to_map(ra["chi2"]),
        "empirical_log_b": logb_emp, "log_bayes_factor_ours_vs_A": log_bf,
        "elevation_m": dem, "slope_deg": slope, "melt_stage": melt,
    }
    maps["d_rf_ours_minus_A"] = maps["ours_rf_algae_mean"] - maps["tierA_rf_algae_mean"]
    io.write_geotiff(os.path.join(geo, "phase4_maps.tif"), maps, scene.transform, scene.crs)

    # ------------------------------------------------------------- validation
    metrics = None
    if truth is not None:
        tq = {"log_b": truth["log_b"], "rf_algae": truth["rf_algae"], "bba": truth["bba"],
              "f_n": truth["f_n"], "r_um": truth["r_um"]}
        results = {
            "Empirical index": {"log_b": (logb_emp,)},
            "Tier A Bayesian": {q: (to_map(ra[f"{q}_mean"]), to_map(ra.get(f"{q}_q025", ra[f"{q}_mean"] * np.nan)),
                                    to_map(ra.get(f"{q}_q975", ra[f"{q}_mean"] * np.nan)))
                                for q in ("log_b", "rf_algae", "bba", "r_um")},
            "Physics-informed Bayesian": {q: (to_map(ro[f"{q}_mean"]),
                                              to_map(ro.get(f"{q}_q025", ro[f"{q}_mean"] * np.nan)),
                                              to_map(ro.get(f"{q}_q975", ro[f"{q}_mean"] * np.nan)))
                                          for q in ("log_b", "rf_algae", "bba", "f_n", "r_um")},
        }
        # derived quantities have no grid quantiles: build 95 % intervals as mean +/- 1.96 SD
        for nm, rr in (("Tier A Bayesian", ra), ("Physics-informed Bayesian", ro)):
            for q in ("rf_algae", "bba"):
                m_, s_ = to_map(rr[f"{q}_mean"]), to_map(rr[f"{q}_sd"])
                results[nm][q] = (m_, m_ - 1.96 * s_, m_ + 1.96 * s_)
        metrics = VAL.compare(tq, results)
        metrics.to_csv(os.path.join(tab, "validation_metrics.csv"), index=False, float_format="%.5g")
        red_a = VAL.error_reduction(metrics, "Physics-informed Bayesian", "Tier A Bayesian")
        red_e = VAL.error_reduction(metrics, "Physics-informed Bayesian", "Empirical index")
        red_a.to_csv(os.path.join(tab, "error_reduction_vs_tierA.csv"), float_format="%.2f")
        red_e.to_csv(os.path.join(tab, "error_reduction_vs_empirical.csv"), float_format="%.2f")
        print(f"\nValidation against synthetic truth (truth optics: {a.truth_model}):")
        print(metrics.round(3).to_string(index=False))
        print("\nError reduction, physics-informed vs Tier A (%):")
        print(red_a.round(1).to_string())
    if a.field_csv:
        fp = VAL.field_points(a.field_csv, scene.transform, scene.crs, scene.shape)
        for col, mp in (("ours", maps["ours_log_b_mean"]), ("tierA", maps["tierA_log_b_mean"]),
                        ("empirical", logb_emp)):
            fp[f"log_b_{col}"] = mp[fp.row, fp.col]
        fp.to_csv(os.path.join(tab, "field_validation_points.csv"), index=False)
        rows = [dict(method=c, **VAL.metrics(fp.log_b_obs, fp[f"log_b_{c}"])) for c in ("ours", "tierA", "empirical")]
        pd.DataFrame(rows).to_csv(os.path.join(tab, "field_validation_metrics.csv"), index=False)
        print("\nField validation (log10 cells/mL):")
        print(pd.DataFrame(rows).round(3).to_string(index=False))

    # ------------------------------------------------------------- MCMC cross-check
    vidx = np.flatnonzero(valid)
    lbm = ro["log_b_mean"]
    picks = [int(np.nanargmin(np.abs(lbm - np.nanpercentile(lbm, q)))) for q in (25, 85)]
    mcmc_rows, corner_figs = [], []
    lp_all, sk, mu_b_all = prior_logpdfs(em_o.axes, melt.ravel()[valid], slope.ravel()[valid], pc)
    for n_pick, i in enumerate(picks):
        pix = vidx[i]
        pp = dict(mu_b=mu_b_all[i], sd_b=pc.sd_b * pc.scale, f_alpha=pc.f_alpha, f_beta=pc.f_beta,
                  mu_r=pc.r0 + pc.r_melt * melt.ravel()[pix], sd_r=pc.sd_r * pc.scale)
        out = INV.mcmc_pixel(em_o, Rp[pix], sigma, pp, sk, n_steps=a.mcmc_steps, burn=a.mcmc_steps // 3,
                             seed=n_pick)
        lp1 = {k: (v[i:i + 1] if v.shape[0] > 1 else v) for k, v in lp_all.items()}
        g1 = INV.GridPosterior(em_o, sigma).run(Rp[pix][None, :], lp1, sk, keep_full=np.array([0]))
        post = g1["full_posteriors"][0]
        marg = {}
        for ax, n in enumerate(em_o.names):
            other = tuple(j for j in range(post.ndim) if j != ax)
            if len(em_o.axes[n]) > 1:
                marg[n] = (em_o.axes[n], post.sum(axis=other))
        samples, names = out["samples"], list(out["names"])
        pg = em_o.meta["pg_per_cell"]
        B = 10 ** samples[:, names.index("log_b")]
        fn = samples[:, names.index("f_n")]
        samples = np.column_stack([samples, np.log10(B * (fn * pg["nordenskioeldii"] +
                                                         (1 - fn) * pg["alaskanum"]) * 1e-3)])
        names.append("log_pig")
        diag = dict(rhat=np.append(out["rhat"], np.nan), ess=np.append(out["ess"], np.nan))
        tr = None
        if truth is not None:
            tr = {"log_b": truth["log_b"].ravel()[pix], "f_n": truth["f_n"].ravel()[pix],
                  "r_um": truth["r_um"].ravel()[pix], "k": truth["k"].ravel()[pix],
                  "log_pig": np.log10(truth["pigment_ug_l"].ravel()[pix])}
        rr, cc = divmod(pix, W)
        corner_figs.append((M.corner(samples, names, marg, tr, diag,
                                     title=f"Pixel ({rr}, {cc}): MCMC (histograms) vs exact grid posterior (lines)"),
                            f"Fig4C_corner_pixel{n_pick + 1}"))
        mcmc_rows.append(dict(pixel=int(pix), row=rr, col=cc, acceptance=out["acceptance"],
                              **{f"rhat_{n}": v for n, v in zip(out["names"], out["rhat"])},
                              **{f"ess_{n}": v for n, v in zip(out["names"], out["ess"])},
                              **{f"mcmc_mean_{n}": samples[:, j].mean() for j, n in enumerate(names)},
                              **{f"grid_mean_{n}": g1[f"{n}_mean"][0] for n in em_o.names}))
    pd.DataFrame(mcmc_rows).to_csv(os.path.join(tab, "mcmc_diagnostics.csv"), index=False, float_format="%.5g")
    print("\nMCMC cross-check (R-hat should be < 1.05; MCMC and grid means should agree):")
    print(pd.DataFrame(mcmc_rows).filter(regex="rhat|mean_log_b|mean_f_n|acceptance").round(3).to_string(index=False))

    # ------------------------------------------------------------- figures
    M.set_style(a.usetex)
    g = M.MapGrid(scene.transform, scene.crs, scene.shape)
    note = "DEMO input - not Phase 1 results" if a.demo else ""
    rgb = None
    if truth is None:
        rgb = np.clip(np.stack([scene.bands[b] for b in ("B4", "B3", "B2")], -1) / 0.8, 0, 1)
        rgb = np.where(np.isfinite(rgb), rgb, 0.85)
    fig = M.fig4a(g, {"empirical index": logb_emp,
                      "Tier A Bayesian": maps["tierA_log_b_mean"],
                      "ours (physics-informed)": maps["ours_log_b_mean"]},
                  maps["ours_log_b_q975"] - maps["ours_log_b_q025"], maps["ours_f_n_mean"],
                  truth=None if truth is None else truth["log_b"], rgb=rgb, demo_note=note, f_info=shrink_f)
    M.save(fig, figd, "Fig4A_biomass_maps")
    fig = M.fig4b(g, maps["ours_bba_mean"], maps["ours_rf_algae_mean"], maps["d_rf_ours_minus_A"],
                  maps["ours_bba_sd"], maps["ours_rf_algae_sd"], log_bf,
                  r"$\ln$ Bayes factor (ours vs. Tier A)", demo_note=note)
    M.save(fig, figd, "Fig4B_albedo_forcing_maps")
    for fig, name in corner_figs:
        M.save(fig, figd, name)
    if metrics is not None:
        M.save(M.fig_s4_validation(truth["log_b"], {"Empirical index": logb_emp,
                                                    "Tier A Bayesian": maps["tierA_log_b_mean"],
                                                    "Physics-informed Bayesian": maps["ours_log_b_mean"]}, metrics),
               figd, "FigS4_validation_scatter")

    summary = dict(source=a.source, scene=scene.item.get("id", a.local_dir), date=scene.date, sza=scene.sza,
                   shape=[H, W], resolution_m=a.resolution, valid_pixels=int(valid.sum()),
                   empirical_coefs=list(map(float, coef)), f_n_identifiability_sd_ratio=float(shrink_f),
                   chi2_fail_frac={k: float(np.mean(v["chi2"] > CHI2_99_DF4)) for k, v in res.items()},
                   median_log_bayes_factor=float(np.nanmedian(log_bf)),
                   median_rf_ours=float(np.nanmedian(maps["ours_rf_algae_mean"])),
                   median_rf_tierA=float(np.nanmedian(maps["tierA_rf_algae_mean"])),
                   args=vars(a), runtime_s=round(time.time() - t0, 1))
    with open(os.path.join(a.outdir, "summary.json"), "w") as fh:
        json.dump(summary, fh, indent=2, default=str)
    print(f"\nMedian ln Bayes factor (ours vs Tier A): {summary['median_log_bayes_factor']:.2f}")
    print(f"Median RF_algae: ours {summary['median_rf_ours']:.1f}, Tier A {summary['median_rf_tierA']:.1f} W m^-2")
    print(f"Done in {time.time() - t0:.0f} s -> {a.outdir}")


if __name__ == "__main__":
    main()
