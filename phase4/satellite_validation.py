#!/usr/bin/env python
"""
Satellite-scale validation: Sentinel-2 pixel retrievals against the counted S6 samples of July 2017.

The field validation (field_validation.py) tests the retrieval on ground spectra under the sampled
ice (plot scale, ~0.5 m footprint). This script applies the SAME retrieval (emulator, priors, sigma,
ice-radius prior and model error tau from the field calibration) to Sentinel-2 L2A pixels and compares
the pixel estimate with the cell counts of the samples that fall in it:

  samples    20 counted S6 samples with GPS positions (15, 21, 22, 23 July 2017): counts from the Cook
             et al. (2020) archive, positions from Tedstone et al. (2020) UK PDC 10.5285/77ca631f-a3a4-
             4f26-bc90-57bb17baa6fc (uav_sb_locations.csv, UTM 23N / EPSG:32623 - the Sentinel-2 tile
             22WEV is in UTM 22N, so positions are reprojected)
  scenes     primary: the nearest clear scene (tile cloud < 10 %): 11 Jul for the 15 Jul samples (-4 d),
             21 Jul for the 21-23 Jul samples (0 to -2 d); sensitivity: the same-day scenes 15 Jul
             (tile cloud 34 %) and 23 Jul (43 %, thin cirrus near the site, masked by SCL)
  masking    Sen2Cor SCL class 11 (snow/ice) only, as in the maps

Scale mismatch. A sample is a ~1 m^2 patch of the top 2 cm; a pixel is 100 m^2, and the 2017 Sentinel-2
geolocation error is ~10 m (before the 2019 global reference image). Each sample is therefore compared
with (a) the pixel containing it and (b) the mean over the 3x3 pixels around it. The irreducible part of
the disagreement - plot-to-plot variability inside one pixel - is estimated from samples sharing a
pixel (pooled SD of their log counts), and the site-level comparison (mean over all samples of a day vs
mean over the pixels they fall in) removes most of it.

    python phase4/satellite_validation.py --phenol williamson2020 --outdir phase4/results/satellite_validation
"""

from __future__ import annotations

import argparse
import dataclasses
import json
import os
import sys
import time

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "phase2"))

import empirical_data as ED  # noqa: E402

SAMPLE_CRS = 32623
# scene -> sample days (day_month prefixes of the sample IDs) it is compared with
SCENES = {
    "S2B_22WEV_20170711_0_L2A": ("15_7",),
    "S2B_22WEV_20170721_0_L2A": ("21_7", "22_7", "23_7"),
    "S2B_22WEV_20170715_0_L2A": ("15_7",),
    "S2A_22WEV_20170723_0_L2A": ("21_7", "22_7", "23_7"),
}
PRIMARY = ("S2B_22WEV_20170711_0_L2A", "S2B_22WEV_20170721_0_L2A")
HALF_WINDOW_M = 150.0


def samples():
    """Counted S6 samples with positions: sample, day, x/y (UTM 23N), cells, cells_counted, log_b_obs."""
    loc = pd.read_csv(ED.path("tedstone2020_s6_2017_sample_locations.csv"))
    cnt = pd.read_csv(ED.path("cook2020_archive_cell_counts.csv"))
    d = loc.merge(cnt[["sample", "cells_per_ml", "cells_counted"]], on="sample", how="inner")
    d["day"] = d["sample"].str.rsplit("_", n=1).str[0]
    d = d.rename(columns={"cells_per_ml": "cells"})
    d["log_b_obs"] = np.where(d.cells > 0, np.log10(d.cells.clip(lower=1e-9)), np.nan)
    d["log_b_obs_sd"] = 1.0 / np.log(10) / np.sqrt(d.cells_counted.clip(lower=1))
    return d


def invert_scene(scene_id, pts, cal, a, cache):
    """Read a 300 m window around the samples, invert it with the field-calibrated retrieval, and return
    (per-sample table, maps dict, scene)."""
    import emulator as E
    import inversion as INV
    import s2_io as io
    from pyproj import Transformer
    from priors import PriorConfig, prior_logpdfs
    from run_phase4 import parse_args

    d = parse_args([])                                   # same physics defaults as the maps
    item = io.item_from_s3(scene_id)
    epsg = int(str(item["properties"].get("proj:epsg") or item["properties"]["proj:code"]).split(":")[-1])
    tr = Transformer.from_crs(SAMPLE_CRS, epsg, always_xy=True)
    x, y = tr.transform(pts.x_utm23n.to_numpy(), pts.y_utm23n.to_numpy())
    # snap the window to the 10 m Sentinel-2 grid (tile origins are multiples of 10 m) so pixels are native
    cx, cy = np.round(np.mean(x), -1), np.round(np.mean(y), -1)
    bounds = (cx - HALF_WINDOW_M, cy - HALF_WINDOW_M, cx + HALF_WINDOW_M, cy + HALF_WINDOW_M)
    sc = io.read_scene(item, bounds, resolution=10.0)
    spacecraft = ED.spacecraft_from_scene(scene_id)
    cfg = E.EmulatorConfig(model="ours", tier=a.tier, phenol=a.phenol, photosynthetic=True,
                           f_n=(0.0, 1.0, d.f_step), ice_mode=d.ice_mode, rho=d.rho, rho_bottom=d.rho_bottom,
                           sza=round(sc.sza), r_um=tuple(d.r_range), spacecraft=spacecraft,
                           dust_ppb=E.DUST_NODES_PPB, day_of_year=int(pd.Timestamp(sc.date).dayofyear))
    em = E.build_emulator(cfg, phase1_l2=a.phase1_l2, biosnicar=a.biosnicar, workers=a.workers,
                          cache=os.path.join(cache, f"emulator_{a.phenol}_{scene_id}.npz")).refine_log_b(d.log_b_step)
    pc = dataclasses.replace(PriorConfig.for_density(d.rho_bottom), mu_lnr=cal["mu_lnr"], sd_lnr=cal["sd_lnr"])
    R = sc.stack().reshape(-1, 4)
    valid = sc.mask.ravel() & np.all(np.isfinite(R), axis=1)
    lp, sk, mk, _ = prior_logpdfs(em.axes, int(valid.sum()), pc)
    r = INV.GridPosterior(em, np.full(4, cal["sigma"])).run(R[valid], lp, sk, mk=mk)
    H, W = sc.shape

    def to_map(v):
        out = np.full(H * W, np.nan)
        out[valid] = v
        return out.reshape(H, W)

    maps = {k: to_map(r[k]) for k in ("log_b_mean", "log_b_sd", "log_b_q025", "log_b_q975", "chi2",
                                      "rf_algae_mean", "bba_mean")}
    maps["B2"] = sc.bands["B2"]
    inv = ~sc.transform
    rows = []
    for (xi, yi), p in zip(zip(x, y), pts.itertuples()):
        col, row = (int(np.floor(v)) for v in inv * (xi, yi))
        blk = (slice(max(row - 1, 0), row + 2), slice(max(col - 1, 0), col + 2))
        m, s = maps["log_b_mean"][row, col], maps["log_b_sd"][row, col]
        half = 1.96 * np.hypot(s, cal["tau"])
        # pixel identity: scene + absolute pixel-centre coordinates (row/col are window-relative and
        # collide between scenes whose windows differ)
        xc, yc = sc.transform * (col + 0.5, row + 0.5)
        rows.append(dict(scene=scene_id, scene_date=sc.date, sample=p.sample, day=p.day, cells=p.cells,
                         log_b_obs=p.log_b_obs, log_b_obs_sd=p.log_b_obs_sd, row=row, col=col,
                         pixel=f"{scene_id}|{xc:.0f}_{yc:.0f}", B2=maps["B2"][row, col], chi2=maps["chi2"][row, col],
                         sat_log_b_mean=m, sat_log_b_sd=s, sat_log_b_q025=maps["log_b_q025"][row, col],
                         sat_log_b_q975=maps["log_b_q975"][row, col], sat_log_b_q025_cal=m - half,
                         sat_log_b_q975_cal=m + half, sat3x3_log_b_mean=np.nanmean(maps["log_b_mean"][blk]),
                         sat_rf_algae=maps["rf_algae_mean"][row, col], sat_bba=maps["bba_mean"][row, col]))
    out = pd.DataFrame(rows)
    out.attrs.update(sza=sc.sza, valid_pixels=int(valid.sum()), n_pixels=H * W,
                     chi2_fail=float(np.mean(r["chi2"] > 13.28)))
    return out, maps, sc


def _metrics(t, e, lo=None, hi=None):
    from scipy.stats import spearmanr
    ok = np.isfinite(t) & np.isfinite(e)
    err = e[ok] - t[ok]
    o = dict(n=int(ok.sum()), bias=float(err.mean()), rmse=float(np.sqrt(np.mean(err ** 2))),
             spearman=float(spearmanr(t[ok], e[ok]).statistic) if ok.sum() > 2 else np.nan)
    if lo is not None:
        o["coverage95"] = float(np.mean((t[ok] >= lo[ok]) & (t[ok] <= hi[ok])))
    return o


def within_pixel_sd(df):
    """Pooled SD of log counts among samples that share a pixel (plot-to-plot variability a pixel
    cannot resolve), with its degrees of freedom and a 95 % chi-square interval. Pixels are identified
    by scene AND absolute position (see invert_scene); zero counts have no log value and are excluded
    (reported by the caller). Pooling assumes one common within-pixel variance."""
    from scipy.stats import chi2
    if df.pixel.str.contains("|", regex=False).mean() < 1:
        raise ValueError("pixel ids must carry the scene and absolute position")
    ss, dof = 0.0, 0
    for _, g in df[np.isfinite(df.log_b_obs)].groupby("pixel"):
        if len(g) > 1:
            ss += float(np.sum((g.log_b_obs - g.log_b_obs.mean()) ** 2))
            dof += len(g) - 1
    if not dof:
        return np.nan, 0, (np.nan, np.nan)
    sd = np.sqrt(ss / dof)
    ci = (float(np.sqrt(ss / chi2.ppf(0.975, dof))), float(np.sqrt(ss / chi2.ppf(0.025, dof))))
    return float(sd), dof, ci


def subpixel_mixing(em, sd_dex, mu=3.5, n=4000, seed=0, r_um=None, f_n=0.5):
    """What a pixel retrieval of an area-mixed surface estimates: draw log10 B ~ N(mu, sd_dex) for the
    sub-pixel patches, average their band reflectances (linear areal mixing), and find the log10 B whose
    reflectance matches the mean (least squares over the emulator's log B axis at fixed radius, f_n and
    the lowest dust node). Returns the retrieved value next to the geometric (mean log) and arithmetic
    (log mean) abundance of the patches."""
    from scipy.interpolate import RegularGridInterpolator
    rng = np.random.default_rng(seed)
    lb = em.axes["log_b"]
    r_um = r_um or float(np.exp(np.mean(np.log(em.axes["r_um"]))))
    pts = []
    for v in lb:
        q = {"log_b": v, "f_n": f_n, "r_um": r_um, "dust_ppb": em.axes.get("dust_ppb", [0])[0]}
        pts.append([q[nm] if len(em.axes[nm]) > 1 else em.axes[nm][0] for nm in em.names])
    I = RegularGridInterpolator(tuple(em.coord(nm) for nm in em.names), em.data["bands"], bounds_error=False,
                                fill_value=None)
    conv = lambda P: np.array([[em.coord(nm)[0] if len(em.axes[nm]) == 1 else  # noqa: E731
                                (np.log10(x + 100.0) if nm == "dust_ppb" else x) for nm, x in zip(em.names, p)]
                               for p in P])
    curve = I(conv(pts))                                      # (n_lb, 4)
    x = np.clip(rng.normal(mu, sd_dex, n), lb[0], lb[-1])
    R = np.array([np.interp(x, lb, curve[:, b]) for b in range(4)]).T.mean(axis=0)
    fine = np.linspace(lb[0], lb[-1], 2001)
    cf = np.array([np.interp(fine, lb, curve[:, b]) for b in range(4)]).T
    retrieved = float(fine[np.argmin(np.sum((cf - R) ** 2, axis=1))])
    return dict(sd_dex=sd_dex, mu=mu, retrieved=retrieved, geometric=float(x.mean()),
                arithmetic=float(np.log10(np.mean(10 ** x))))


def figure(mapdf, df, maps, sc, pts_xy, plot_scale, path):
    """a: retrieval map of one scene with that scene's samples; b: pixel vs count for all rows of df."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sys.path.insert(0, os.path.join(HERE, "..", "phase3"))
    fig, ax = plt.subplots(1, 2, figsize=(7.2, 3.3), gridspec_kw=dict(width_ratios=[1, 1.15]))
    t = sc.transform
    ext = (t.c, t.c + t.a * sc.shape[1], t.f + t.e * sc.shape[0], t.f)
    im = ax[0].imshow(maps["log_b_mean"], extent=ext, cmap="viridis", vmin=3, vmax=5.5)
    ax[0].scatter(*pts_xy, c=mapdf.log_b_obs, cmap="viridis", vmin=3, vmax=5.5, s=16, edgecolor="white", lw=0.6)
    ax[0].set_xticks([]), ax[0].set_yticks([])
    ax[0].plot([ext[0] + 20, ext[0] + 120], [ext[2] + 20] * 2, color="white", lw=2)
    ax[0].text(ext[0] + 70, ext[2] + 28, "100 m", color="white", ha="center", fontsize=7)
    plt.colorbar(im, ax=ax[0], fraction=0.046, label=r"$\log_{10}$ cells mL$^{-1}$")
    ax[0].set_title(f"a  Sentinel-2 retrieval, {sc.date} (points: counts)", loc="left", fontsize=8)
    ok = np.isfinite(df.log_b_obs)
    d = df[ok]
    ax[1].errorbar(d.log_b_obs, d.sat_log_b_mean, yerr=[d.sat_log_b_mean - d.sat_log_b_q025_cal,
                                                         d.sat_log_b_q975_cal - d.sat_log_b_mean],
                   fmt="o", ms=4, color="#1b6ca8", ecolor="#9ec3df", lw=0.8, label="S2 pixel (95 % incl. $\\tau$)")
    ax[1].scatter(d.log_b_obs, d.sat3x3_log_b_mean, marker="s", s=14, color="#e07b39", label="S2 3x3 mean", zorder=3)
    if plot_scale is not None:
        ps = plot_scale.reindex(d["sample"]).to_numpy()
        ax[1].scatter(d.log_b_obs, ps, marker="^", s=16, facecolor="none", edgecolor="k", lw=0.7,
                      label="ground spectrum (plot scale)", zorder=3)
    lim = (2.0, 5.8)
    ax[1].plot(lim, lim, color="0.4", lw=0.8)
    ax[1].set_xlim(lim), ax[1].set_ylim(lim)
    ax[1].set_xlabel(r"counted $\log_{10}$ cells mL$^{-1}$")
    ax[1].set_ylabel(r"retrieved $\log_{10}$ cells mL$^{-1}$")
    ax[1].legend(fontsize=6, loc="upper left", frameon=False)
    ax[1].set_title("b  pixel vs. sample, primary scenes", loc="left", fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=600)
    fig.savefig(path.replace(".png", ".pdf"))
    plt.close(fig)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--phase1-l2", default=os.path.join(HERE, "..", "phase1", "results", "level2"))
    p.add_argument("--phenol", default="tddft", choices=["tddft", "williamson2020"])
    p.add_argument("--tier", default=None, choices=["C", "D"], help="default D for tddft, C otherwise")
    p.add_argument("--biosnicar", default=None)
    p.add_argument("--scenes", nargs="+", default=list(SCENES))
    p.add_argument("--workers", type=int, default=os.cpu_count() or 1)
    p.add_argument("--outdir", default=os.path.join(HERE, "results", "satellite_validation"))
    a = p.parse_args(argv)
    a.tier = a.tier or ("D" if a.phenol == "tddft" else "C")
    t0 = time.time()
    cache = os.path.join(a.outdir, "cache")
    os.makedirs(cache, exist_ok=True)

    # ---------------------------------------------------------------- field calibration (no counts used)
    import field_validation as FV
    # the samples compared with pixels below are EXCLUDED from the field calibration (sigma, radius prior,
    # tau), so the satellite comparison is held out (earlier versions fitted tau on all S6 counts)
    pts_all = samples()
    fdf, fmet, fres = FV.run(phase1_l2=a.phase1_l2, biosnicar=a.biosnicar, workers=a.workers, cache_dir=cache,
                             tier=a.tier, phenol=a.phenol, models=("ours",), dust=True, verbose=False,
                             exclude=set(pts_all["sample"]))
    r = fres["ours"]
    cal = dict(sigma=r["sigma_all"], mu_lnr=r["mu_lnr"], sd_lnr=r["sd_lnr"], tau=r["tau"])
    print(f"Field calibration: sigma {cal['sigma']:.3f}, ice radius median {np.exp(cal['mu_lnr']):.0f} um "
          f"(ln-SD {cal['sd_lnr']:.2f}), model error tau {cal['tau']:.2f} dex", flush=True)
    plot_scale = fdf[fdf.dataset == "s6_2017"].set_index("sample")["ours_log_b_mean"]

    pts = samples()
    print(f"{len(pts)} counted samples with GPS positions; days {sorted(set(pts.day))}", flush=True)
    tabs, keep = [], {}
    for sid in a.scenes:
        sub = pts[pts.day.isin(SCENES[sid])].reset_index(drop=True)
        t1 = time.time()
        tab, maps, sc = invert_scene(sid, sub, cal, a, cache)
        tab["primary"] = sid in PRIMARY
        tab["plot_log_b_mean"] = plot_scale.reindex(tab["sample"]).to_numpy()
        tabs.append(tab)
        keep[sid] = (maps, sc, sub)
        print(f"[{sid}] SZA {sc.sza:.1f}, {tab.attrs['valid_pixels']}/{tab.attrs['n_pixels']} snow/ice pixels, "
              f"chi2 fail {100 * tab.attrs['chi2_fail']:.0f} %, {time.time() - t1:.0f} s", flush=True)
    df = pd.concat(tabs, ignore_index=True)
    os.makedirs(a.outdir, exist_ok=True)
    df.to_csv(os.path.join(a.outdir, "satellite_validation_samples.csv"), index=False, float_format="%.5g")

    # ---------------------------------------------------------------- metrics
    rows = []
    for name, sel in [("primary (nearest clear scene)", df.primary)] + [(s, df.scene == s) for s in a.scenes]:
        d = df[sel & np.isfinite(df.log_b_obs)]
        t = d.log_b_obs.to_numpy()
        for est, lo, hi in (("sat_log_b_mean", "sat_log_b_q025_cal", "sat_log_b_q975_cal"),
                            ("sat3x3_log_b_mean", None, None), ("plot_log_b_mean", None, None)):
            rows.append(dict(scenes=name, estimate=est, **_metrics(
                t, d[est].to_numpy(), None if lo is None else d[lo].to_numpy(), None if hi is None else d[hi].to_numpy())))
        # site level: mean log count of each day's samples vs mean retrieval over the pixels they occupy
        # E[log B] (geometric) and log E[B] (arithmetic, zeros included) are both reported: a pixel sees
        # area-mixed reflectance, which corresponds to neither exactly (see subpixel_mixing)
        dall = df[sel]
        for day, g in d.groupby("day"):
            px = g.drop_duplicates("pixel")
            ga = dall[dall.day == day]
            arith = float(np.log10(ga.cells.mean()))
            rows.append(dict(scenes=name, estimate=f"site mean, {day}", n=len(g), n_zero=int((ga.cells <= 0).sum()),
                             bias=float(px.sat_log_b_mean.mean() - g.log_b_obs.mean()),
                             bias_vs_arithmetic=float(px.sat_log_b_mean.mean() - arith),
                             obs_mean=float(g.log_b_obs.mean()), obs_log_arith_mean=arith,
                             obs_se=float(g.log_b_obs.std(ddof=1) / np.sqrt(len(g))),
                             sat_mean=float(px.sat_log_b_mean.mean()), n_pixels=len(px)))
    met = pd.DataFrame(rows)
    met.to_csv(os.path.join(a.outdir, "satellite_validation_metrics.csv"), index=False, float_format="%.4g")
    prim = df[df.primary]
    sd_wp, dof, sd_ci = within_pixel_sd(prim)
    summary = dict(calibration=cal, within_pixel_sd_log_counts=sd_wp, within_pixel_dof=dof, within_pixel_sd_ci95=sd_ci,
                   n_samples=int(np.isfinite(prim.log_b_obs).sum()), n_pixels=int(prim.pixel.nunique()),
                   scenes={s: dict(sza=keep[s][1].sza, date=keep[s][1].date) for s in keep},
                   phenol=a.phenol, runtime_s=round(time.time() - t0))
    with open(os.path.join(a.outdir, "satellite_validation_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1, default=float)
    pd.set_option("display.width", 200)
    print(met.round(3).to_string(index=False))
    print(f"Plot-to-plot SD of log counts within one 10 m pixel: {sd_wp:.2f} dex ({dof} dof) - the RMSE floor "
          "for any pixel-vs-sample comparison")

    # ---------------------------------------------------------------- figure (21 Jul scene, all primary points)
    from pyproj import Transformer
    maps, sc, sub = keep.get(PRIMARY[1], next(iter(keep.values())))
    epsg = sc.crs.to_epsg()
    tr = Transformer.from_crs(SAMPLE_CRS, epsg, always_xy=True)
    xy = tr.transform(sub.x_utm23n.to_numpy(), sub.y_utm23n.to_numpy())
    figure(df[df.scene == PRIMARY[1]].reset_index(drop=True), df[df.primary], maps, sc, xy, plot_scale, os.path.join(a.outdir, "FigS6_satellite_validation.png"))
    print(f"Done in {time.time() - t0:.0f} s -> {a.outdir}")


if __name__ == "__main__":
    main()
