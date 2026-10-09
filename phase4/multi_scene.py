#!/usr/bin/env python
"""
Seasonal series of Sentinel-2 retrievals over the S6 dark-zone AOI (2019), so the map-based forcing is
not a single-day number. Runs run_phase4.py (same physics, priors and field calibration) on several
clear scenes and aggregates each map:

  per scene: median and 5-95 % range over bare-ice pixels of log10 B, BBA, instantaneous algal forcing
             at overpass (ours and Tier A), plus a daily-mean forcing estimate
  daily mean: RF_daily = RF_overpass x mean_24h(SW) / SW(overpass hour), with the measured hourly SW at
             PROMICE KAN_M on that day - i.e. the albedo reduction is held at its overpass value (it varies
             < 10 % with solar zenith over the day: literature_comparison.py) and scaled with the real
             irradiance cycle, clouds included; melt potential = RF_daily x 86400 / 3.34e6 cm w.e. d^-1
  bare ice:  Sen2Cor SCL class 11 does not separate snow from ice; pixels with retrieved BBA > 0.65
             (fresh/old snow; bare ice is below ~0.6 even when clean - field CI 0.53 +/- 0.03, Cook et al.
             2020 Table 2) are excluded, and the excluded fraction is reported.

Scenes (tile 22WEV, cloud < 5 %, no nodata over the AOI; found with s2_io.search_scenes):
  08 Jul, 15 Jul, 23 Jul, 02 Aug, 12 Aug, 29 Aug 2019

    python phase4/multi_scene.py --phenol tddft --outdir phase4/results/multi_scene
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "phase2"))

import empirical_data as ED  # noqa: E402

SCENES = ["S2B_22WEV_20190708_0_L2A", "S2B_22WEV_20190715_0_L2A", "S2A_22WEV_20190723_0_L2A",
          "S2A_22WEV_20190802_0_L2A", "S2A_22WEV_20190812_0_L2A", "S2A_22WEV_20190829_0_L2A"]
SNOW_BBA = 0.65
L_F = 3.34e6                      # J per cm w.e. per m^2


def sw_ratio(date, hour_utc):
    """mean_24h(SW) / SW(overpass hour) from PROMICE KAN_M hourly dsr_cor (NaN if the day is incomplete)."""
    d = pd.read_csv(ED.path("promice_KAN_M_hour_JJA_radiation.csv"), parse_dates=["time"])
    d = d[d.time.dt.strftime("%Y-%m-%d") == date].set_index("time").sort_index()
    if len(d) < 20:
        return np.nan, np.nan
    full = pd.date_range(pd.Timestamp(date), periods=24, freq="h")
    # tilt-corrected SW where available; 2019 has only the uncorrected 'dsr' (the ratio below is
    # insensitive to a tilt calibration that scales the whole day)
    col = d.dsr_cor if d.dsr_cor.notna().sum() >= 20 else d.dsr
    sw = col.reindex(full).interpolate(limit_direction="both").clip(lower=0).to_numpy()
    return float(sw.mean() / sw[int(hour_utc)]), float(sw.mean())


def summarise(run_dir):
    import rasterio
    s = json.load(open(os.path.join(run_dir, "summary.json")))
    with rasterio.open(os.path.join(run_dir, "geotiff", "phase4_maps.tif")) as src:
        m = {src.descriptions[i]: src.read(i + 1) for i in range(src.count)}
    valid = np.isfinite(m["ours_log_b_mean"])
    ice = valid & (m["ours_bba_mean"] <= SNOW_BBA)
    ratio, sw_mean = sw_ratio(s["date"], 15)     # tile 22WEV is imaged at 15:00-15:25 UTC (STAC datetime)
    row = dict(scene=s["scene"], date=s["date"], sza=s["sza"], valid_pixels=int(valid.sum()),
               snow_excluded_frac=float(1 - ice.sum() / max(valid.sum(), 1)), sw_daily_mean=sw_mean,
               sw_ratio_daily_to_overpass=ratio, tau=s.get("model_error_tau_dex", {}).get("ours", np.nan))
    for k, lab in (("ours_log_b_mean", "log_b"), ("ours_bba_mean", "bba"), ("ours_rf_algae_mean", "rf_ours"),
                   ("tierA_rf_algae_mean", "rf_tierA"), ("tierA_log_b_mean", "log_b_tierA")):
        v = m[k][ice]
        row.update({f"{lab}_median": float(np.nanmedian(v)), f"{lab}_p5": float(np.nanpercentile(v, 5)),
                    f"{lab}_p95": float(np.nanpercentile(v, 95)), f"{lab}_mean": float(np.nanmean(v))})
    for lab in ("rf_ours", "rf_tierA"):
        row[f"{lab}_daily_mean"] = row[f"{lab}_mean"] * ratio
        row[f"melt_{lab}_cm_we_d"] = row[f"{lab}_daily_mean"] * 86400.0 / L_F
    return row


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--scenes", nargs="+", default=SCENES)
    p.add_argument("--phenol", default="tddft", choices=["tddft", "williamson2020"])
    p.add_argument("--phase1-l2", default=os.path.join(HERE, "..", "phase1", "results", "level2"))
    p.add_argument("--tier", default=None, choices=["C", "D"], help="default D for tddft, C otherwise")
    p.add_argument("--resolution", type=float, default=20.0)
    p.add_argument("--workers", type=int, default=os.cpu_count() or 1)
    p.add_argument("--outdir", default=os.path.join(HERE, "results", "multi_scene"))
    p.add_argument("--aggregate-only", action="store_true")
    a = p.parse_args(argv)
    os.makedirs(a.outdir, exist_ok=True)
    a.outdir = os.path.abspath(a.outdir)
    cache = os.path.join(a.outdir, "shared_cache")      # absolute: symlinks below must resolve from any folder
    rows = []
    for sid in a.scenes:
        out = os.path.join(a.outdir, sid)
        if not a.aggregate_only and not os.path.isfile(os.path.join(out, "summary.json")):
            os.makedirs(os.path.join(out, "cache"), exist_ok=True)
            # share field-validation emulators between scenes (they do not depend on the scene)
            for f in os.listdir(cache) if os.path.isdir(cache) else []:
                dst = os.path.join(out, "cache", f)
                if f.startswith("field_emulator_") and not os.path.exists(dst):
                    os.symlink(os.path.join(cache, f), dst)
            cmd = [sys.executable, os.path.join(HERE, "run_phase4.py"), "--source", "s2", "--scene-id", sid,
                   "--phenol", a.phenol, "--tier", a.tier or ("D" if a.phenol == "tddft" else "C"), "--phase1-l2", a.phase1_l2, "--resolution", str(a.resolution),
                   "--workers", str(a.workers), "--no-dem", "--mcmc-steps", "1500", "--outdir", out]
            print(" ".join(cmd), flush=True)
            subprocess.run(cmd, check=True)
            os.makedirs(cache, exist_ok=True)
            for f in os.listdir(os.path.join(out, "cache")):
                src = os.path.join(out, "cache", f)
                if f.startswith("field_emulator_") and not os.path.islink(src) and not os.path.exists(os.path.join(cache, f)):
                    os.replace(src, os.path.join(cache, f))
                    os.symlink(os.path.join(cache, f), src)
        if os.path.isfile(os.path.join(out, "summary.json")):
            rows.append(summarise(out))
    tab = pd.DataFrame(rows).sort_values("date")
    tab.to_csv(os.path.join(a.outdir, "multi_scene_summary.csv"), index=False, float_format="%.4g")
    pd.set_option("display.width", 220)
    cols = ["date", "valid_pixels", "snow_excluded_frac", "log_b_median", "bba_median", "rf_ours_median",
            "rf_tierA_median", "rf_ours_daily_mean", "melt_rf_ours_cm_we_d", "melt_rf_tierA_cm_we_d"]
    print(tab[[c for c in cols if c in tab]].round(3).to_string(index=False))
    if len(tab) > 1:
        season = dict(n_scenes=len(tab), rf_ours_daily_mean_season=float(tab.rf_ours_daily_mean.mean()),
                      rf_ours_daily_mean_range=[float(tab.rf_ours_daily_mean.min()), float(tab.rf_ours_daily_mean.max())],
                      melt_ours_mean_cm_we_d=float(tab.melt_rf_ours_cm_we_d.mean()),
                      rf_tierA_daily_mean_season=float(tab.rf_tierA_daily_mean.mean()))
        json.dump(season, open(os.path.join(a.outdir, "multi_scene_season.json"), "w"), indent=1)
        print(json.dumps(season, indent=1))
        figure(tab, os.path.join(a.outdir, "FigS7_seasonal_series.png"))


def figure(tab, path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    d = pd.to_datetime(tab.date)
    fig, ax = plt.subplots(1, 3, figsize=(7.2, 2.4))
    ax[0].fill_between(d, tab.log_b_p5, tab.log_b_p95, color="#9ec3df")
    ax[0].plot(d, tab.log_b_median, "o-", color="#1b6ca8", ms=3)
    ax[0].set_ylabel(r"$\log_{10}$ cells mL$^{-1}$")
    ax[1].fill_between(d, tab.bba_p5, tab.bba_p95, color="#d9d9d9")
    ax[1].plot(d, tab.bba_median, "o-", color="0.3", ms=3)
    ax[1].set_ylabel("broadband albedo")
    ax[2].plot(d, tab.rf_ours_daily_mean, "o-", color="#1b6ca8", ms=3, label="ours")
    ax[2].plot(d, tab.rf_tierA_daily_mean, "s--", color="#e07b39", ms=3, label="Tier A")
    ax[2].set_ylabel(r"daily-mean algal RF (W m$^{-2}$)")
    ax[2].legend(fontsize=6, frameon=False)
    for i, x in enumerate(ax):
        x.tick_params(axis="x", labelrotation=45, labelsize=6)
        x.set_title("abc"[i], loc="left", fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=600)
    fig.savefig(path.replace(".png", ".pdf"))
    plt.close(fig)


if __name__ == "__main__":
    main()
