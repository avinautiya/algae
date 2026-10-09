#!/usr/bin/env python
"""
Seasonal series of Sentinel-2 retrievals over the S6 dark-zone AOI (2019), so the map-based forcing is
not a single-day number. Runs run_phase4.py (same physics, priors and field calibration) on several
clear scenes and aggregates each map:

  per scene: median and 5-95 % range over bare-ice pixels of log10 B, BBA, instantaneous algal forcing
             at overpass (ours and Tier A), plus a daily-mean forcing estimate
  daily mean: the map forcing is (albedo reduction) x (model clear-sky SW at the overpass SZA), so
             RF_daily = [RF_overpass / SW_model(overpass)] x mean_24h(SW_measured at PROMICE KAN_M).
             The albedo reduction is held at its overpass value (it varies < 10 % with solar zenith over
             the day: literature_comparison.py). Earlier versions divided by the MEASURED SW at a hard-coded
             15 UTC instead of the model SW (column *_legacy_formula keeps that value for comparison).
             Potential melt = RF_daily x 86400 / 3.34e6 cm w.e. d^-1: an upper bound with all extra energy
             to melt, not a surface-energy-balance melt.
  aggregate: a mean over the sampled clear-sky dates, NOT a seasonal mean.
  completion: a scene run counts only with its COMPLETE.json (output checksums).
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


def _promice_day(date):
    d = pd.read_csv(ED.path("promice_KAN_M_hour_JJA_radiation.csv"), parse_dates=["time"])
    d = d[d.time.dt.strftime("%Y-%m-%d") == date].set_index("time").sort_index()
    if len(d) < 20:
        return None
    full = pd.date_range(pd.Timestamp(date), periods=24, freq="h")
    # tilt-corrected SW where available; 2019 has only the uncorrected 'dsr'
    col = d.dsr_cor if d.dsr_cor.notna().sum() >= 20 else d.dsr
    return col.reindex(full).interpolate(limit_direction="both").clip(lower=0).to_numpy()


def sw_at(date, when_utc, label_offset_h: float = 0.5):
    """(daily mean SW, SW interpolated to the overpass time). PROMICE hourly values are averages over the
    hour starting at the time stamp, so they are placed at the hour centre (label_offset_h = 0.5);
    label_offset_h = 1.0 gives the hour-ending convention (sensitivity)."""
    sw = _promice_day(date)
    if sw is None:
        return np.nan, np.nan
    t = pd.Timestamp(when_utc)
    t = t.tz_convert(None) if t.tzinfo else t
    h = (t - pd.Timestamp(date)).total_seconds() / 3600.0
    return float(sw.mean()), float(np.interp(h, np.arange(24) + label_offset_h, sw))


def sw_model(summary):
    """Irradiance used for the emulator forcing: recorded, or recomputed exactly as the emulator does
    (clear-sky at the integer SZA node and day of year)."""
    v = summary.get("sw_down_model_w_m2")
    if v is not None and np.isfinite(v):
        return float(v)
    import biosnicar_bridge as bb
    return bb.sw_down_clear_sky(round(summary["sza"]), ED.clear_sky_transmissivity()[0],
                                int(pd.Timestamp(summary["date"]).dayofyear))


def run_complete(run_dir, args_fp=None):
    """A scene run is consumed only through its COMPLETE.json (artifact checksums; args if given)."""
    import provenance as PV
    p = os.path.join(run_dir, "COMPLETE.json")
    if not os.path.isfile(p):
        return False, "no COMPLETE.json"
    rec = json.load(open(p))
    for name, sha in rec.get("artifacts", {}).items():
        f = os.path.join(run_dir, name)
        if not os.path.isfile(f) or PV.file_sha256(f) != sha:
            return False, f"{name} missing or changed"
    if args_fp is not None and rec.get("args_fingerprint") != args_fp:
        return False, "run settings differ"
    return True, "ok"


def summarise(run_dir, overpass_fallback="15:12"):
    import rasterio
    s = json.load(open(os.path.join(run_dir, "summary.json")))
    with rasterio.open(os.path.join(run_dir, "geotiff", "phase4_maps.tif")) as src:
        m = {src.descriptions[i]: src.read(i + 1) for i in range(src.count)}
    valid = np.isfinite(m["ours_log_b_mean"])
    ice = valid & (m["ours_bba_mean"] <= SNOW_BBA)
    when = s.get("overpass_datetime") or ED_overpass(s["scene"]) or f"{s['date']} {overpass_fallback}"
    sw_mean, sw_ov = sw_at(s["date"], when)
    sw_mod = sw_model(s)
    row = dict(scene=s["scene"], date=s["date"], overpass_utc=str(when), sza=s["sza"], valid_pixels=int(valid.sum()),
               snow_excluded_frac=float(1 - ice.sum() / max(valid.sum(), 1)), sw_daily_mean=sw_mean,
               sw_measured_overpass=sw_ov, sw_model_overpass=sw_mod,
               clearness_overpass=sw_ov / sw_mod if sw_mod else np.nan,
               tau=s.get("model_error_tau_dex", {}).get("ours", np.nan))
    for k, lab in (("ours_log_b_mean", "log_b"), ("ours_bba_mean", "bba"), ("ours_rf_algae_mean", "rf_ours"),
                   ("tierA_rf_algae_mean", "rf_tierA"), ("tierA_log_b_mean", "log_b_tierA"),
                   ("tierA_bba_mean", "bba_tierA")):
        if k not in m:
            continue
        v = m[k][ice]
        row.update({f"{lab}_median": float(np.nanmedian(v)), f"{lab}_p5": float(np.nanpercentile(v, 5)),
                    f"{lab}_p95": float(np.nanpercentile(v, 95)), f"{lab}_mean": float(np.nanmean(v))})
    for lab in ("rf_ours", "rf_tierA"):
        # the map forcing is (albedo reduction) x (model clear-sky SW at overpass): divide by THAT SW to get
        # the albedo reduction, then multiply by the measured 24 h mean SW (albedo reduction held at its
        # overpass value; its diurnal SZA dependence is the stated approximation, see literature_comparison)
        dalpha = row[f"{lab}_mean"] / sw_mod
        row[f"dalpha_{lab}"] = dalpha
        row[f"{lab}_daily_mean"] = dalpha * sw_mean
        row[f"{lab}_daily_mean_legacy_formula"] = row[f"{lab}_mean"] * sw_mean / sw_ov if sw_ov else np.nan
        # POTENTIAL melt: all extra absorbed energy assumed to go into melt (no SEB, no cold content,
        # no refreezing); an upper bound, not a melt estimate
        row[f"melt_potential_{lab}_cm_we_d"] = row[f"{lab}_daily_mean"] * 86400.0 / L_F
    return row


def ED_overpass(scene_id):
    """Overpass time of a scene from its STAC item (network); None if unavailable."""
    try:
        import s2_io as io
        return io.item_from_s3(scene_id, retries=1)["properties"]["datetime"]
    except Exception:  # noqa: BLE001
        return None


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
    p.add_argument("--accept-legacy", action="store_true",
                   help="with --aggregate-only: also aggregate runs made before completion markers existed "
                        "(reported per scene as legacy)")
    a = p.parse_args(argv)
    os.makedirs(a.outdir, exist_ok=True)
    a.outdir = os.path.abspath(a.outdir)
    cache = os.path.join(a.outdir, "shared_cache")      # absolute: symlinks below must resolve from any folder
    rows = []
    for sid in a.scenes:
        out = os.path.join(a.outdir, sid)
        done, why = run_complete(out)
        if not a.aggregate_only and not done:
            print(f"{sid}: running ({why})", flush=True)
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
        done, why = run_complete(out)
        if done or (a.aggregate_only and a.accept_legacy and os.path.isfile(os.path.join(out, "summary.json"))):
            rows.append(summarise(out))
            rows[-1]["completion"] = "COMPLETE.json" if done else f"legacy (no marker: {why})"
        else:
            print(f"{sid}: skipped from the aggregate ({why})", flush=True)
    if not rows:
        raise SystemExit("no completed scene runs")
    tab = pd.DataFrame(rows).sort_values("date")
    import provenance as PV
    PV.atomic_to_csv(tab, os.path.join(a.outdir, "multi_scene_summary.csv"), index=False, float_format="%.4g")
    pd.set_option("display.width", 220)
    cols = ["date", "valid_pixels", "snow_excluded_frac", "log_b_median", "bba_median", "rf_ours_median",
            "rf_tierA_median", "clearness_overpass", "rf_ours_daily_mean", "rf_ours_daily_mean_legacy_formula",
            "melt_potential_rf_ours_cm_we_d", "melt_potential_rf_tierA_cm_we_d"]
    print(tab[[c for c in cols if c in tab]].round(3).to_string(index=False))
    if len(tab) > 1:
        # a mean over the SAMPLED clear-sky dates, not a seasonal mean (cloudy days are not represented)
        season = dict(n_dates=len(tab), dates=tab.date.tolist(),
                      rf_ours_daily_mean_sampled_dates=float(tab.rf_ours_daily_mean.mean()),
                      rf_ours_daily_mean_range=[float(tab.rf_ours_daily_mean.min()), float(tab.rf_ours_daily_mean.max())],
                      melt_potential_ours_mean_cm_we_d=float(tab.melt_potential_rf_ours_cm_we_d.mean()),
                      rf_tierA_daily_mean_sampled_dates=float(tab.rf_tierA_daily_mean.mean()),
                      note="mean over clear-sky sampled dates only; not a seasonal mean; melt is potential "
                           "(no surface energy balance)")
        PV_write(os.path.join(a.outdir, "multi_scene_season.json"), season)
        print(json.dumps(season, indent=1))
        figure(tab, os.path.join(a.outdir, "FigS7_seasonal_series.png"))


def PV_write(path, obj):
    import provenance as PV
    PV.atomic_write_text(path, json.dumps(obj, indent=1, default=float))


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
