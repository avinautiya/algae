#!/usr/bin/env python
"""
H5 (protocol amendment A2/A2.1): station broadband albedo near solar noon, predicted from the ice-algae
count measured UNDER the station radiometer (PROMBIO), optics-only with known abundance.

  * Models: M0 no algae, M1 tierA_empirical, M1b measured_mac_C, M3 tddft_D, M3-alt tddft_C/_iid.
    Frozen held-out emulators (posterior-mean optics, D1).
  * Nuisance: radius population (EB grid) + albedo discrepancy SD, chosen on DEVELOPMENT station-days
    only (PROMBIO 2021/2023), then frozen for the final test (PROMBIO 2024).
  * Abundance: replicates averaged in log space. Left-censored counts (< 1000 cells/mL) are integrated
    over log B in [1, 3) with the abundance prior.
  * Development scores are leave-one-station-day-out and are reported as DEVELOPMENT, not validation.

    python phase4/h5_station_biology.py --manifest records/prombio_matching/manifest.csv \
        --outdir phase4/results/h5 [--split development|test]
"""
from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd
from scipy.stats import norm

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "phase2"))

import h1_albedo as H1  # noqa: E402
import seb  # noqa: E402

RAW = os.path.join(HERE, "..", "data", "promice_raw")
MAX_NODE_OFFSET = 1.0


def noon_target(station, date, lat, lon):
    """Albedo over the 3 hours nearest local solar noon (hour-start stamps, midpoints); None + reason."""
    d = pd.read_csv(os.path.join(RAW, f"{station}_hour.csv"), usecols=["time", "dsr_cor", "usr_cor", "snow_height"],
                    parse_dates=["time"])
    day = d[d.time.dt.date == pd.Timestamp(date).date()].copy()
    if day.empty:
        return None, "no station data"
    day["sza"] = seb.solar_zenith_hourly(day.time, lat, lon)
    sel = day.nsmallest(3, "sza")
    ok = sel.dsr_cor.notna() & sel.usr_cor.notna()
    if ok.sum() < 2:
        return None, "fewer than 2 valid noon hours (dsr_cor/usr_cor)"
    if (sel.snow_height.isna() | (sel.snow_height >= 0.02)).any():
        return None, "snow_height >= 0.02 m or missing at noon"
    s = sel[ok]
    return dict(albedo_obs=float(s.usr_cor.sum() / s.dsr_cor.sum()), sza_noon=float(sel.sza.min()),
                hours=",".join(t.strftime("%H") for t in sel.time)), None


def station_days(manifest, split):
    m = pd.read_csv(manifest)
    m = m[(m.role == ("final_test" if split == "test" else "development")) & (m.sampling_stratum == "under radiometer")]
    rows, excl = [], []
    for (st, date), g in m.groupby(["station", "date"]):
        g_ok = g[g.qc_exclusion.isna() | (g.qc_exclusion == "")]
        if g_ok.empty:
            excl.append(dict(station=st, date=date, reason="; ".join(sorted(set(g.qc_exclusion.dropna())))))
            continue
        lat, lon = float(g_ok.match_lat.iloc[0]), float(g_ok.match_lon.iloc[0])
        tg, why = noon_target(st, date, lat, lon)
        if tg is None:
            excl.append(dict(station=st, date=date, reason=why))
            continue
        obs = g_ok[g_ok.obs_state == "observed"].ice_cells_ml.astype(float)
        cens = (g_ok.obs_state == "left_censored").any()
        rows.append(dict(station=st, date=date, n_samples=len(g_ok), n_censored=int((g_ok.obs_state == "left_censored").sum()),
                         log_b=float(np.mean(np.log10(np.clip(obs, 10.0, None)))) if len(obs) else np.nan,
                         censored_only=bool(cens and obs.empty), **tg))
    return pd.DataFrame(rows), pd.DataFrame(excl)


def grids_for(df, opt, ems):
    """Per station-day BBA over the nuisance grid at its abundance (censored-only: prior-weighted mix over [1, 3))."""
    from priors import PriorConfig
    pc = PriorConfig.for_density(690.0)
    out = []
    for r in df.itertuples():
        em = ems[r.node]
        if r.censored_only:
            lbs = em.axes["log_b"][em.axes["log_b"] < 3.0]
            w = norm.pdf(lbs, pc.mu_b, pc.sd_b)
            w = w / w.sum()
            out.append(sum(wi * H1.node_bba(em, lb, opt == "no_algae") for wi, lb in zip(w, lbs)))
        else:
            out.append(H1.node_bba(em, r.log_b, opt == "no_algae"))
    return out


def run(manifest, outdir, split):
    dev, dev_ex = station_days(manifest, "development")
    data = {"development": (dev, dev_ex)}
    if split == "test":
        data["test"] = station_days(manifest, "test")
    nodes = sorted(int(f.split("sza")[1][:-4]) for f in os.listdir(H1.EMU_DIR) if f.startswith("heldout_tddft_D_sza"))
    for k, (df, ex) in data.items():
        if df.empty:
            continue
        df["node"] = [min(nodes, key=lambda n: abs(n - z)) for z in df.sza_noon]
        bad = (df.node - df.sza_noon).abs() > MAX_NODE_OFFSET
        data[k] = (df[~bad].reset_index(drop=True),
                   pd.concat([ex, df[bad].assign(reason="noon SZA outside the emulator domain")[["station", "date", "reason"]]]))
    preds, settings = [], {}
    hyp = [(m, s, d) for m in H1.MU_LNR for s in H1.SD_LNR for d in H1.DISC_SD]
    for opt in H1.OPTICS:
        dv = data["development"][0]
        all_nodes = sorted(set(dv.node) | (set(data["test"][0].node) if "test" in data and not data["test"][0].empty else set()))
        ems = H1.load_emulators(opt, all_nodes)
        Wc = {}

        def W(z, m, s):
            if (z, m, s) not in Wc:
                Wc[(z, m, s)] = H1.prior_weights(ems[z], m, s)
            return Wc[(z, m, s)]
        gd = grids_for(dv, opt, ems)
        ll = np.array([[H1.plot_logp(g, W(z, m, s), y, d) for (m, s, d) in hyp]
                       for g, z, y in zip(gd, dv.node, dv.albedo_obs)])
        # development: leave-one-station-day-out
        for i in range(len(dv)):
            tr = np.arange(len(dv)) != i
            m, s, d = hyp[int(np.argmax(ll[tr].sum(0)))]
            x, cdf, mean, (q05, q50, q95) = H1.predictive(gd[i], W(dv.node[i], m, s), d)
            y = dv.albedo_obs[i]
            preds.append(dict(split="development_LOO", model=opt, station=dv.station[i], date=dv.date[i], albedo_obs=y,
                              pred=mean, err=mean - y, abs_err=abs(mean - y), crps=H1.crps(x, cdf, y),
                              covered90=float(q05 <= y <= q95), width90=q95 - q05))
        m, s, d = hyp[int(np.argmax(ll.sum(0)))]
        settings[opt] = dict(radius_median_um=float(np.exp(m)), radius_ln_sd=float(s), disc_sd=float(d), n_dev=len(dv))
        if "test" in data and not data["test"][0].empty:
            te = data["test"][0]
            gt = grids_for(te, opt, ems)
            for i in range(len(te)):
                x, cdf, mean, (q05, q50, q95) = H1.predictive(gt[i], W(te.node[i], m, s), d)
                y = te.albedo_obs[i]
                preds.append(dict(split="final_test", model=opt, station=te.station[i], date=te.date[i], albedo_obs=y,
                                  pred=mean, err=mean - y, abs_err=abs(mean - y), crps=H1.crps(x, cdf, y),
                                  covered90=float(q05 <= y <= q95), width90=q95 - q05))
        print(f"[H5] {opt} done", flush=True)
    P = pd.DataFrame(preds)
    S = P.groupby(["split", "model"]).agg(n=("station", "size"), mae=("abs_err", "mean"), bias=("err", "mean"),
                                          crps=("crps", "mean"), coverage90=("covered90", "mean"),
                                          width90=("width90", "mean")).reset_index()
    import provenance as PV
    os.makedirs(outdir, exist_ok=True)
    PV.atomic_to_csv(P, os.path.join(outdir, f"h5_predictions_{split}.csv"), index=False, float_format="%.5g")
    PV.atomic_to_csv(S, os.path.join(outdir, f"h5_summary_{split}.csv"), index=False, float_format="%.5g")
    for k, (df, ex) in data.items():
        PV.atomic_to_csv(ex, os.path.join(outdir, f"h5_excluded_{k}.csv"), index=False)
        PV.atomic_to_csv(df.drop(columns=["log_b"]), os.path.join(outdir, f"h5_station_days_{k}.csv"), index=False)
    PV.atomic_write_text(os.path.join(outdir, f"h5_settings_{split}.json"), json.dumps(dict(
        settings=settings, manifest_sha256=PV.file_sha256(manifest), environment=PV.environment()), indent=1))
    return P, S, data


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--manifest", default=os.path.join(HERE, "..", "records", "prombio_matching", "manifest.csv"))
    p.add_argument("--outdir", default=os.path.join(HERE, "results", "h5"))
    p.add_argument("--split", choices=["development", "test"], default="development")
    a = p.parse_args(argv)
    P, S, data = run(a.manifest, a.outdir, a.split)
    pd.set_option("display.width", 200)
    print(S.round(4).to_string(index=False))
    for k, (df, ex) in data.items():
        print(k, "included", len(df), "excluded", len(ex))
        print(ex.to_string(index=False))


if __name__ == "__main__":
    main()
