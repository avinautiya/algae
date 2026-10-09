#!/usr/bin/env python
"""Automated pre-scoring audit of the H2/H4 population against the frozen protocol (§3, §5, amendment A1).
Re-derives every inclusion rule from the RAW station file and the stored pixel record, independently of
station_pixels.py. Exit code 1 on any violation.

    python phase4/audit_h2h4.py --pixels phase4/results/station_pixels/station_pixels.csv
"""
import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
RAW = os.path.join(HERE, "..", "data", "promice_raw")


def audit(pixels):
    px = pd.read_csv(pixels)
    v = px[px.status == "valid"].copy()
    problems, checks = [], {}
    checks["years_in_A1"] = bool(v.year.between(2016, 2023).all())
    checks["months_JJA"] = bool(pd.to_datetime(v.datetime).dt.month.isin([6, 7, 8]).all())
    checks["one_per_station_day"] = bool(not v.duplicated(["station", "day"]).any())
    scl_ok = v.scl.map(lambda s: all(c == 11 for c in json.loads(s)) and len(json.loads(s)) == 9)
    checks["scl_all_11_3x3"] = bool(scl_ok.all())
    checks["reflectance_finite_positive"] = bool((v[["B2", "B3", "B4", "B8"]] > 0).all().all())
    for st in v.station.unique():
        raw = pd.read_csv(os.path.join(RAW, f"{st}_hour.csv"), usecols=["time", "dsr_cor", "usr_cor", "snow_height"],
                          parse_dates=["time"]).set_index("time")
        for r in v[v.station == st].itertuples():
            t0 = pd.Timestamp(r.datetime).tz_localize(None).floor("h")
            sub = raw.reindex([t0 + pd.Timedelta(hours=k) for k in (-1, 0, 1)])
            if sub[["dsr_cor", "usr_cor"]].isna().any().any():
                problems.append((r.scene_id, "dsr_cor/usr_cor missing"))
            if (sub.snow_height.isna() | (sub.snow_height >= 0.02)).any():
                problems.append((r.scene_id, "snow rule"))
            alb = sub.usr_cor.sum() / sub.dsr_cor.sum()
            if not np.isclose(alb, r.albedo_obs, rtol=1e-5):
                problems.append((r.scene_id, f"target albedo mismatch {alb} vs {r.albedo_obs}"))
            if not np.isclose((sub.dsr_cor - sub.usr_cor).mean(), r.absorbed_obs, rtol=1e-5):
                problems.append((r.scene_id, "absorbed target mismatch"))
            if not np.isclose(sub.dsr_cor.mean(), r.sw_down, rtol=1e-5):
                problems.append((r.scene_id, "forcing SW mismatch"))
    checks["station_rules_and_targets"] = not problems
    # station SW-up / albedo must never be a prediction input: inspect the scorer's input columns
    src = open(os.path.join(HERE, "h2h4_station.py")).read()
    pred_block = src.split("def predict")[1].split("def block_contrast")[0]
    checks["no_station_target_in_inputs"] = ("usr" not in pred_block) and ("albedo_obs" not in pred_block.split("for m in")[0])
    blocks = v[v.year != 2019].groupby(["station", "year"]).size()
    out = dict(n_valid=int(len(v)), n_excluded=int((px.status != "valid").sum()),
               excluded_reasons=px[px.status != "valid"].reason.str[:60].value_counts().to_dict(),
               primary_blocks=int(len(blocks)), per_block={f"{k[0]}_{k[1]}": int(n) for k, n in blocks.items()}, checks=checks, problems=problems[:20],
               passed=all(checks.values()))
    return out


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("--pixels", default=os.path.join(HERE, "results", "station_pixels", "station_pixels.csv"))
    p.add_argument("--out", default=os.path.join(HERE, "..", "records", "h2h4", "audit.json"))
    a = p.parse_args()
    res = audit(a.pixels)
    os.makedirs(os.path.dirname(a.out), exist_ok=True)
    json.dump(res, open(a.out, "w"), indent=1, default=str)
    print(json.dumps(res, indent=1, default=str))
    sys.exit(0 if res["passed"] else 1)
