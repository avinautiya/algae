#!/usr/bin/env python
"""
H3 prerequisites (docs/glacier_model_validation_protocol.md §4), checked BEFORE any algae model is run
for H3. Per station-year, on the bare-ice days where both ablation records and complete forcing exist:

  P1 the two independent ablation records (pressure transducer z_pt_cor, stake/sonic z_stake_cor) agree
     within 25 %: |total_pt / total_stake - 1| <= 0.25;
  P2 the reference SEB (phase4/seb.py, chi = 0, MEASURED daily albedo) reproduces the observed lowering
     within the sensor disagreement: the model total lies between the two sensor totals.

If either fails, H3 is UNTESTED for that station-year. No record is selected.

    python phase4/h3_prerequisites.py --outdir phase4/results/h3_prerequisites
"""

from __future__ import annotations

import argparse
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "phase2"))

import seb  # noqa: E402

AGREE = 0.25


def check(station, year):
    f = seb.load_station_forcing(station, year)
    if f.empty or f.dsr.notna().sum() < 24 * 30:
        return dict(station=station, year=year, status="UNTESTED", reason="insufficient forcing")
    out = {}
    for rec in ("z_pt_cor", "z_stake_cor"):
        if f[rec].notna().sum() < 24 * 10:
            return dict(station=station, year=year, status="UNTESTED", reason=f"{rec} missing")
        r, dd = seb.validate(year, ablation=rec, forcing=f)
        out[rec] = dd
    days = out["z_pt_cor"].index.intersection(out["z_stake_cor"].index)
    if len(days) < 10:
        return dict(station=station, year=year, status="UNTESTED", reason=f"only {len(days)} common bare-ice days")
    pt = float(out["z_pt_cor"].loc[days, "obs"].sum())
    sk = float(out["z_stake_cor"].loc[days, "obs"].sum())
    model = float(out["z_pt_cor"].loc[days, "model"].sum())
    p1 = sk > 0 and abs(pt / sk - 1) <= AGREE
    p2 = min(pt, sk) <= model <= max(pt, sk)
    reason = "" if (p1 and p2) else "; ".join(x for x, ok in (
        (f"P1 sensors disagree (pt/stake {pt / sk:.2f})" if sk > 0 else "P1 stake total <= 0", p1),
        (f"P2 model {model:.3f} outside sensors [{min(pt, sk):.3f}, {max(pt, sk):.3f}] m w.e.", p2)) if not ok)
    return dict(station=station, year=year, n_days=len(days), sw_source=str(f.sw_source.mode().iloc[0]),
                obs_pt_mwe=pt, obs_stake_mwe=sk, model_mwe=model, ratio_model_pt=model / pt if pt else np.nan,
                ratio_model_stake=model / sk if sk else np.nan, P1=bool(p1), P2=bool(p2),
                status="PREREQUISITES MET" if (p1 and p2) else "UNTESTED", reason=reason)


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--outdir", default=os.path.join(HERE, "results", "h3_prerequisites"))
    p.add_argument("--stations", nargs="+", default=["KAN_L", "KAN_M"])
    p.add_argument("--years", nargs="+", type=int, default=list(range(2016, 2024)))
    a = p.parse_args(argv)
    rows = [check(s, y) for s in a.stations for y in a.years]
    df = pd.DataFrame(rows)
    import provenance as PV
    os.makedirs(a.outdir, exist_ok=True)
    PV.atomic_to_csv(df, os.path.join(a.outdir, "h3_prerequisites.csv"), index=False, float_format="%.4g")
    PV.atomic_write_text(os.path.join(a.outdir, "settings.json"), json.dumps(dict(
        agree=AGREE, chi=0.0, z0=1e-3, slab_m=0.1, raw_sha256=open(os.path.join(HERE, "..", "data", "promice_raw", "SHA256SUMS")).read(),
        environment=PV.environment()), indent=1))
    pd.set_option("display.width", 220)
    print(df.round(3).to_string(index=False))


if __name__ == "__main__":
    main()
