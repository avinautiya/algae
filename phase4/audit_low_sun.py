#!/usr/bin/env python
"""Audit of the low-sun SW replacement rule in seb.load_station_forcing (missing SW -> 0 when sun elevation < 5 deg).

Per station-year (JJA):
  * replaced hours by elevation bin (< 0, 0-2, 2-5 deg) and the longest run;
  * an UPPER BOUND on the shortwave energy those hours could have carried: clear-sky SW at the hour's
    elevation (seb's clear-sky parameterisation), compared with the measured SW energy of the period;
  * H3 prerequisite verdicts under the strict rule (replace only when the sun is below the horizon).
Replaced values are modelled replacements, never measurements.
"""
import json
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, "..", "phase2"))
import seb  # noqa: E402
import biosnicar_bridge as bb  # noqa: E402
import h3_prerequisites as H3  # noqa: E402


def audit(station, year):
    f = seb.load_station_forcing(station, year, months=(6, 7, 8))
    rep = f.sw_source == "low_sun_zero"
    sza = seb.solar_zenith_hourly(f.time, *seb.STATION_LATLON[station])
    elev = 90 - sza
    ub = np.array([bb.sw_down_clear_sky(z) if z < 90 else 0.0 for z in sza])
    runs = (rep != rep.shift()).cumsum()[rep]
    meas = f.sw_down.where(f.sw_source.isin(["dsr_cor", "dsr"]))
    return dict(station=station, year=year, hours=len(f), replaced=int(rep.sum()),
                replaced_elev_lt0=int((rep & (elev < 0)).sum()), replaced_elev_0_2=int((rep & (elev >= 0) & (elev < 2)).sum()),
                replaced_elev_2_5=int((rep & (elev >= 2)).sum()),
                longest_run_h=int(runs.value_counts().max()) if rep.any() else 0,
                upper_bound_energy_MJ_m2=float(ub[rep.to_numpy()].sum() * 3600 / 1e6),
                measured_energy_MJ_m2=float(np.nansum(meas) * 3600 / 1e6),
                still_missing_daylight_h=int((f.sw_source == "missing").sum()))


if __name__ == "__main__":
    rows = [audit(s, y) for s in ("KAN_L", "KAN_M") for y in range(2016, 2024)]
    A = pd.DataFrame(rows)
    A["upper_bound_fraction"] = A.upper_bound_energy_MJ_m2 / A.measured_energy_MJ_m2
    base = pd.DataFrame([H3.check(s, y) for s, y in [("KAN_L", 2016), ("KAN_L", 2022), ("KAN_L", 2019), ("KAN_L", 2020)]])
    seb.LOW_SUN_ELEV_DEG = 0.0
    strict = pd.DataFrame([H3.check(s, y) for s, y in [("KAN_L", 2016), ("KAN_L", 2022), ("KAN_L", 2019), ("KAN_L", 2020)]])
    out = os.path.join(HERE, "..", "records", "low_sun_audit")
    os.makedirs(out, exist_ok=True)
    A.to_csv(os.path.join(out, "replacements.csv"), index=False, float_format="%.4g")
    cmp_ = base[["station", "year", "status", "n_days", "model_mwe", "obs_pt_mwe", "obs_stake_mwe"]].merge(
        strict[["station", "year", "status", "n_days", "model_mwe"]], on=["station", "year"], suffixes=("_5deg", "_strict0deg"))
    cmp_.to_csv(os.path.join(out, "h3_prerequisites_rule_sensitivity.csv"), index=False, float_format="%.4g")
    pd.set_option("display.width", 220)
    print(A.round(3).to_string(index=False))
    print(cmp_.round(3).to_string(index=False))
