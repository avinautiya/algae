#!/usr/bin/env python
"""Host-coupling demonstration on the named host: the reference point SEB (phase4/seb.py), KAN_L 2022.

KAN_L 2022 is the one station-year whose measured-albedo SEB lies within the ablation-sensor range. That
verdict depends on the low-sun SW rule (P6-SW-1), and H3 remains UNTESTED.

All runs use IDENTICAL hourly forcing (seb.load_station_forcing, 2022-07-24..2022-08-26):

  0. host as is: measured station albedo (daily mean of hourly usr_cor/dsr_cor; days without albedo
     get missing SW, never filled);
  Mode A (forward, prescribed abundance): host algae-free albedo alpha_0(t) and abundance B(t) from the
     satellite retrieval of the 8 valid S2 station-days (TD-DFT D posterior means of r, dust, log B),
     linearly interpolated in time. Established optics (Tier A) and the proposed optics (TD-DFT D) are
     evaluated at the IDENTICAL state (r, dust, B, SZA node); host albedo = alpha_0 - dalpha;
  Mode B (state estimation): observed station albedo used unchanged; counterfactual = observed + dalpha
     (signed). NOT an independent validation.

Outputs: modelled increments (m w.e.) and potential melt, labelled as modelled. Also runtime, energy
closure and n_sub convergence (numerical error). The optics come from frozen held-out emulators by
multilinear lookup (unqualified surrogate use, flagged), with posterior-mean calibration (single draw
'posterior_mean', D1).

    python common/run_budgeted.py --short --mem-mb 400 --name coupling_demo -- \
        python3 -m model_integration.example_reference_seb
"""
from __future__ import annotations

import json
import os
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "phase4"), str(ROOT / "phase2")]

from model_integration import host_modes as HM  # noqa: E402

STATION, YEAR = "KAN_L", 2022
OUT = ROOT / "records" / "host_coupling"


def retrieved_states():
    import scene_product as SP
    import emulator as E
    px = pd.read_csv(ROOT / "records" / "station_pixels" / "station_pixels.csv")
    v = px[(px.status == "valid") & (px.station == STATION) & (px.year == YEAR)].reset_index(drop=True)
    nodes = sorted(int(f.split("sza")[1][:-4]) for f in os.listdir(SP.EMU_DIR) if f.startswith("heldout_tddft_D_sza"))
    v["node"] = [min(nodes, key=lambda n: abs(n - z)) for z in v.sza]
    v = v[(v.node - v.sza).abs() <= 1.0].reset_index(drop=True)
    cal = json.load(open(SP.SETTINGS))["settings"][f"{SP.CAL_FOLD}/tddft_D"]
    rows = []
    for n in sorted(v.node.unique()):
        idx = np.flatnonzero(v.node == n)
        em = SP.add_indicator(E.Emulator.load(os.path.join(SP.EMU_DIR, f"heldout_tddft_D_sza{n}.npz")).refine_log_b(0.05))
        r, pc = SP.posterior(em, v.loc[idx, ["B2", "B3", "B4", "B8"]].to_numpy(float), cal)
        for k, i in enumerate(idx):
            rows.append(dict(day=v.day[i], node=n, log_b=r["log_b_mean"][k], r_um=r["r_um_mean"][k],
                             dust_ppb=r["dust_ppb_mean"][k], f_n=float(pc.f_alpha / (pc.f_alpha + pc.f_beta))))
    return pd.DataFrame(rows).sort_values("day").reset_index(drop=True)


def optics_at_state(st):
    """alpha_0 (algae-free) and dalpha for Tier A and TD-DFT D at the identical state."""
    import scene_product as SP
    import emulator as E
    out = []
    cache = {}
    for r in st.itertuples():
        row = dict(day=r.day)
        for opt in ("tierA_empirical", "tddft_D"):
            key = (opt, r.node)
            if key not in cache:
                em = E.Emulator.load(os.path.join(SP.EMU_DIR, f"heldout_{opt}_sza{r.node}.npz"))
                cache[key] = (em, em.interpolator("bba"), em.interpolator("rf_algae"), float(em.meta["sw_down"]))
            em, Ib, Ir, sw = cache[key]
            coords = {"log_b": r.log_b, "f_n": r.f_n, "r_um": r.r_um, "dust_ppb": np.log10(r.dust_ppb + 100.0)}
            x = np.array([[coords[n] for n in em.active]])
            bba, rf = float(Ib(x)[0]), float(Ir(x)[0])
            row[f"dalpha_{opt}"] = rf / sw
            row[f"alpha0_{opt}"] = bba + rf / sw
        out.append(row)
    d = pd.DataFrame(out)
    d["alpha0"] = d.alpha0_tddft_D
    d["alpha0_consistency"] = (d.alpha0_tddft_D - d.alpha0_tierA_empirical).abs()
    return d


def hourly(series_days, values, times):
    x = (pd.to_datetime(series_days) - times.iloc[0]).dt.total_seconds().to_numpy()
    t = (times - times.iloc[0]).dt.total_seconds().to_numpy()
    return np.interp(t, x + 43200.0, values)


def main():
    import seb
    t0 = time.time()
    st = retrieved_states()
    op = optics_at_state(st)
    f = seb.load_station_forcing(STATION, YEAR, months=(7, 8))
    d0, d1 = pd.Timestamp(op.day.iloc[0]), pd.Timestamp(op.day.iloc[-1]) + pd.Timedelta(days=1)
    f = f[(f.time >= d0) & (f.time < d1)].reset_index(drop=True)
    alpha0 = hourly(op.day, op.alpha0.to_numpy(), f.time)
    da = {o: hourly(op.day, op[f"dalpha_{o}"].to_numpy(), f.time) for o in ("tierA_empirical", "tddft_D")}
    host = HM.HostContract("reference point SEB phase4/seb.py", "absent", "prognostic",
                           note="algae-free host albedo from the satellite-retrieved ice/dust state")
    obs = HM.HostContract("reference SEB with measured station albedo", "implicit_observed", "observed")
    day = f.time.dt.floor("D")
    alb_meas = f.albedo.groupby(day).transform("mean").to_numpy()
    res = {}
    for o in da:
        res[f"modeA_{o}"] = HM.reference_seb_paired(host, f, alpha0 - da[o], alpha0, mode="A", method="anomaly")
    fB = f.copy()
    fB.loc[np.isnan(alb_meas), "sw_down"] = np.nan                  # days without measured albedo: SW missing
    ab = np.where(np.isnan(alb_meas), 0.5, alb_meas)                 # fill never used (SW missing there)
    for o in da:
        res[f"modeB_{o}"] = HM.reference_seb_paired(obs, fB, ab, np.clip(ab + da[o], 0, 1), mode="B")
    # numerical checks on the Mode A TD-DFT run: energy closure and timestep convergence
    m1 = seb.run(f, alpha0 - da["tddft_D"], n_sub=1, return_budget=True)
    m4 = seb.run(f, alpha0 - da["tddft_D"], n_sub=4)
    closure = float(np.nanmax(np.abs(m1.attrs.get("budget", {}).get("closure_residual", [0.0])))) if hasattr(m1, "attrs") else np.nan
    num = dict(melt_n_sub1_mwe=float(m1.melt_mwe[m1.valid].sum()), melt_n_sub4_mwe=float(m4.melt_mwe[m4.valid].sum()),
               max_abs_alpha0_tierA_vs_tddftD=float(op.alpha0_consistency.max()), closure_residual=closure,
               runtime_s=round(time.time() - t0, 1))
    OUT.mkdir(parents=True, exist_ok=True)
    st.merge(op, on="day").to_csv(OUT / "states_and_optics.csv", index=False, float_format="%.5g")
    tab = pd.DataFrame([dict(run=k, **{kk: (vv if not isinstance(vv, dict) else json.dumps(vv)) for kk, vv in v.items()})
                        for k, v in res.items()])
    tab.to_csv(OUT / "coupling_runs.csv", index=False, float_format="%.5g")
    json.dump(dict(numerics=num, window=[str(d0.date()), str((d1 - pd.Timedelta(days=1)).date())], station=STATION,
                   calibration_draw_ids=["posterior_mean"], optics_lookup="frozen held-out emulators, multilinear; "
                   "unqualified surrogate use (benchmark qualification pending)", h3_status="UNTESTED",
                   quantities="modelled increments and potential melt; not observed ablation; runoff not modelled"),
              open(OUT / "numerics.json", "w"), indent=1)
    pd.set_option("display.width", 220)
    print(st.merge(op, on="day").round(4).to_string(index=False))
    print(tab[["run", "modelled_increment_mwe", "potential_melt_mwe", "valid_hours", "missing_hours", "validation_status"]].to_string(index=False))
    print(json.dumps(num, indent=1))


if __name__ == "__main__":
    main()
