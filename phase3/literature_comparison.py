#!/usr/bin/env python
"""
Headline forcing and melt next to published estimates (Cook et al. 2020; Williamson et al. 2020).

Published numbers are DAILY values for abundance classes at S6, on specific days, so the comparison is
made on the same footing:

  1. abundance: each published class (mean cells mL^-1; Monte Carlo over its reported spread, log-normal
     with the class mean and SD) - all other inputs drawn from the Phase 3 parameter PDFs (molecular,
     cellular, ice, dust);
  2. albedo reduction d_alpha(SZA) from the Phase 3 forward model at six solar zenith angles (the
     spectral weighting is the clear-sky BioSNICAR illumination at each angle);
  3. hourly forcing  RF_h = d_alpha(SZA_h) x SW_h with SW_h the MEASURED hourly downwelling shortwave
     (PROMICE KAN_M, tilt-corrected 'dsr_cor') on the published day, SZA_h at S6 (67.08 N, 49.38 W);
     daily mean RF = mean over 24 h;
  4. melt potential  M = RF_daily x 86400 s / (334 J g^-1 x 1 g cm^-3 x 10^4 cm^2 m^-2)  [cm w.e. d^-1],
     i.e. all absorbed energy goes to melt (the published methods assume the same, minus 1-5 % to
     photochemistry - within the uncertainty).

Published values (quotes in data/empirical/SOURCES.md, 'Literature comparison'):
  Cook et al. 2020 (TC 14:309), S6, 21 Jul 2017
    Hbio 2.9e4 +/- 2.01e4 cells/mL (1 SD): daily-mean RF 116 W m^-2; melt 1.35 +/- 0.01 (SE) cm w.e. (RF
         method), 1.37 +/- 0.48 (SE) (energy balance); 26.15 +/- 3.77 % of local melt
    Lbio 4.73e3 +/- 2.57e3: RF 65 W m^-2; melt 1.01 +/- 0.01, 0.95 +/- 0.41 cm w.e.; 21.62 +/- 5.07 %
  Williamson et al. 2020 (PNAS 117:5694), S6 area, 26 Jul 2016
    high 8,989 +/- 4,773 cells/mL (n = 103): 1.86 +/- 0.99 (SE) cm w.e. d^-1
    low 186 +/- 276 (n = 27): 0.03 +/- 0.00 cm w.e. d^-1
Note: 116 W m^-2 as a 24 h mean would melt 3.0 cm w.e. d^-1 with the same latent heat, not 1.35; the
published RF and melt of Cook et al. are therefore not mutually consistent under a 24 h-mean reading
(their RF is probably a daytime mean). Both are listed; the melt is the like-for-like comparison.

    python phase3/literature_comparison.py --phase1-l2 phase1/results/level2 --n 200 --outdir ...
"""

from __future__ import annotations

import argparse
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

LAT, LON = 67.08, -49.38                       # S6
L_F = 334.0 * 1e4                              # J per (cm w.e. over 1 m^2)
SZA_NODES = np.array([40.0, 50.0, 60.0, 70.0, 80.0, 88.0])
TIERS = ("A", "C", "D")

CASES = [
    dict(case="Cook2020_Hbio", source="Cook et al. 2020", date="2017-07-21", conc=2.9e4, conc_sd=2.01e4,
         rf_pub=116.0, melt_pub=1.35, melt_pub_se=0.01, melt_eb=1.37, melt_eb_se=0.48),
    dict(case="Cook2020_Lbio", source="Cook et al. 2020", date="2017-07-21", conc=4.73e3, conc_sd=2.57e3,
         rf_pub=65.0, melt_pub=1.01, melt_pub_se=0.01, melt_eb=0.95, melt_eb_se=0.41),
    dict(case="Williamson2020_high", source="Williamson et al. 2020", date="2016-07-26", conc=8989.0,
         conc_sd=4773.0, rf_pub=np.nan, melt_pub=1.86, melt_pub_se=0.99, melt_eb=np.nan, melt_eb_se=np.nan),
    dict(case="Williamson2020_low", source="Williamson et al. 2020", date="2016-07-26", conc=186.0,
         conc_sd=276.0, rf_pub=np.nan, melt_pub=0.03, melt_pub_se=0.0, melt_eb=np.nan, melt_eb_se=np.nan),
]


def hourly_sw(date):
    """Measured hourly SW down (PROMICE KAN_M, W m^-2) and solar zenith at S6 for one UTC day."""
    d = pd.read_csv(ED.path("promice_KAN_M_hour_JJA_radiation.csv"), parse_dates=["time"])
    d = d[d.time.dt.strftime("%Y-%m-%d") == date].set_index("time").sort_index()
    if len(d) < 20:
        raise RuntimeError(f"PROMICE KAN_M has only {len(d)} hours on {date}")
    # fill isolated missing hours by linear interpolation in time
    full = pd.date_range(pd.Timestamp(date), periods=24, freq="h")
    sw = d["dsr_cor"].reindex(full).interpolate(limit_direction="both").to_numpy(float)
    t = pd.Series(full + pd.Timedelta(minutes=30))           # hour-centred
    mu, _ = ED._solar_mu(t, LAT, LON)
    sza = np.degrees(np.arccos(np.clip(np.asarray(mu, float), -1, 1)))
    return np.clip(sw, 0, None), sza


def lognormal_conc(mean, sd, u):
    s2 = np.log(1.0 + (sd / mean) ** 2)
    from scipy.stats import norm
    return np.exp(np.log(mean) - 0.5 * s2 + np.sqrt(s2) * norm.ppf(u))


_FM = None


def _init(kw):
    global _FM
    from forward_model import ForwardModel
    _FM = ForwardModel(**kw)


def _dalpha(p):
    """Albedo reduction weighted by the normalised illumination spectrum (= rf for SW_down = 1; the model
    is built with sw_down=1.0) at every SZA node, per tier: array (tier, sza)."""
    out = np.full((len(TIERS), len(SZA_NODES)), np.nan)
    for j, z in enumerate(SZA_NODES):
        _FM.sza = float(z)
        r = _FM.evaluate(p)
        for i, t in enumerate(TIERS):
            out[i, j] = r[f"rf_{t}"]
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--phase1-l2", default=os.path.join(HERE, "..", "phase1", "results", "level2"))
    p.add_argument("--functional", default="B3LYP")
    p.add_argument("--demo", action="store_true")
    p.add_argument("--biosnicar", default=None)
    p.add_argument("--n", type=int, default=200, help="Monte Carlo samples per case")
    p.add_argument("--workers", type=int, default=os.cpu_count() or 1)
    p.add_argument("--seed", type=int, default=2024)
    p.add_argument("--outdir", default=os.path.join(HERE, "results", "literature_comparison"))
    a = p.parse_args(argv)
    t0 = time.time()
    os.makedirs(a.outdir, exist_ok=True)

    import biosnicar_bridge as bb
    import cell_optics as co
    import tddft_calibration as TC
    from param_space import ParameterSpace, default_parameters
    lig = co.demo_spectrum(bb.locate_biosnicar(a.biosnicar)) if a.demo else co.load_phase1(a.phase1_l2, "level2",
                                                                                          a.functional)
    cal = TC.calibrate(lig, verbose=False)
    space = ParameterSpace(default_parameters(include_tier_d=True, calibration=cal.summary()))
    kw = dict(phase1_l2=None if a.demo else a.phase1_l2, functional=a.functional, demo=a.demo,
              biosnicar=a.biosnicar, sw_down=1.0, qtable_cache=os.path.join(a.outdir, "qstar_table.npz"),
              calibration_point=cal.point())
    _init(kw)                                 # build once (also fills the Q* table cache) before forking

    rows, draws = [], []
    rng = np.random.default_rng(a.seed)
    for c in CASES:
        sw, sza = hourly_sw(c["date"])
        X = space.as_dicts(space.lhs(a.n, seed=int(rng.integers(1 << 30))))
        conc = lognormal_conc(c["conc"], c["conc_sd"], (np.arange(a.n) + rng.random(a.n)) / a.n)
        rng.shuffle(conc)
        for x, cc in zip(X, conc):
            x["conc_cells_ml"] = float(cc)
        if a.workers > 1:
            import multiprocessing as mp
            with mp.get_context("fork").Pool(a.workers) as pool:
                da = np.array(pool.map(_dalpha, X, chunksize=4))
        else:
            da = np.array([_dalpha(x) for x in X])                # (n, tier, sza)
        zh = np.clip(sza, SZA_NODES[0], SZA_NODES[-1])
        for i, t in enumerate(TIERS):
            dah = np.array([np.interp(zh, SZA_NODES, d[i]) for d in da])        # (n, 24)
            dah[:, sza > 90] = 0.0
            rf = (dah * sw[None, :]).mean(axis=1)                               # daily mean W m^-2
            melt = rf * 86400.0 / L_F
            for k in range(a.n):
                draws.append(dict(case=c["case"], tier=t, conc=conc[k], rf_daily=rf[k], melt=melt[k]))
            rows.append(dict(case=c["case"], source=c["source"], date=c["date"], tier=t, conc_mean=c["conc"],
                             rf_daily_median=np.median(rf), rf_daily_p2_5=np.percentile(rf, 2.5),
                             rf_daily_p97_5=np.percentile(rf, 97.5), melt_median=np.median(melt),
                             melt_mean=melt.mean(), melt_p2_5=np.percentile(melt, 2.5),
                             melt_p97_5=np.percentile(melt, 97.5), rf_pub=c["rf_pub"], melt_pub=c["melt_pub"],
                             melt_pub_se=c["melt_pub_se"], melt_eb_pub=c["melt_eb"], melt_eb_pub_se=c["melt_eb_se"],
                             sw_daily_mean=float(sw.mean())))
        print(f"[{c['case']}] done ({time.time() - t0:.0f} s)", flush=True)
    tab = pd.DataFrame(rows)
    tab.to_csv(os.path.join(a.outdir, "literature_comparison.csv"), index=False, float_format="%.4g")
    pd.DataFrame(draws).to_csv(os.path.join(a.outdir, "literature_comparison_draws.csv"), index=False,
                               float_format="%.5g")
    pd.set_option("display.width", 220)
    print(tab.drop(columns=["source", "date"]).round(3).to_string(index=False))
    with open(os.path.join(a.outdir, "literature_comparison_meta.json"), "w") as fh:
        json.dump(dict(n=a.n, sza_nodes=SZA_NODES.tolist(), demo=a.demo, functional=a.functional,
                       runtime_s=round(time.time() - t0)), fh, indent=1)


if __name__ == "__main__":
    main()
