#!/usr/bin/env python
"""
Field-site bias study: which physically motivated change removes the abundance over-estimate at S6
without degrading the independent southern-Greenland site? Every variant is scored with the same
leave-one-out field validation (field_validation.run), so hyper-parameters (sigma, ice-radius prior)
are re-selected per variant and nothing is tuned on the scored sample.

Variants (all with the meltwater-unit correction of cell counts, biosnicar_bridge.MELTWATER_TO_BIOSNICAR):
  baseline        algae uniform over the 2 cm sampling layer, no mineral dust
  dust            + mineral-dust axis with the MEASURED S6 loading as prior (342 ug/g, Cook et al. 2020)
  film_uniform    3-layer column (2 mm film + 18 mm crust), algae uniform over both  (control for 'film')
  film            same column, all algae in the top 2 mm (same cells per m^2)
  dust_film       dust + film

Site calibration factor (post hoc, on the best variant): a multiplicative abundance correction
estimated at S6 by leave-one-out (the mean log error of the OTHER S6 samples) is applied (a) at S6, and
(b) unchanged at southern Greenland - the transfer test a satellite application would face.

    python phase4/bias_study.py --phase1-l2 phase1/results/level2 --phenol williamson2020 --outdir ...
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

import emulator as E  # noqa: E402
import field_validation as FV  # noqa: E402

FILM_DZ = 0.002
VARIANTS = {
    "baseline": {},
    "dust": dict(dust_ppb=E.DUST_NODES_PPB),
    "film_uniform": dict(film_dz=FILM_DZ, film_only=False),
    "film": dict(film_dz=FILM_DZ, film_only=True),
    "dust_film": dict(dust_ppb=E.DUST_NODES_PPB, film_dz=FILM_DZ, film_only=True),
}


def site_calibration(df: pd.DataFrame, model="ours"):
    """Leave-one-out multiplicative correction estimated on S6, applied at S6 and transferred."""
    s6 = (df.dataset == "s6_2017").to_numpy() & np.isfinite(df.log_b_obs.to_numpy())
    err = (df[f"{model}_log_b_mean"] - df.log_b_obs).to_numpy()
    shift = np.full(len(df), np.nan)
    for i in range(len(df)):
        others = s6.copy()
        others[i] = False
        shift[i] = np.nanmean(err[others])
    out = {}
    for ds in ("s6_2017", "sgris_2021"):
        m = (df.dataset == ds).to_numpy() & np.isfinite(df.log_b_obs.to_numpy())
        est = df[f"{model}_log_b_mean"].to_numpy()[m] - shift[m]
        lo = df[f"{model}_log_b_q025"].to_numpy()[m] - shift[m]
        hi = df[f"{model}_log_b_q975"].to_numpy()[m] - shift[m]
        t = df.log_b_obs.to_numpy()[m]
        e = est - t
        out[ds] = dict(n=int(m.sum()), bias=float(e.mean()), rmse=float(np.sqrt(np.mean(e ** 2))),
                       coverage95=float(np.mean((t >= lo) & (t <= hi))), mean_shift_dex=float(np.mean(shift[m])))
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--phase1-l2", default=os.path.join(HERE, "..", "phase1", "results", "level2"))
    p.add_argument("--phenol", default="williamson2020", choices=["tddft", "williamson2020"])
    p.add_argument("--tier", default="C", choices=["C", "D"],
                   help="C: calibrated uncomplexed pigment; D: + measured Fe-complex absorption (phenol=tddft only)")
    p.add_argument("--biosnicar", default=None)
    p.add_argument("--variants", nargs="+", default=list(VARIANTS))
    p.add_argument("--workers", type=int, default=os.cpu_count() or 1)
    p.add_argument("--outdir", default=os.path.join(HERE, "results", "bias_study"))
    a = p.parse_args(argv)
    os.makedirs(a.outdir, exist_ok=True)
    rows, extra = [], {}
    for v in a.variants:
        t0 = time.time()
        ov = dict(VARIANTS[v], f_n=(0.0, 1.0, 0.2))
        df, m, res = FV.run(phase1_l2=a.phase1_l2, biosnicar=a.biosnicar, workers=a.workers, cache_dir=a.outdir,
                            phenol=a.phenol, tier=a.tier, verbose=False, emu_overrides=ov, models=("ours",), dust=False)
        df.to_csv(os.path.join(a.outdir, f"samples_{v}.csv"), index=False, float_format="%.5g")
        mo = m[m.method.str.contains("ours")]
        for r in mo.itertuples():
            rows.append(dict(variant=v, dataset=r.dataset, n=r.n, bias=r.bias, rmse=r.rmse, spearman=r.spearman,
                             coverage95=r.coverage95, coverage95_obs=r.coverage95_obs,
                             coverage95_cal=r.coverage95_cal, tau=res["ours"]["tau"],
                             log_evidence=res["ours"]["total_log_evidence"], sigma=res["ours"]["sigma_all"],
                             r_median_um=float(np.exp(res["ours"]["mu_lnr"])), r_lnsd=res["ours"]["sd_lnr"]))
        extra[v] = dict(site_calibration=site_calibration(df), runtime_s=round(time.time() - t0))
        print(f"[{v}] done in {time.time() - t0:.0f} s", flush=True)
    tab = pd.DataFrame(rows)
    tab.to_csv(os.path.join(a.outdir, "bias_study_metrics.csv"), index=False, float_format="%.4g")
    with open(os.path.join(a.outdir, "bias_study_site_calibration.json"), "w") as fh:
        json.dump(extra, fh, indent=1)
    pd.set_option("display.width", 200)
    print(tab.round(3).to_string(index=False))
    print(json.dumps(extra, indent=1))


if __name__ == "__main__":
    main()
