#!/usr/bin/env python
"""
Dust prior at the independent site. Mineral dust was measured at S6 (Cook et al. 2020: 342 ug/g in the
algal ice) but not at the southern Greenland site of Chevrollier et al. (2023), who measured no dust and
argue from their spectral fits that dust does not lower the albedo significantly above 350 nm. The maps
and the field validation borrow the S6 dust prior for every site; this script quantifies what that
borrowing does to the retrieved abundance at the southern site.

Variants (same leave-one-out field validation as field_validation.run; only the dust prior changes):
  s6_prior     measured S6 log-normal (median ~3.1e5 ppb, ln-SD 0.46) - the default everywhere
  broad        same median, ln-SD 1.5 (dust essentially free over two orders of magnitude)
  low_dust     median / 10, ln-SD 1.0 (a cleaner site than S6)
  no_dust      dust axis removed

For each it reports the per-site bias / RMSE / coverage and the posterior median dust of the southern
samples (what their own spectra say about dust, conditional on the model).

    python phase4/dust_sensitivity.py --phenol williamson2020 --outdir phase4/results/dust_sensitivity
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

import empirical_data as ED  # noqa: E402
import field_validation as FV  # noqa: E402

MU, SD = ED.dust_prior()
VARIANTS = {
    "s6_prior": dict(dust=True, prior={}),
    "broad": dict(dust=True, prior=dict(sd_lndust=1.5)),
    "low_dust": dict(dust=True, prior=dict(mu_lndust=MU - np.log(10.0), sd_lndust=1.0)),
    "no_dust": dict(dust=False, prior={}),
}


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--phase1-l2", default=os.path.join(HERE, "..", "phase1", "results", "level2"))
    p.add_argument("--phenol", default="tddft", choices=["tddft", "williamson2020"])
    p.add_argument("--biosnicar", default=None)
    p.add_argument("--workers", type=int, default=os.cpu_count() or 1)
    p.add_argument("--outdir", default=os.path.join(HERE, "results", "dust_sensitivity"))
    a = p.parse_args(argv)
    os.makedirs(a.outdir, exist_ok=True)
    rows, dust_post = [], {}
    for v, cfg in VARIANTS.items():
        done = os.path.join(a.outdir, f"variant_{v}.json")
        if os.path.isfile(done):                         # finished before an interruption
            d = json.load(open(done))
            rows += d["rows"]
            if d["dust_post"] is not None:
                dust_post[v] = d["dust_post"]
            print(f"[{v}] loaded", flush=True)
            continue
        n0 = len(rows)
        df, m, res = FV.run(phase1_l2=a.phase1_l2, biosnicar=a.biosnicar, workers=a.workers, cache_dir=a.outdir,
                            phenol=a.phenol, verbose=False, models=("ours",), dust=cfg["dust"],
                            prior_overrides=cfg["prior"])
        df.to_csv(os.path.join(a.outdir, f"samples_{v}.csv"), index=False, float_format="%.5g")
        for r in m[m.method.str.contains("ours")].itertuples():
            rows.append(dict(variant=v, dataset=r.dataset, n=r.n, bias=r.bias, rmse=r.rmse, spearman=r.spearman,
                             coverage95=r.coverage95, coverage95_cal=r.coverage95_cal, tau=res["ours"]["tau"],
                             log_evidence=res["ours"]["total_log_evidence"]))
        if "ours_dust_ppb_q50" in df:
            dust_post[v] = {ds: float(np.nanmedian(df.loc[df.dataset == ds, "ours_dust_ppb_q50"]))
                            for ds in ("s6_2017", "sgris_2021")}
        with open(done, "w") as fh:
            json.dump(dict(rows=rows[n0:], dust_post=dust_post.get(v)), fh, default=float)
        print(f"[{v}] done", flush=True)
    tab = pd.DataFrame(rows)
    tab.to_csv(os.path.join(a.outdir, "dust_sensitivity_metrics.csv"), index=False, float_format="%.4g")
    out = dict(s6_prior_median_ppb=float(np.exp(MU)), s6_prior_ln_sd=float(SD),
               posterior_median_dust_ppb=dust_post, variants={k: v["prior"] for k, v in VARIANTS.items()})
    with open(os.path.join(a.outdir, "dust_sensitivity_summary.json"), "w") as fh:
        json.dump(out, fh, indent=1, default=float)
    pd.set_option("display.width", 200)
    print(tab.round(3).to_string(index=False))
    print(json.dumps(out, indent=1, default=float))


if __name__ == "__main__":
    main()
