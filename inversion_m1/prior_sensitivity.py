"""Prior-vs-data contribution for log10 cells, from the stored Baseline-A profiles (no refit).

For a prior N(mu, sd) on log10 cells the posterior-mode profile objective is chi2_profile(g) + ((g-mu)/sd)^2 exactly
(other parameters re-optimised, as in the stored profile). Priors are NOT derived from any count in the evaluation
data: mu is the midpoint of the imposed bounds (3.75) and sd is varied over 0.5/1/2 dex as a SENSITIVITY family.
Reported per setting: contraction = (posterior 95 % width)/(no-prior width), the share of the shrinkage attributable to the
prior = 1 - (posterior width / prior-only width is not used; see 'prior_only_width'), shift of the mode, and the (post hoc,
non-tuning) agreement with counts. A narrower interval is not evidence of better inference.
"""
import json, os
import numpy as np
from scipy.interpolate import PchipInterpolator
import pandas as pd
R = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "records", "inversion_m1", "baseline_A")
fits = {}
for ln in open(os.path.join(R, "per_plot_fits.jsonl")):
    d = json.loads(ln)
    if d.get("status") == "ok": fits[d["sample"]] = d
sc = json.load(open(os.path.join(R, "error_scale_lodo.json")))["per_day_error_scale"]
ct = pd.read_csv(os.path.join(R, "..", "..", "..", "data", "empirical", "cook2020_archive_cell_counts.csv")).set_index("sample")
day = lambda s: "S6:2017-%02d-%02d" % (int(s.split("_")[1]), int(s.split("_")[0]))
fine = np.linspace(1.0, 6.5, 1101)
MU = 3.75
res = {}
for sig_name, sig_of in {"sigma_total (no kappa)": lambda s: sc[day(s)]["sigma_total"], "sigma_obs 0.015": lambda s: 0.015}.items():
    for sd in [None, 2.0, 1.0, 0.5]:
        W, W0, shift, err, cov, n = [], [], [], [], 0, 0
        for s, d in fits.items():
            g = np.array([x for x, v in zip(d["profile_grid"], d["profile_sse"]) if v is not None]); v = np.array([x for x in d["profile_sse"] if x is not None])
            chi = PchipInterpolator(g, (v) / sig_of(s) ** 2)(fine)
            chi0 = chi - chi.min()
            ins0 = fine[chi0 <= 3.84]; w0 = ins0[-1] - ins0[0]
            pen = chi + (0 if sd is None else ((fine - MU) / sd) ** 2)
            i = int(np.argmin(pen)); ins = fine[(pen - pen.min()) <= 3.84]; w = ins[-1] - ins[0]
            W.append(w); W0.append(w0); shift.append(fine[i] - fine[np.argmin(chi)])
            if ct.loc[s, "cells_per_ml"] > 0:
                c = np.log10(ct.loc[s, "cells_per_ml"]); sdc = 1 / np.log(10) / np.sqrt(max(ct.loc[s, "cells_counted"], 1))
                err.append(fine[i] - c); cov += (ins[0] - 2 * sdc <= c <= ins[-1] + 2 * sdc); n += 1
        err = np.array(err)
        res[f"{sig_name} | prior sd={sd}"] = dict(
            median_width_dex=float(np.median(W)), median_width_noprior=float(np.median(W0)),
            median_contraction=float(np.median(np.array(W) / np.array(W0))), prior_95pct_width_dex=None if sd is None else 2 * 1.96 * sd,
            median_abs_mode_shift_dex=float(np.median(np.abs(shift))), bias_vs_counts_dex=float(err.mean()), rmse_vs_counts_dex=float(np.sqrt((err ** 2).mean())),
            counts_within_interval=cov / n)
json.dump(res, open(os.path.join(R, "prior_sensitivity.json"), "w"), indent=1)
for k, v in res.items(): print(k, {a: (round(b, 3) if isinstance(b, float) else b) for a, b in v.items()})
