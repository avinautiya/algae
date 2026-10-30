"""Audit of the Baseline-A profile-likelihood intervals for log10 cells (reads stored profile SSEs; no new RT runs).

Construction (baseline_fixed_optics.fit_plot / profile_interval):
  * profile parameter: log10 cells per mL meltwater, grid 1.0..6.5 step 0.5 (the imposed bounds are the grid ends);
  * at each grid value the other three parameters (dust, dark-particle proxy, bubble radius) are re-optimised by
    bounded least squares from ONE warm start (the global optimum) -> a non-global nuisance optimum can only raise the
    SSE and therefore narrow the interval; this is checked separately on a subset (see profile_restart_check);
  * threshold: (SSE(g) - SSE_min) / (sigma*kappa)^2 <= 3.84 (chi2, 1 dof, 95 %);
  * sigma, kappa: leave-one-day-out error scale (TRAINING days only).
A side is 'bound-limited' when the 3.84 threshold is not crossed before the grid end (the data do not close the interval
inside the imposed bounds); 'closed' when it is crossed on that side.
"""
import json, os
import numpy as np

R = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "records", "inversion_m1", "baseline_A")
fits = {}
for ln in open(os.path.join(R, "per_plot_fits.jsonl")):
    d = json.loads(ln)
    if d.get("status") == "ok": fits[d["sample"]] = d
sc = json.load(open(os.path.join(R, "error_scale_lodo.json")))
sc = sc["per_day_error_scale"]

def interval(grid, sse, smin, scale):
    g = np.array([x for x, s in zip(grid, sse) if s is not None], float)
    d = (np.array([s for s in sse if s is not None]) - smin) / scale
    ins = d <= 3.84
    if not ins.any():                       # narrower than the 0.5-dex grid: closed, width reported as the grid step (upper bound)
        return 0.0, 0.5, False, False
    lo_i, hi_i = np.flatnonzero(ins)[[0, -1]]
    lo_open, hi_open = lo_i == 0, hi_i == len(g) - 1
    lo = g[lo_i] if lo_open else float(np.interp(3.84, [d[lo_i], d[lo_i - 1]], [g[lo_i], g[lo_i - 1]]))
    hi = g[hi_i] if hi_open else float(np.interp(3.84, [d[hi_i], d[hi_i + 1]], [g[hi_i], g[hi_i + 1]]))
    return lo, hi, lo_open, hi_open

def day_of(s):
    d, m = (int(x) for x in s.split("_")[:2]); return f"S6:2017-{m:02d}-{d:02d}"

scen = {
  "as_reported (sigma_total*kappa)": lambda s: (sc[day_of(s)]["sigma_total"] * sc[day_of(s)]["kappa"]) ** 2,
  "sigma_total, no kappa": lambda s: sc[day_of(s)]["sigma_total"] ** 2,
  "sigma_obs 0.015, no inflation": lambda s: 0.015 ** 2,
  "Dark Zone sigma 0.003, no inflation": lambda s: 0.003 ** 2,
}
out = {}
for name, f in scen.items():
    w, lo_o, hi_o, both_o, closed = [], 0, 0, 0, 0
    for s, d in fits.items():
        lo, hi, a, b = interval(d["profile_grid"], d["profile_sse"], d["sse"], f(s))
        w.append(hi - lo); lo_o += a; hi_o += b; both_o += (a and b); closed += (not a and not b)
    n = len(fits)
    out[name] = dict(median_width_dex=float(np.median(w)), p10=float(np.percentile(w, 10)), p90=float(np.percentile(w, 90)),
                     frac_lower_end_bound_limited=lo_o / n, frac_upper_end_bound_limited=hi_o / n,
                     frac_both_ends_bound_limited=both_o / n, frac_closed_both_sides=closed / n)
out["_note"] = ("widths are on a 0.5-dex grid with linear interpolation of the threshold crossing; 'bound-limited' = threshold not "
                "reached before the imposed bound (1.0 or 6.5 log10 cells/mL)")
json.dump(out, open(os.path.join(R, "interval_audit.json"), "w"), indent=1)
print(json.dumps(out, indent=1))

# ---- calibration of the scenarios against held-out spectra and against counts (reported, not tuned on) ----
import pandas as pd
ct = pd.read_csv(os.path.join(R, "..", "..", "..", "data", "empirical", "cook2020_archive_cell_counts.csv")).set_index("sample")
cal = {}
for name, f in scen.items():
    cov_spec, cov_cnt, n_cnt, wid = [], 0, 0, []
    for s, d in fits.items():
        sig = np.sqrt(f(s))                       # (sigma*kappa) or sigma
        cov_spec.append(np.mean(np.abs(d["holdout_resid"]) <= 1.96 * sig))
        if ct.loc[s, "cells_per_ml"] > 0:
            lo, hi, *_ = interval(d["profile_grid"], d["profile_sse"], d["sse"], f(s))
            c = np.log10(ct.loc[s, "cells_per_ml"]); sd = 1 / np.log(10) / np.sqrt(max(ct.loc[s, "cells_counted"], 1))
            cov_cnt += (lo - 2 * sd <= c <= hi + 2 * sd); n_cnt += 1
    cal[name] = dict(heldout_window_band_coverage_1p96sigma=float(np.mean(cov_spec)), counts_within_interval=cov_cnt / n_cnt, n_positive_counts=n_cnt)
out["_calibration_vs_heldout_and_counts"] = cal
json.dump(out, open(os.path.join(R, "interval_audit.json"), "w"), indent=1)
print(json.dumps(cal, indent=1))
