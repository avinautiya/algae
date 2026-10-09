"""
Statistics for Phase 3: descriptive summaries, bootstrap intervals, Sobol' indices
(Saltelli 2010 / Jansen estimators via SALib), convergence, and one-at-a-time
(tornado) effects.

Terminology used in outputs (keep it straight in the paper):
  * "95 % interval (P2.5-P97.5)": spread of the forcing itself across the input
    uncertainty - an uncertainty/prediction interval, NOT a confidence interval.
  * "95 % CI of the mean/median": bootstrap confidence interval of that statistic,
    i.e. how precisely the Monte Carlo sample pins it down.
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def replicate_ci(x, replicate, stat=np.mean, level=0.95):
    """CI of a Monte Carlo statistic from R independent replicates: stat over all, +/- t * SE where
    SE = SD(replicate statistics) / sqrt(R). Valid for LHS designs, unlike an iid bootstrap."""
    from scipy.stats import t as tdist
    x, rep = np.asarray(x, float), np.asarray(replicate)
    ok = np.isfinite(x)
    vals = np.array([stat(x[ok & (rep == r)]) for r in np.unique(rep[ok])])
    R = vals.size
    if R < 2:
        return np.nan, np.nan
    se = vals.std(ddof=1) / np.sqrt(R)
    h = tdist.ppf(0.5 + level / 2, R - 1) * se
    m = stat(x[ok])
    return m - h, m + h


def describe(x, n_boot: int = 2000, seed: int = 0, replicate=None) -> dict:
    """Summary of a Monte Carlo output. With `replicate` (independent LHS replicate labels) the CIs of the
    mean and median come from the replicate spread; otherwise from an iid bootstrap (valid only for
    simple random samples)."""
    x_all = np.asarray(x, dtype=float)
    if replicate is not None:
        rep = np.asarray(replicate)
        out = describe(x_all, n_boot, seed)
        out["mean_ci_lo"], out["mean_ci_hi"] = replicate_ci(x_all, rep, np.mean)
        out["median_ci_lo"], out["median_ci_hi"] = replicate_ci(x_all, rep, np.median)
        out["ci_method"] = f"{len(np.unique(rep))} LHS replicates"
        return out
    x = x_all[np.isfinite(x_all)]
    rng = np.random.default_rng(seed)
    idx = rng.integers(0, x.size, size=(n_boot, x.size))
    boot = x[idx]
    bm, bmed = boot.mean(axis=1), np.median(boot, axis=1)
    q = np.percentile(x, [2.5, 5, 25, 50, 75, 95, 97.5])
    return dict(n=x.size, mean=x.mean(), sd=x.std(ddof=1), median=q[3],
                p2_5=q[0], p5=q[1], q25=q[2], q75=q[4], p95=q[5], p97_5=q[6], iqr=q[4] - q[2],
                mean_ci_lo=np.percentile(bm, 2.5), mean_ci_hi=np.percentile(bm, 97.5),
                median_ci_lo=np.percentile(bmed, 2.5), median_ci_hi=np.percentile(bmed, 97.5),
                cv=x.std(ddof=1) / abs(x.mean()) if x.mean() != 0 else np.nan)


def describe_frame(df: pd.DataFrame, columns, **kw) -> pd.DataFrame:
    return pd.DataFrame({c: describe(df[c].to_numpy(), **kw) for c in columns}).T


# --------------------------------------------------------------------------- #
# Sobol'                                                                       #
# --------------------------------------------------------------------------- #
def sobol(problem: dict, Y, second_order: bool = True, n_resamples: int = 1000, seed: int = 0):
    """SALib Sobol' analysis -> (first/total table, second-order matrix or None).

    S1  first-order index: share of Var(Y) explained by the parameter alone.
    ST  total-effect index: share including all its interactions; ST - S1 = interaction share.
    Confidence intervals: bootstrap (n_resamples), 95 %.
    """
    from SALib.analyze import sobol as sobol_analyze

    Y = np.asarray(Y, dtype=float)
    res = sobol_analyze.analyze(problem, Y, calc_second_order=second_order,
                                num_resamples=n_resamples, conf_level=0.95, seed=seed)
    names = sorted(set(problem["groups"]), key=problem["groups"].index) if "groups" in problem \
        else problem["names"]
    tab = pd.DataFrame({"S1": res["S1"], "S1_conf": res["S1_conf"],
                        "ST": res["ST"], "ST_conf": res["ST_conf"]}, index=names)
    tab["interaction"] = tab.ST - tab.S1
    s2 = None
    if second_order:
        s2 = pd.DataFrame(res["S2"], index=names, columns=names)
    return tab, s2


def sobol_convergence(problem: dict, Y, n_base: int, second_order: bool = True,
                      steps=(64, 128, 256, 512, 1024, 2048, 4096), n_resamples: int = 200):
    """Indices recomputed on the first n Saltelli blocks (valid because the Sobol'
    sequence prefix of length 2^k is itself balanced). Returns a long DataFrame."""
    D = len(set(problem["groups"])) if "groups" in problem else problem["num_vars"]
    block = 2 * D + 2 if second_order else D + 2
    rows = []
    for n in steps:
        if n > n_base:
            break
        tab, _ = sobol(problem, np.asarray(Y)[: n * block], second_order, n_resamples)
        for name, r in tab.iterrows():
            rows.append(dict(n_base=n, parameter=name, S1=r.S1, ST=r.ST, ST_conf=r.ST_conf))
    return pd.DataFrame(rows)


# --------------------------------------------------------------------------- #
# One-at-a-time (tornado)                                                      #
# --------------------------------------------------------------------------- #
def tornado_design(space, lo_q: float = 0.05, hi_q: float = 0.95):
    """For each parameter: all others at their median, this one at its lo/hi quantile."""
    base = space.medians()
    lo, hi = space.quantiles(lo_q), space.quantiles(hi_q)
    rows = [dict(base, _param="__base__", _level="base")]
    for n in space.names:
        rows.append(dict(base, **{n: lo[n]}, _param=n, _level="lo"))
        rows.append(dict(base, **{n: hi[n]}, _param=n, _level="hi"))
    return rows


def tornado_table(df: pd.DataFrame, output: str) -> pd.DataFrame:
    base = df.loc[df._param == "__base__", output].iloc[0]
    out = []
    for n, g in df[df._param != "__base__"].groupby("_param"):
        lo = g.loc[g._level == "lo", output].iloc[0]
        hi = g.loc[g._level == "hi", output].iloc[0]
        out.append(dict(parameter=n, base=base, at_low=lo, at_high=hi,
                        d_low=lo - base, d_high=hi - base, swing=abs(hi - lo)))
    return pd.DataFrame(out).sort_values("swing", ascending=False).reset_index(drop=True)
