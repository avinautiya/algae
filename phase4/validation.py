"""
Error metrics against ground truth (synthetic truth or field measurements).

For each method and quantity: bias (mean error), RMSE, MAE, R^2, and - for Bayesian
methods - the empirical coverage of the nominal 95 % credible interval and its mean
width (a calibrated posterior covers ~95 % of truths; much less = overconfident).
"""

from __future__ import annotations

import numpy as np
import pandas as pd


def metrics(truth, est, lo=None, hi=None):
    t = np.asarray(truth, dtype=float).ravel()
    e = np.asarray(est, dtype=float).ravel()
    ok = np.isfinite(t) & np.isfinite(e)
    t, e = t[ok], e[ok]
    err = e - t
    out = dict(n=int(ok.sum()), bias=err.mean(), rmse=np.sqrt(np.mean(err ** 2)), mae=np.abs(err).mean(),
               r2=1 - np.sum(err ** 2) / np.sum((t - t.mean()) ** 2))
    if lo is not None and hi is not None:
        l = np.asarray(lo, dtype=float).ravel()[ok]
        h = np.asarray(hi, dtype=float).ravel()[ok]
        out.update(coverage95=float(np.mean((t >= l) & (t <= h))), ci_width=float(np.mean(h - l)))
    return out


def compare(truth: dict, results: dict) -> pd.DataFrame:
    """truth: quantity -> array; results: method -> quantity -> (est, lo, hi)."""
    rows = []
    for method, qs in results.items():
        for q, triple in qs.items():
            if q not in truth:
                continue
            est, lo, hi = (triple + (None, None))[:3] if isinstance(triple, tuple) else (triple, None, None)
            rows.append(dict(method=method, quantity=q, **metrics(truth[q], est, lo, hi)))
    df = pd.DataFrame(rows)
    return df


def error_reduction(df: pd.DataFrame, ours: str, baseline: str) -> pd.DataFrame:
    """Percent reduction of |bias| and RMSE of `ours` relative to `baseline`, per quantity."""
    a = df[df.method == ours].set_index("quantity")
    b = df[df.method == baseline].set_index("quantity")
    q = a.index.intersection(b.index)
    return pd.DataFrame({"abs_bias_reduction_pct": 100 * (1 - a.loc[q, "bias"].abs() / b.loc[q, "bias"].abs()),
                         "rmse_reduction_pct": 100 * (1 - a.loc[q, "rmse"] / b.loc[q, "rmse"])})


def field_points(csv_path, transform, crs, shape):
    """Read field validation data (columns: lon, lat, cells_per_ml [, any others]) and return
    a DataFrame with pixel row/col indices on the scene grid (points outside are dropped)."""
    from pyproj import Transformer
    df = pd.read_csv(csv_path)
    x, y = Transformer.from_crs(4326, crs, always_xy=True).transform(df.lon.values, df.lat.values)
    col = np.floor((np.asarray(x) - transform.c) / transform.a).astype(int)
    row = np.floor((np.asarray(y) - transform.f) / transform.e).astype(int)
    keep = (row >= 0) & (row < shape[0]) & (col >= 0) & (col < shape[1])
    df = df[keep].copy()
    df["row"], df["col"] = row[keep], col[keep]
    df["log_b_obs"] = np.log10(df.cells_per_ml)
    return df
