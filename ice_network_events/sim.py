"""Run a method on the reduced model (RTM_seq) and compute event metrics. Wet threshold as in x09_diagnostics: 0.5 % of the mean inflow (chatter-free)."""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import methods as M, series as S, rtm_seq as RS, forcing as F
from scipy.interpolate import PchipInterpolator


def fns_tr(ser):
    t = ser["t"]; Qf = PchipInterpolator(t, ser["Q"]); Hf = PchipInterpolator(t, ser["H"])
    return (lambda x: np.atleast_1d(Qf(x)) * 1.0 if np.ndim(x) else float(Qf(x))), (lambda x: (Hf(x) / np.maximum(Qf(x), 1e-12)))


def event_metrics(t, q1, q2, qbar):
    wet = q2 > 0.005 * qbar; dt = np.gradient(t); d1 = float(np.trapezoid(q1, t)); d2 = float(np.trapezoid(q2, t)); tr = int(np.sum(np.abs(np.diff(wet.astype(int)))))
    wh = float(np.sum(dt[wet]) / 3600); tot = t[-1] / 3600
    cls = "dry" if wh < 0.05 * tot else "wet" if wh > 0.95 * tot else "intermittent"
    return dict(wet_h=wh, cls=cls, trans=tr, cum1=d1, cum2=d2, share_deliv=d2 / max(d1 + d2, 1e-12), share_end=float(q2[-1] / max(q1[-1] + q2[-1], 1e-12)))


def fixed_tabs(p, stations=17):
    import rtm
    q1 = dict(rtm.full(p), T=1.0); return (rtm.tabulate(q1, 1, stations, "fixed"), rtm.tabulate(q1, 2, stations, "fixed"))


def se_run(p, ser, ftabs):
    """water-storage-only model (fixed geometry, no heat): the cheap correction B; returns the time series too"""
    m = S.means(ser); dur = m["duration"]; r = RS.run(p, [(dur, 0, 0)], "fixed", tabs=ftabs, t_end=dur, lake_T0=ser["H"][0] / ser["Q"][0], fns=fns_tr(ser))
    return r, event_metrics(r["t"], r["Q1"], r["Q2"], m["Q_bar"])


def rtm_run(p, ser, method, th, ftabs=None, se=None):
    m = S.means(ser); dur = m["duration"]
    if method == "SE": return se_run(p, ser, ftabs)[1], None
    segs = M.forcing_for(method, ser, th, dur, se)
    fns = fns_tr(ser) if segs is None else F.make(segs, 900.0)
    r = RS.run(p, [(dur, 0, 0)], "thermal", tabs=th["tabs"], t_end=dur, lake_T0=ser["H"][0] / ser["Q"][0], fns=fns)
    return event_metrics(r["t"], r["Q1"], r["Q2"], m["Q_bar"]), segs


def all_methods(p, ser, th=None, ftabs=None):
    th = th or M.thresholds(p); ftabs = ftabs or fixed_tabs(p); r, _ = se_run(p, ser, ftabs); m = S.means(ser)
    se = dict(t=r["t"], wet=r["Q2"] > 0.005 * m["Q_bar"]); out = {}
    for meth in ("TR", "SE", "AVG", "PW6", "PW3", "EA", "PWn"):
        out[meth], segs = rtm_run(p, ser, meth, th, ftabs, se)
        if segs is not None: out[meth]["n_blocks"] = len(segs)
    return out, th, ftabs
