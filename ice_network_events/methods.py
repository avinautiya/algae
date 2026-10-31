"""Threshold-aware averaging for the two-outlet lake.
Notation: outlet 1 (low sill s1) is always the first to overtop; outlet 2 (sill s2 = s1 + gap) is "wet" when stage > s2.
  I*  = Q1(stage = s2)            steady inflow at which outlet 2 just activates (exact for fixed geometry; drifts with sill lowering, see event_theory.md)
  g1  = dQ1/dH at s2, tau_c = A/g1  stage relaxation time at the threshold (linearised lake: A dH/dt = I - Q1(H)).
  I_f = exponentially filtered inflow, I_f' = (I - I_f)/tau_c.
Criterion A (a priori, from the inflow series, I*, tau_c): U = min(f_above, 1 - f_above), f_above = fraction of time I_f > I*.  Ordinary averaging predicts wet iff mean(I) > I* (f_above in {0,1}); U >> 0 means the mean misplaces the threshold crossing.
Methods (all give a forcing for the same reduced model rtm_seq / for the network):
  AVG  constant mean forcing;  PW6  6 h block means;  PWn  equal blocks, same number as EA;  EA  blocks bounded by the crossings of I_f with I*, runs shorter than tau_c merged into the neighbour;  TR  the series itself."""
import os, sys
import numpy as np
from scipy.interpolate import PchipInterpolator
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "ice_network_capture")); sys.path.insert(0, os.path.join(HERE, "..", "ice_network_sequence"))
import rtm, series as S


def thresholds(p, tabs=None, stations=17):
    p = rtm.full(p)
    if tabs is None: q1 = dict(p, T=1.0); tabs = (rtm.tabulate(q1, 1, stations, "thermal"), rtm.tabulate(q1, 2, stations, "thermal"))
    t1 = tabs[0]; gap = p["offset"]
    Istar = float(np.interp(gap, t1["u"], t1["Q"])); g1 = float(np.interp(gap, t1["u"], t1["g"])); return dict(Istar=Istar, g1=g1, tau_c=p["area"] / g1, tabs=tabs)


def filtered(ser, tau_c, I0=None, dt=60.0):
    t = np.array(ser["t"]); f = PchipInterpolator(t, ser["Q"]); tt = np.arange(t[0], t[-1] + 1e-9, dt); I = f(tt); out = np.empty_like(I); a = np.exp(-dt / tau_c)
    out[0] = I[0] if I0 is None else I0
    for k in range(1, len(tt)): out[k] = a * out[k - 1] + (1 - a) * 0.5 * (I[k] + I[k - 1])
    return tt, I, out


def steady_margin(Ibar, tabs, gap):
    """stage margin of the steady state at constant inflow Ibar relative to sill 2: positive = outlet 2 wet (u2 with Q1+Q2 = Ibar), negative = dry (s2 - H1 with Q1(H1) = Ibar)."""
    t1, t2 = tabs; Istar = float(np.interp(gap, t1["u"], t1["Q"]))
    if Ibar <= Istar:
        u = np.interp(Ibar, t1["Q"], t1["u"]); return float(u - gap)
    f = lambda u2: float(np.interp(u2 + gap, t1["u"], t1["Q"]) + np.interp(u2, t2["u"], t2["Q"])) - Ibar
    lo, hi = 0.0, 1.0
    for _ in range(60):
        mid = 0.5 * (lo + hi); lo, hi = (mid, hi) if f(mid) < 0 else (lo, mid)
    return 0.5 * (lo + hi)


def cusum_max(tt, I, Istar, sign):
    """largest water-volume excursion on the 'wrong' side of the threshold: sign=+1 deficit sum of (I* - I)+ , sign=-1 excess (I - I*)+ ; reflected at zero (Page's CUSUM)"""
    c = 0.0; m = 0.0; dt = np.diff(tt)
    for k in range(1, len(tt)):
        c = max(0.0, c + sign * (Istar - 0.5 * (I[k] + I[k - 1])) * dt[k - 1]); m = max(m, c)
    return m


def criterion(ser, th, p):
    """a-priori screen from the inflow series and the rating tables only (no time integration of the lake model).
    R = (largest wrong-side volume excursion of the inflow about I*) / (water volume that separates the steady stage from sill 2 at the mean inflow).  R >= 1: the threshold can be crossed by the fluctuations although the mean says otherwise."""
    p = rtm.full(p); tt, I, If = filtered(ser, th["tau_c"]); m = S.means(ser); Is = th["Istar"]; A = p["area"]
    um = steady_margin(m["Q_bar"], th["tabs"], p["offset"]); wet = m["Q_bar"] > Is
    D = cusum_max(tt, I, Is, +1 if wet else -1); V = A * max(abs(um), 1e-6); R = D / V
    return dict(Istar=Is, tau_c_h=th["tau_c"] / 3600, Q_bar=m["Q_bar"], ratio_mean=m["Q_bar"] / Is, I_min=float(I.min()), I_max=float(I.max()), margin_m=um, cusum_vol=D, storage_vol=V, R=R, mean_pred="wet" if wet else "dry",
                crosses=bool(I.min() < Is < I.max()), unsafe=bool(R >= 1.0),
                guaranteed=("wet" if I.min() > Is else "dry" if I.max() < Is else "indeterminate"))     # comparison theorem: inflow always above I* -> outlet 2 is eventually and permanently wet (once reached); always below -> never activates


def ea_edges(se, th, dur_s):
    """se: dict(t, wet) wet/dry indicator of outlet 2 from the water-storage model driven by the full series.  Blocks = maximal wet/dry runs; runs shorter than tau_c merged into the longer neighbour."""
    t, w = se["t"], se["wet"]; runs = []; i = 0
    while i < len(t):
        j = i
        while j < len(t) and w[j] == w[i]: j += 1
        runs.append([t[i], t[j - 1] if j < len(t) else t[-1], bool(w[i])]); i = j
    runs[0][0] = 0.0; runs[-1][1] = dur_s
    while len(runs) > 1:
        L = [r[1] - r[0] for r in runs]; k = int(np.argmin(L))
        if L[k] >= th["tau_c"]: break
        nb = [j for j in (k - 1, k + 1) if 0 <= j < len(runs)]; j = max(nb, key=lambda j: L[j]); lo, hi = min(j, k), max(j, k)
        runs[lo] = [runs[lo][0], runs[hi][1], runs[j][2]]; del runs[hi]
    return [0.0] + [float(r[1]) for r in runs[:-1]] + [float(dur_s)]


def forcing_for(method, ser, th, dur_s, se=None):
    """returns list of segments (duration, Q, T) or None for TR"""
    if method == "TR": return None
    if method == "AVG": edges = [0.0, dur_s]
    elif method == "PW6": edges = list(np.arange(0.0, dur_s + 1, 6 * 3600.0)); edges[-1] = dur_s
    elif method == "PW3": edges = list(np.arange(0.0, dur_s + 1, 3 * 3600.0)); edges[-1] = dur_s
    elif method == "EA": edges = ea_edges(se, th, dur_s)
    elif method == "PWn": n = len(ea_edges(se, th, dur_s)) - 1; edges = list(np.linspace(0.0, dur_s, n + 1))
    return S.block_series(ser, edges)
