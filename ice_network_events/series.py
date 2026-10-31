"""Hourly melt-water forcing series (diurnal cycle x multi-hour spell) with enthalpy-consistent temperature, and its matched means.
Q(t) = Qm [1 + a sin(2 pi (t - phase)/24 h)] [1 - depth * spell(t)], spell = smooth boxcar (start, length); T(t) = Tm (Q/Qm)^gamma (warm water when flow is high).
Samples are hourly (t, Q, H = Q T); the simulators join them with a monotone C1 (PCHIP) interpolant, so water and heat budgets are exactly what the series states."""
import numpy as np
from scipy.interpolate import PchipInterpolator
HR = 3600.0


def make_series(f, dur_h=36.0):
    t = np.arange(0.0, dur_h + 1e-9, 1.0)
    sp = np.zeros_like(t)
    if f.get("depth", 0) > 0:
        a, b = f["spell_start_h"], f["spell_start_h"] + f["spell_len_h"]; w = 1.0
        sp = 0.5 * (np.tanh((t - a) / w) - np.tanh((t - b) / w))
    Q = f["Qm"] * (1 + f["amp"] * np.sin(2 * np.pi * (t - f["phase_h"]) / 24.0)) * (1 - f.get("depth", 0) * sp)
    Q = np.maximum(Q, 1e-4); T = f["Tm"] * (Q / f["Qm"]) ** f.get("gamma", 0.5)
    return dict(t=(t * HR).tolist(), Q=Q.tolist(), H=(Q * T).tolist())


def means(ser):
    t = np.array(ser["t"]); Q = np.array(ser["Q"]); H = np.array(ser["H"])
    f = PchipInterpolator(t, Q); g = PchipInterpolator(t, H); V = float(f.integrate(t[0], t[-1])); Hq = float(g.integrate(t[0], t[-1])); D = t[-1] - t[0]
    return dict(Q_bar=V / D, T_bar=Hq / V, duration=D, volume=V)


def block_series(ser, edges_s):
    """replace the series by constant block means between edges (s), water- and heat-conserving per block; blocks joined with 900 s raised-cosine ramps as in the sequence study.
    Returns segments [(duration, Q, T)]"""
    t = np.array(ser["t"]); f = PchipInterpolator(t, ser["Q"]); g = PchipInterpolator(t, ser["H"]); seg = []
    for a, b in zip(edges_s[:-1], edges_s[1:]):
        V = float(f.integrate(a, b)); Hq = float(g.integrate(a, b)); seg.append((b - a, V / (b - a), Hq / max(V, 1e-12)))
    return seg
