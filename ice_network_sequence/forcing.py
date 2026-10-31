"""Forcing episodes with matched water and sensible-heat budgets.
An experiment is a list of constant-forcing segments [(duration_s, Q m3/s, T C), ...] joined by raised-cosine ramps of width ramp_s centred on each boundary (smooth RHS for the integrators;
the ramp is symmetric, so reversing the episode order leaves int Q dt and int Q T dt unchanged). The averaged control uses discharge-weighted (enthalpy-conserving) temperature
T_bar = int Q T dt / int Q dt, constant discharge Q_bar = int Q dt / duration. Budgets are verified numerically on the actual functions (verify_budget)."""
import numpy as np


def _step(x):    # smooth 0 -> 1 over x in [-0.5, 0.5]
    x = np.clip(x + 0.5, 0.0, 1.0); return 0.5 - 0.5 * np.cos(np.pi * x)


def make(segments, ramp_s=900.0):
    ends = np.cumsum([s[0] for s in segments])[:-1]
    Q = [s[1] for s in segments]; H = [s[1] * s[2] for s in segments]      # ramp the water flux Q and the heat flux H = Q T linearly in the smooth step: both budgets are then exactly order-independent
    def f(vals):
        def g(t):
            t = np.asarray(t, float); out = np.full(t.shape, vals[0], float)
            for k, tb in enumerate(ends): out = out + (vals[k + 1] - vals[k]) * _step((t - tb) / ramp_s)
            return out
        return g
    Qf, Hf = f(Q), f(H)
    return Qf, (lambda t: Hf(t) / Qf(t))


def budgets(Qf, Tf, t_end, n=200001):
    t = np.linspace(0.0, t_end, n); q = Qf(t); T = Tf(t)
    return dict(water_m3=float(np.trapezoid(q, t)), heat_J_over_rho_c=float(np.trapezoid(q * T, t)))   # rho_w c_w * this = incoming sensible heat (J)


def sequences(A, B, tau_A, tau_B, tail, ramp_s=900.0):
    """A, B = (Q, T). Returns dict of segment lists for AB, BA and the averaged control, plus the common continuation (tail duration, mean forcing).
    The mean over the two episodes is used for the continuation and the control so that the full-experiment budgets agree."""
    tot = tau_A + tau_B + tail
    # mean forcing of the two episodes (including the ramp): water mean and discharge-weighted temperature
    seg0 = [(tau_A, *A), (tau_B, *B)]; Qf, Tf = make(seg0, ramp_s); b = budgets(Qf, Tf, tau_A + tau_B)
    qm = b["water_m3"] / (tau_A + tau_B); tm = b["heat_J_over_rho_c"] / b["water_m3"]
    seq = dict(AB=[(tau_A, *A), (tau_B, *B), (tail, qm, tm)], BA=[(tau_B, *B), (tau_A, *A), (tail, qm, tm)], AVG=[(tot, qm, tm)])
    return seq, dict(Q_bar=qm, T_bar=tm, duration=tot)


def verify_budget(seq, ramp_s=900.0):
    out = {}
    for k, s in seq.items():
        Qf, Tf = make(s, ramp_s); out[k] = budgets(Qf, Tf, sum(x[0] for x in s))
    return out
