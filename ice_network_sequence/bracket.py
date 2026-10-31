"""Lie-bracket expansion of the sequence difference for the reduced state x = (V, T_L, d1, d2) of RTM_seq with constant-forcing vector fields f_A, f_B:
   Phi_B(Phi_A(x0)) - Phi_A(Phi_B(x0)) = tau_A tau_B [f_A, f_B](x0) + O(tau^3),   [f_A, f_B] = Df_B f_A - Df_A f_B,
valid for smooth fields (no wet/dry or activation crossing within the episodes), instantaneous switches, tau much smaller than every state relaxation time. Jacobians by central differences.
f(x; Q, T) = ( Q - Q1(u1) - Q2(u2),  Q (T - T_L)/V,  T_L E1^(1)(u1),  T_L E2^(1)(u2) ),  u_i = H - s_i0 + d_i, H = z_b + V/A."""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "ice_network_capture"))
import rtm, rtm_seq, forcing


def make_field(p, tabs):
    p = rtm.full(p); t1, t2 = tabs; A = p["area"]
    ip = lambda tab, k, u: float(np.interp(u, tab["u"], tab[k])) if u > 0 else 0.0
    def f(x, Q, T):
        V, TL, d1, d2 = x; H = p["bottom"] + V / A; u1 = H - t1["sill"] + d1; u2 = H - t2["sill"] + d2; q1, q2 = ip(t1, "Q", u1), ip(t2, "Q", u2)
        return np.array([Q - q1 - q2, Q * (T - TL) / V, TL * ip(t1, "E", u1) if q1 > 0 else 0.0, TL * ip(t2, "E", u2) if q2 > 0 else 0.0])
    return f


def jac(f, x, Q, T, h=None):
    x = np.asarray(x, float); n = x.size; J = np.zeros((n, n)); h = h if h is not None else np.array([1e-3 * x[0], 1e-4, 1e-6, 1e-6])
    for j in range(n):
        e = np.zeros(n); e[j] = h[j]; J[:, j] = (f(x + e, Q, T) - f(x - e, Q, T)) / (2 * h[j])
    return J


def bracket(f, x0, A, B):
    fA, fB = f(x0, *A), f(x0, *B); JA, JB = jac(f, x0, *A), jac(f, x0, *B)
    return JB @ fA - JA @ fB, dict(fA=fA, fB=fB, JA=JA, JB=JB)


def exact_difference(p, tabs, A, B, tau, x0, ramp_s=30.0):
    """RTM_seq run of A then B (and B then A) with constant episodes of length tau and short ramps; returns the state difference (V, T_L, d1, d2) at t = 2 tau"""
    out = {}
    for name, seg in (("AB", [(tau, *A), (tau, *B)]), ("BA", [(tau, *B), (tau, *A)])):
        r = rtm_seq.run(dict(p), seg, tabs=tabs, ramp_s=ramp_s, lake_T0=x0[1]); out[name] = np.array([r["V"][-1], r["TL"][-1], r["d1"][-1], r["d2"][-1]])
    return out["AB"] - out["BA"]


if __name__ == "__main__":
    p = dict(offset=0.05, deficit=-0.10); q1 = rtm.full(dict(p, T=1.0)); tabs = (rtm.tabulate(q1, 1, 17, "thermal"), rtm.tabulate(q1, 2, 17, "thermal"))
    A, B = (0.03, 2.0), (0.09, 0.2); pf = rtm.full(p); V0 = pf["area"] * (pf["sill1"] - pf["deficit"] - pf["bottom"]); x0 = np.array([V0, 0.65, 0.0, 0.0])
    f = make_field(p, tabs); br, parts = bracket(f, x0, A, B)
    print("bracket [fA,fB] =", br)
    for tau in (300, 600, 1200, 2400, 4800, 9600):
        d = exact_difference(p, tabs, A, B, float(tau), x0); pred = tau * tau * br
        print("tau %5d s  exact dV %.3e dTL %.3e dd1 %.3e dd2 %.3e | predicted %.3e %.3e %.3e %.3e | ratio d2 %.3f" % (tau, d[0], d[1], d[2], d[3], pred[0], pred[1], pred[2], pred[3], d[3] / pred[3]))
