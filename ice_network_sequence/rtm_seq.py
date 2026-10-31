"""RTM with time-varying inflow/temperature and a well-mixed lake that stores water AND heat (state H, S = V T_L, d1, d2, s1, s2, Qheat1, Qheat2).
Rates are tabulated once at T = 1 C and scaled linearly with the lake temperature T_L (checked against the network in the development set): E_i(u, T_L) = T_L E_i^(1)(u). Dissipation (about 0.3 % of the thermal melt at 1 C) is dropped; meltwater input to the lake is not represented (same as RTM)."""
import os, sys
import numpy as np
from scipy.integrate import solve_ivp
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "ice_network_capture"))
import rtm, forcing


def run(p, segments, mode="thermal", stations=17, tabs=None, ramp_s=900.0, t_end=None, lake_T0=None, fns=None):
    p = rtm.full(p); T_end = t_end or sum(s[0] for s in segments); Qf, Tf = fns if fns is not None else forcing.make(segments, ramp_s)
    if tabs is None:
        q1 = dict(p, T=1.0); tabs = (rtm.tabulate(q1, 1, stations, "fixed" if mode == "fixed" else "thermal"), rtm.tabulate(q1, 2, stations, "fixed" if mode == "fixed" else "thermal"))
    t1, t2 = tabs; A = p["area"]; V0 = A * (p["sill1"] - p["deficit"] - p["bottom"]); T0 = Tf(0.0) if lake_T0 is None else lake_T0
    ip = lambda tab, k, u: float(np.interp(u, tab["u"], tab[k])) if u > 0 else 0.0
    def rhs(t, y):
        V, S, d1, d2, s1, s2, h1, h2, c1, c2 = y; V = max(V, 1e-6); TL = S / V; H = p["bottom"] + V / A
        u1 = H - t1["sill"] + d1; u2 = H - t2["sill"] + d2; q1 = ip(t1, "Q", u1); q2 = ip(t2, "Q", u2); I = float(Qf(t)); Ti = float(Tf(t))
        e1 = TL * ip(t1, "E", u1) if q1 > 0 else 0.0; e2 = TL * ip(t2, "E", u2) if q2 > 0 else 0.0
        m1 = TL * ip(t1, "ms", u1) if q1 > 0 else 0.0; m2 = TL * ip(t2, "ms", u2) if q2 > 0 else 0.0
        return [I - q1 - q2, I * Ti - (q1 + q2) * TL, e1, e2, m1, m2, TL * ip(t1, "qice", u1) if q1 > 0 else 0.0, TL * ip(t2, "qice", u2) if q2 > 0 else 0.0, q1, q2]
    sol = solve_ivp(rhs, (0, T_end), [V0, V0 * T0, 0, 0, 0, 0, 0, 0, 0, 0], method="LSODA", rtol=1e-8, atol=1e-9, max_step=300.0, dense_output=True)
    t = np.arange(0.0, T_end + 1e-9, 60.0); Y = sol.sol(t); V = np.maximum(Y[0], 1e-6); H = p["bottom"] + V / A
    u1 = H - t1["sill"] + Y[2]; u2 = H - t2["sill"] + Y[3]
    Q1 = np.array([ip(t1, "Q", u) for u in u1]); Q2 = np.array([ip(t2, "Q", u) for u in u2])
    return dict(t=t, H=H, V=Y[0], TL=Y[1] / V, d1=Y[2], d2=Y[3], s1=Y[4], s2=Y[5], heat1=Y[6], heat2=Y[7], cum1=Y[8], cum2=Y[9], Q1=Q1, Q2=Q2, success=bool(sol.success), tabs=tabs)
