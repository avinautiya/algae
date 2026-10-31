"""Reduced translation-invariant capture model (RTM) for a lake feeding two evolving ice channels.

Exact identity used (no closure): for outlet i with quasi-steady rating Q_i(H, geometry_i(t)), dQ_i/dt = g_i (dH/dt + E_i),
g_i = dQ_i/dH, E_i = (dQ_i/dt at fixed H)/g_i   ("rating-effective incision rate", m/s),
so  d/dt ln(Q2/Q1) = (lam2 - lam1) dH/dt + (r2 - r1),  lam_i = d ln Q_i/dH,  r_i = lam_i E_i  (fractional capacity growth rate).
RTM approximation (controlled, see docs/capture_theory.md): each route's rating and E_i are tabulated ONCE from the production hydraulics and heat-transfer kernels at the
initial geometry as functions of the head u_i over its sill; geometry then evolves only as a translation, u_i = H - s_i0 + d_i, dd_i/dt = E_i. State (H, d1, d2): no stations, no thermal march in the loop.
"""
import copy, os, sys
import numpy as np
from scipy.integrate import solve_ivp
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, "..")
sys.path.insert(0, os.path.join(ROOT, "ice_network", "src"))
from ice_network.evolution.network import Network
from ice_network.constants import LF, RHO_I

DEFAULT = dict(I=0.05, T=1.0, area=2000.0, bottom=99.0, deficit=0.05, offset=0.05, amp=0.0, period_h=3.0, n_man=0.04, ke=0.5, nu_factor=1.0, bank=1.0, crest=2.0, sill1=100.0,
               w1=0.4, L1=100.0, S1=0.005, w2=0.4, L2=60.0, S2=0.02, duration_h=8.0)


def full(p):
    q = dict(DEFAULT); q.update(p); return q


def route_net(p, k, stations=17, melting=True, melt_mode="thermal"):
    p = full(p); w, L, S = p[f"w{k}"], p[f"L{k}"], p[f"S{k}"]; sill = p["sill1"] + (p["offset"] if k == 2 else 0.0)
    x = np.linspace(0.0, L, stations); z = sill - S * np.maximum(x - p["crest"], 0.0)
    e = dict(a=0, b=1, z=z, w=np.full(stations, w), zs=z + p["bank"], length=L, msub=16, n_man=p["n_man"], ke=p["ke"])
    net = Network([dict(z_bottom=p["bottom"], a_ref=p["area"], h_ref=1.0, q=0.0)], [dict(h_fixed=float(z[-1]) - 5.0)], [e], melting=melting, nu_factor=p["nu_factor"], melt_mode=melt_mode)
    return net, sill


def tabulate(p, k, stations=17, mode="thermal", ugrid=None, dt_pert=10.0):
    """rating Q(u), g(u)=dQ/dH, E(u) for route k at its initial geometry; mode: thermal | dissipation | fixed"""
    p = full(p)
    net, sill = route_net(p, k, stations, melting=(mode != "fixed"), melt_mode="dissipation_all" if mode == "dissipation" else "thermal")
    h_out = net.outlets[0]["h_fixed"]; T = p["T"]
    ug = np.concatenate([[0.0], np.geomspace(2e-4, 1.2, 70)]) if ugrid is None else ugrid
    Q = np.zeros_like(ug); g = np.zeros_like(ug); E = np.zeros_like(ug); ms = np.zeros_like(ug)   # ms: inlet-sill (station 0) bed-lowering rate, what the network and the reference report as sill lowering
    for i, u in enumerate(ug):
        H = sill + u; st = net.edge_state(net.edges[0], H, h_out, T); Q[i] = abs(st["Q"])
        if Q[i] <= 0.0: continue
        dh = 1e-6 * max(1.0, u)
        Qp = abs(net.edge_state(net.edges[0], H + dh, h_out, T)["Q"]); Qm = abs(net.edge_state(net.edges[0], max(H - dh, sill), h_out, T)["Q"]) if u > dh else Q[i]
        g[i] = (Qp - Qm) / (dh + (dh if u > dh else 0.0))
        if mode != "fixed":
            e0 = net.edges[0]; wm = 0.5 * (e0["w"][0] + e0["w"][1]); ms[i] = st["q_ice"][0] * (1.0 - st["wall_frac"][0]) / (RHO_I * LF * wm * e0["dx"])
        if mode == "fixed" or g[i] <= 0: continue
        n2 = copy.deepcopy(net); hs = st["q_ice"] * dt_pert; ws = hs * st["wall_frac"]
        n2.update_geometry(hs, ws)
        Q2 = abs(n2.edge_state(n2.edges[0], H, h_out, T)["Q"]); E[i] = (Q2 - Q[i]) / dt_pert / g[i]
    return dict(u=ug, Q=Q, g=g, E=E, ms=ms, sill=sill)


def _interp(tab, key, u):
    return float(np.interp(u, tab["u"], tab[key])) if u > 0 else 0.0


def rtm_run(p, mode="thermal", stations=17, tabs=None, horizon_h=None, max_step=300.0):
    p = full(p); T_end = (horizon_h or p["duration_h"]) * 3600.0
    tabs = tabs or (tabulate(p, 1, stations, mode), tabulate(p, 2, stations, mode)); t1, t2 = tabs
    inflow = (lambda t: p["I"]) if p["amp"] == 0 else (lambda t: p["I"] * (1 + p["amp"] * np.sin(2 * np.pi * t / (p["period_h"] * 3600.0))))
    H0 = p["sill1"] - p["deficit"]
    def rhs(t, y):
        H, d1, d2, _, __ = y; u1 = H - t1["sill"] + d1; u2 = H - t2["sill"] + d2
        q1, q2 = _interp(t1, "Q", u1), _interp(t2, "Q", u2)
        return [(inflow(t) - q1 - q2) / p["area"], _interp(t1, "E", u1) if q1 > 0 else 0.0, _interp(t2, "E", u2) if q2 > 0 else 0.0, _interp(t1, "ms", u1) if q1 > 0 else 0.0, _interp(t2, "ms", u2) if q2 > 0 else 0.0]
    sol = solve_ivp(rhs, (0, T_end), [H0, 0.0, 0.0, 0.0, 0.0], method="LSODA", rtol=1e-8, atol=1e-10, max_step=max_step, dense_output=True)
    t = np.arange(0.0, T_end + 1e-9, 60.0); Y = sol.sol(t)
    u1 = Y[0] - t1["sill"] + Y[1]; u2 = Y[0] - t2["sill"] + Y[2]
    Q1 = np.array([_interp(t1, "Q", u) for u in u1]); Q2 = np.array([_interp(t2, "Q", u) for u in u2])
    return dict(t=t, H=Y[0], d1=Y[1], d2=Y[2], s1=Y[3], s2=Y[4], Q1=Q1, Q2=Q2, tabs=tabs, success=bool(sol.success))


def events(res, I):
    t, Q1, Q2 = res["t"], res["Q1"], res["Q2"]; thr = 0.01 * I
    def cross(y, th):
        i = np.flatnonzero((y[:-1] <= th) & (y[1:] > th)); return None if not i.size else float(t[i[0]] + (th - y[i[0]]) / (y[i[0] + 1] - y[i[0]]) * 60.0)
    dom = cross(Q2 - Q1, 0.0)
    return dict(t_act1_s=cross(Q1, thr), t_act2_s=cross(Q2, thr), t_dom2_s=dom, d1_m=float(res["d1"][-1]), d2_m=float(res["d2"][-1]), sill1_m=float(res["s1"][-1]), sill2_m=float(res["s2"][-1]), gap_drift_mm=float(1e3 * ((res["s2"][-1] - res["s1"][-1]) - (res["s2"][len(res["t"]) // 2] - res["s1"][len(res["t"]) // 2])) * -1.0),
                delivered1_m3=float(np.trapezoid(Q1, t)), delivered2_m3=float(np.trapezoid(Q2, t)), phi2_end=float(Q2[-1] / max(Q1[-1] + Q2[-1], 1e-12)))


def partition_at_gap(tabs, I, gap):
    """quasi-steady partition for sill gap = s2 - s1 (>= or <= 0): returns u1, Q1, Q2, E1, E2, g1, g2"""
    t1, t2 = tabs
    f = lambda u1: np.interp(u1, t1["u"], t1["Q"]) + (np.interp(u1 - gap, t2["u"], t2["Q"]) if u1 - gap > 0 else 0.0) - I
    lo, hi = max(gap, 0.0) + 1e-9, 1.2
    if f(hi) < 0: return None
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0: hi = mid
        else: lo = mid
    u1 = 0.5 * (lo + hi); u2 = u1 - gap
    q1, q2 = np.interp(u1, t1["u"], t1["Q"]), (np.interp(u2, t2["u"], t2["Q"]) if u2 > 0 else 0.0)
    return dict(u1=u1, u2=u2, Q1=q1, Q2=q2, E1=np.interp(u1, t1["u"], t1["E"]), E2=np.interp(u2, t2["u"], t2["E"]) if u2 > 0 else 0.0,
                g1=np.interp(u1, t1["u"], t1["g"]), g2=np.interp(u2, t2["u"], t2["g"]) if u2 > 0 else 0.0)


def gap_dynamics(tabs, I, gaps=None):
    """phase line of the sill gap Delta = s2 - s1 in the quasi-steady-lake limit: dDelta/dt = E1 - E2 (exact kinematics). Returns Delta grid, G = dDelta/dt,
    kappa (gap at which Q1 = Q2), fixed points (zero crossings of G) with stability, and the switching-time quadrature."""
    gaps = np.linspace(-0.12, 0.30, 421) if gaps is None else gaps
    G = np.full(gaps.size, np.nan); phi2 = np.full(gaps.size, np.nan)
    for i, d in enumerate(gaps):
        s = partition_at_gap(tabs, I, d)
        if s: G[i] = s["E1"] - s["E2"]; phi2[i] = s["Q2"] / (s["Q1"] + s["Q2"])
    ok = np.isfinite(G); kappa = None
    if ok.any():
        j = np.flatnonzero(ok[:-1] & ok[1:] & ((phi2[:-1] - 0.5) * (phi2[1:] - 0.5) < 0))
        if j.size: i = j[0]; kappa = float(gaps[i] + (0.5 - phi2[i]) / (phi2[i + 1] - phi2[i]) * (gaps[i + 1] - gaps[i]))
    fps = []
    for i in np.flatnonzero(ok[:-1] & ok[1:] & (G[:-1] * G[1:] < 0)):
        d0 = float(gaps[i] + G[i] / (G[i] - G[i + 1]) * (gaps[i + 1] - gaps[i])); fps.append(dict(gap=d0, stable=bool((G[i + 1] - G[i]) < 0), slope=float((G[i + 1] - G[i]) / (gaps[i + 1] - gaps[i]))))
    return dict(gaps=gaps, G=G, phi2=phi2, kappa=kappa, fixed_points=fps)


def switch_time_qs(pd, delta0, kappa):
    """time for the sill gap to close from delta0 to kappa under dDelta/dt = G(Delta); None when G >= 0 anywhere on the path (no overtaking in the quasi-steady-lake limit)"""
    if kappa is None or delta0 <= kappa: return None
    g, G = pd["gaps"], pd["G"]; m = (g >= kappa) & (g <= delta0) & np.isfinite(G)
    if m.sum() < 3 or np.any(G[m] >= 0): return None
    return float(np.trapezoid(1.0 / -G[m], g[m]))
