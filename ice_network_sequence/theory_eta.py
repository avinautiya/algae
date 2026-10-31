"""Closed-form sequence effect on the heat delivered to the channels (thermal-storage mechanism).
Assumptions (stated in docs/sequence_theory.md): quasi-steady lake level (outflow = inflow), constant lake volume V, heat to ice at rate T_L kappa(Q) with kappa tabulated from the frozen-geometry quasi-steady partition for the discharge Q,
lake temperature obeys V dT_L/dt = Q (T_in - T_L) (exact for a well-mixed lake). For piecewise-constant episodes the integral of T_L over each episode is explicit.
eta = (H_AB - H_BA) / mean(H_AB, H_BA)."""
import os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "ice_network_capture"))
import rtm


def kappa_of_Q(p, tabs, Q):
    """heat delivered to the ice per degree of lake temperature (W/K x rho_w c_w normalisation absorbed: units of tab 'qice') at the quasi-steady partition of discharge Q (frozen geometry)"""
    s = rtm.partition_at_gap(tabs, Q, p["offset"])
    if s is None: return 0.0
    t1, t2 = tabs; return float(np.interp(s["u1"], t1["u"], t1["qice"]) + (np.interp(s["u2"], t2["u"], t2["qice"]) if s["u2"] > 0 else 0.0))


def heat_sequence(order, A, B, tau, V, T0, kA, kB, qm, tm, km, tail):
    """order 'AB' or 'BA'; A, B = (Q, T_in); returns total heat integral over episodes + tail and the lake temperature at the end"""
    TL = T0; H = 0.0
    segs = [(tau, *A, kA), (tau, *B, kB)] if order == "AB" else [(tau, *B, kB), (tau, *A, kA)]
    segs.append((tail, qm, tm, km))
    for dur, Q, Tin, k in segs:
        th = V / Q; e = np.exp(-dur / th); H += k * (Tin * dur + (TL - Tin) * th * (1 - e)); TL = Tin + (TL - Tin) * e
    return H, TL


def eta(p, pair, tau, tail, V, T0=None):
    A, B = pair; qm = (A[0] + B[0]) / 2; tm = (A[0] * A[1] + B[0] * B[1]) / (A[0] + B[0]); T0 = tm if T0 is None else T0
    q1 = dict(rtm.full(p), T=1.0); tabs = (rtm.tabulate(q1, 1, 17, "thermal"), rtm.tabulate(q1, 2, 17, "thermal")); pp = rtm.full(p)
    kA, kB, km = (kappa_of_Q(pp, tabs, Q) for Q in (A[0], B[0], qm))
    HAB, TAB = heat_sequence("AB", A, B, tau, V, T0, kA, kB, qm, tm, km, tail); HBA, TBA = heat_sequence("BA", A, B, tau, V, T0, kA, kB, qm, tm, km, tail)
    return dict(eta=float((HAB - HBA) / (0.5 * (HAB + HBA))), H_AB=float(HAB), H_BA=float(HBA), TL_end_AB=float(TAB), TL_end_BA=float(TBA), kappa=[kA, kB, km], mean=[qm, tm])
