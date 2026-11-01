"""Reduced comparators for the CFD campaign, parameterised to the ACTUAL CFD geometry (Q=0.07, T=0.5, lane width 0.5 m, reach 2 m, toe floor -0.08, bypass 0, deposit volume 0.1 m3 over 0.5 m2 => equivalent uniform depth 0.2 m, cover 0.5, h_eff=200, k_debris=2).
Imported package code (pinned commit 3e51620) is used unmodified for the connectivity-aware model (topology.run). Two further comparators are implemented here with the SAME equations:
  fractional: one toe control with mean crest (floor + cover*depth) and area-weighted shield (the 'uniform/fractional-cover approximation');
  fixed: geometry frozen at t=0 (no incision; contact never changes).
Inputs identical across methods; no information from CFD runs is used (only the geometry spec)."""
import math, os, sys
import numpy as np
PKG = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "imported", "ice_stream_reactivation_3e51620", "compute")
sys.path.insert(0, os.path.abspath(PKG))
from model import Case          # noqa: E402  (pinned import)
import topology as T            # noqa: E402
def case(): return Case(discharge=.07, allow_supply=False)
def _contact(hist, Q, key_toe):
    return [dict(day=h["day"], f_toe=float(sum(h["flow"][i] for i in key_toe) / Q)) for h in hist]
def connectivity(layout, removal_day, days=21.0, c=None):
    c = c or case(); r = T.run(c, "corridor" if layout == "corridor" else "spanning_strip", cover_fraction=.5, grain_depth=.2, removal_day=removal_day, days=days, dt=225.)
    toe = [0, 1] if layout == "corridor" else [0]; series = _contact(r["history"], c.discharge, toe)
    return dict(method="connectivity_reduced", series=series, toe_melt_m3=r["toe_melt_m3"], energy_resid=r["energy_relative_residual"], final_ice_bed=r["final_ice_bed"])
def _weir_eta(Q, crest, K):
    lo, hi = float(np.min(crest)), float(np.max(crest) + (Q / np.min(K)) ** (2 / 3))
    for _ in range(60):
        eta = (lo + hi) / 2
        if float(np.sum(K * np.maximum(eta - crest, 0) ** 1.5)) < Q: lo = eta
        else: hi = eta
    return (lo + hi) / 2
def fractional(removal_day, days=21.0, evolve=True, cover=.5, depth=.2, dt=225.0, c=None):
    """toe control (width W) with mean-height crest + area-weighted shield; bypass bare."""
    c = c or case(); W, L = c.width, c.length; K = (2 / 3) ** 1.5 * math.sqrt(c.gravity) * np.array([W, W]); bed = np.array([c.initial_toe_ice_z, c.initial_bypass_z]); Q = c.discharge; lat = c.ice_density * c.latent; melt = np.zeros(2); hist = []; A = W * L
    n = math.ceil(days * 86400 / dt)
    for k in range(n + 1):
        t = k * dt; present = removal_day is None or t < removal_day * 86400
        crest = bed + np.array([cover * depth * present, 0.0]); eta = _weir_eta(Q, crest, K); q = K * np.maximum(eta - crest, 0) ** 1.5
        if k % int(21600 / dt) == 0: hist.append(dict(day=t / 86400, f_toe=float(q[0] / Q)))
        if k == n: break
        he = np.array([(1 - cover * present) * c.h_water + cover * present / (1 / c.h_water + depth / c.debris_k), c.h_water]); cap = c.water_density * c.water_cp * q
        ntu = np.divide(he * A, cap, out=np.zeros(2), where=cap > 0); P = cap * c.water_temperature * (-np.expm1(-ntu)); m = P / (lat * A)
        if evolve: bed = bed - dt * m
        melt += dt * m * A
    return dict(method="fractional_cover" if evolve else "fixed_geometry", series=hist, toe_melt_m3=float(melt[0]), final_ice_bed=bed.tolist())
def fixed(removal_day, days=21.0, **kw):
    r = fractional(removal_day, days=days, evolve=False, **kw); return r
def all_methods(layout, removal_day, days=21.0):
    out = [connectivity(layout, removal_day, days), fractional(removal_day, days), fixed(removal_day, days)]; return out
