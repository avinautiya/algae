"""Verification gate for the plan-view solver (M3). `python -m ice_cliff_survival.reactivation.verify_swe OUT.json`.
Tolerances (declared BEFORE judging production output): rest: max|u| < 1e-10 m/s and stage deviation < 1e-12 m (well-balancedness is exact in floating point); mass/energy closure (with storage): < 1e-10 of inflow;
flow partition vs the independent critical-weir formula: |f_toe - f_ref| < 0.05 absolute (the formula ignores friction/approach velocity; the friction head-loss scale is reported) and the dx-series error must not grow;
heat plug-flow limit: relative error of lane melt heat < 5 % and decreasing with dx; bed-lowering vs heat removed from the water (independent boundary-flux measurement): < 1e-6 relative; spatial refinement: report observed order, require |f(dx)-f(dx/2)| to decrease."""
import json, math, sys, time
import numpy as np
from . import cases
from .swe2d import Sim, advance, RHO_W, CP, RHO_I, LF, G
def rest():
    s = cases.build("strip", dx=0.05, tail_slope=0.0, Q=0.0); s.closed_out = True; s.reset_state(0.05); eta0 = (s.h + s.zs() - 0)[1:-1, 1:-1].copy(); mask = s.h[1:-1, 1:-1] > 0
    qacc = np.zeros_like(s.zb); fa = np.zeros(5); V0 = s.h.sum()
    advance(s.h, s.hu, s.hv, s.hT, s.zs(), s.pdep, s.heff_field(), s.ice, s.manning, 0.0, 0.0, s.dx, s.nx, s.ny, 30.0, 0.4, 0.0, 0.4, 0.2, 1e-6, False, qacc, fa, True)
    eta1 = (s.h + s.zs())[1:-1, 1:-1]; u = np.abs(s.hu[1:-1, 1:-1]).max() + np.abs(s.hv[1:-1, 1:-1]).max()
    return dict(test="lake at rest with partly dry deposit", max_momentum=float(u), stage_dev=float(np.abs(eta1 - eta0)[mask].max()), volume_change_rel=float(abs(s.h.sum() - V0) / V0), pass_=bool(u < 1e-10 and np.abs(eta1 - eta0)[mask].max() < 1e-12))
def weir_ref(Q, crests, widths):
    K = np.array(widths) * (2 / 3) ** 1.5 * math.sqrt(G); crests = np.array(crests); lo, hi = crests.min(), crests.max() + (Q / K.min()) ** (2 / 3) + 1.0
    for _ in range(80):
        eta = 0.5 * (lo + hi)
        if (K * np.maximum(eta - crests, 0) ** 1.5).sum() < Q: lo = eta
        else: hi = eta
    q = K * np.maximum(eta - crests, 0) ** 1.5; return eta, q / Q
def gvf_ref(Q, lanes, L=2.0, n=0.02):
    """independent reference WITH friction: each lane is a prismatic flat channel (floor z, width W) fed from a pool at stage eta (no entrance loss) and ending in a critical brink; Manning gradually varied flow integrated upstream from the brink. lanes=[(z,W),...]. Returns eta, flow fractions."""
    from scipy.integrate import solve_ivp
    from scipy.optimize import brentq
    def head(q, z):          # specific energy + z at the lane entrance for unit discharge q
        hc = (q * q / G) ** (1 / 3)
        def rhs(x, h): return -(0.0 - n * n * q * q / h ** (10 / 3)) / (1 - q * q / (G * h ** 3))      # integrate in x measured upstream from the brink: dh/dxu = -dh/dx
        sol = solve_ivp(rhs, [0, L], [hc * 1.0005], rtol=1e-10, atol=1e-13); h = sol.y[0, -1]; return z + h + q * q / (2 * G * h * h)
    def lane_q(eta, z, W):
        if eta <= z: return 0.0
        f = lambda q: head(q, z) - eta
        qmax = 1.2 * math.sqrt(G) * (eta - z) ** 1.5
        return W * brentq(f, 1e-9, qmax) if f(1e-9) < 0 < f(qmax) else W * qmax
    tot = lambda eta: sum(lane_q(eta, z, W) for z, W in lanes) - Q
    zmin = min(z for z, _ in lanes); eta = brentq(tot, zmin + 1e-4, zmin + 2.0, xtol=1e-12); qs = np.array([lane_q(eta, z, W) for z, W in lanes]); return float(eta), qs / Q
def column_ref(s, Q=0.07):
    """independent reference: every y-column of both lanes is a parallel broad-crested weir with crest = max over the reach of (floor + solid deposit height) in that column, width dx; shared pool stage; lateral exchange and friction ignored (friction checked separately by gvf_ref)."""
    zs = s.zs(); cols = list(range(1, s.nw + 1)) + list(range(s.nw + 2, s.ny + 1)); crests = [float(zs[s.i0:s.i1 + 1, j].max()) for j in cols]
    eta, f = weir_ref(Q, crests, [s.dx] * len(cols)); ft = float(f[:s.nw].sum()); return eta, ft
def partition(dxs=(0.1, 0.05, 0.025), layouts=("bare", "strip", "corridor"), **kw):
    out = []
    for lay in layouts:
        for dx in dxs:
            s = cases.build(lay, dx=dx, **kw); r = s.hydraulics(tmax=4000.0, cool=False); d = s.diagnostics(); eta_ref, fr = column_ref(s)
            out.append(dict(layout=lay, dx=dx, f_toe=d["f_toe"], f_ref=fr, err=float(d["f_toe"] - fr), eta_ref=eta_ref, converged=r["converged"], hyd_time=r["hyd_time"], dt_limit=r["dt_limit"], water_resid=float(d["water_resid"]), energy_resid=float(d["energy_resid_rel"])))
    return out
def plug(dxs=(0.1, 0.05, 0.025)):
    out = []
    for dx in dxs:
        s = Sim(dx=dx, z_toe=0.0, z_byp=0.0, heff=20000.0, Tw=0.5, tail_slope=0.5); r = s.hydraulics(tmax=4000.0); d = s.diagnostics(); rows = []
        for ys, nm in ((s.toe_cols(), "toe"), (s.byp_cols(), "byp")):
            Qj = s.Qsec(ys); ntu = s.heff0 * 0.5 * 2.0 / (RHO_W * CP * Qj); P = RHO_W * CP * Qj * 0.5 * (1 - math.exp(-ntu)); Pm = d["melt_heat_" + nm + "_W"]
            # bed lowering profile after dt: analytic dz(x) = heff T(x) dt/(rho_i L), T = Tw exp(-heff x/(rho cp q_w))
            dt = 3600.0; qw = Qj / 0.5; xs = (np.arange(s.i0, s.i1 + 1) - s.i0 + 0.5) * dx; dz_an = s.heff0 * 0.5 * np.exp(-s.heff0 * xs / (RHO_W * CP * qw)) * dt / (RHO_I * LF)
            dzn = (s.last_qacc[s.i0:s.i1 + 1, ys] * dt / (RHO_I * LF)).mean(axis=1); rows.append(dict(lane=nm, NTU=ntu, melt_W=Pm, analytic_W=P, rel_err=(Pm - P) / P, lowering_profile_max_rel_err=float(np.max(np.abs(dzn / dz_an - 1)))))
        out.append(dict(dx=dx, lanes=rows))
    return out
def bed_vs_water(dx=0.05):
    s = cases.build("strip", dx=dx); s.hydraulics(tmax=4000.0); d = s.diagnostics(); dtm = 3600.0; vol, vt, vb = s.lower_bed(dtm); latent_W = RHO_I * LF * vol / dtm
    water_W = d["heat_in_W"] - d["heat_out_W"] - d["storage_rate_W"]; return dict(latent_W=latent_W, water_side_W=water_W, rel_diff=(latent_W - water_W) / water_W, storage_W=d["storage_rate_W"], sensible_in_W=d["heat_in_W"], sensible_out_W=d["heat_out_W"])
def outlet():
    out = []
    for Lt in (0.5, 1.0, 1.5):
        s = cases.build("strip", dx=0.05, Ltail=Lt); s.hydraulics(tmax=4000.0); d = s.diagnostics(); out.append(dict(Ltail=Lt, f_toe=d["f_toe"], melt_toe_W=d["melt_heat_toe_W"], melt_byp_W=d["melt_heat_byp_W"]))
    return out
def order(a, b, c): return float(math.log2(abs(a - b) / abs(b - c))) if abs(b - c) > 0 and abs(a - b) > 0 else float("nan")
if __name__ == "__main__":
    t = time.time(); R = dict(rest=rest(), partition=partition(), plug=plug(), bed_vs_water=bed_vs_water(), outlet=outlet())
    P = R["partition"]
    for lay in ("bare", "strip", "corridor"):
        f = [p["f_toe"] for p in P if p["layout"] == lay]; R.setdefault("order_f_toe", {})[lay] = order(*f)
    R["wall_s"] = time.time() - t; json.dump(R, open(sys.argv[1], "w"), indent=1, default=float); print(json.dumps({k: R[k] for k in ("rest", "order_f_toe", "bed_vs_water", "outlet")}, indent=1, default=float))
    for p in P: print({k: (round(v, 5) if isinstance(v, float) else v) for k, v in p.items()})
    for p in R["plug"]: print(p)
