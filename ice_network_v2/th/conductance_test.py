"""Does the CFD-derived conductance satisfy the common-shape assumption of ice_incision_forecast (commit 471415c, THEOREMS.md Thm 3: a(x,t) = gamma(t) a0(x))?
CFD passive-temperature solution (no melting, Ste -> 0): s(x) = ln(T_b(0)/T_b(x)) = G(x)/C0 with G the cumulative effective conductance [W/K] and C0 = rho c Q (constant); P = P0 exp(-s) is exactly the theorem's Theorem-1 map in the Ste -> 0 limit.
Per bin i (3 H wide over 0..30 H): g_i = s(x_hi) - s(x_lo) (conductance in units of C0). Histories A and B* (matched bulk T at x_c) have the SAME flow, so a common shape would give g_B,i / g_A,i = constant (any constant, the multiplier gamma).
Reports: ratio spread, total-conductance ratio, regional-heat forecast error when the shape of A is rescaled to B's total (the common-shape forecast), band widening needed per cell, for step and straight geometries, n = 16 and 32."""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); R = lambda n: json.load(open(os.path.join(HERE, "results", n + ".json")))
def s_of(res):
    st, tb = res["Tb_x"]; st = np.array(st); tb = np.array(tb); i0 = int(np.argmin(np.abs(st - 0.0))); return st, np.log(tb[i0] / tb)
def bins(st, s, edges): return np.array([np.interp(b, st, s) - np.interp(a, st, s) for a, b in zip(edges[:-1], edges[1:])])
out = {}
for kind in ("step", "straight"):
    for n in (16, 32):
        d = R(f"{kind}_n{n}_base"); A = d["A_dt10"]; B = d["B_dt10"]; stA, sA = s_of(A); stB, sB = s_of(B); edges = np.arange(0.0, 30.01, 3.0)
        gA = bins(stA, sA, edges); gB = bins(stB, sB, edges); r = gB / gA; rbar = float(np.sum(gB) / np.sum(gA))
        # regional heat (fraction of the 0-section power) in thirds of the 0-30 H reach with the common-shape forecast: shape of A, total of B
        def regional(g, cuts): s = np.concatenate([[0], np.cumsum(g)]); P = np.exp(-s); return [float(P[c0] - P[c1]) for c0, c1 in cuts]
        cuts = [(0, 3), (3, 6), (6, 10)]; Hb = regional(gB, cuts); gAscaled = gA * (np.sum(gB) / np.sum(gA)); Hf = regional(gAscaled, cuts)
        out[f"{kind}_n{n}"] = dict(bin_edges_H=edges.tolist(), g_A=gA.tolist(), g_B=gB.tolist(), ratio_B_over_A=r.tolist(), total_ratio=rbar, ratio_min=float(r.min()), ratio_max=float(r.max()), ratio_rel_spread=float((r.max() - r.min()) / rbar),
                                  band_halfwidth_needed_rel=float(np.max(np.abs(r / rbar - 1))), regional_heat_true_B=Hb, regional_heat_common_shape_forecast=Hf, regional_rel_error=[(f - t) / t for f, t in zip(Hf, Hb)])
json.dump(out, open(os.path.join(HERE, "conductance_test.json"), "w"), indent=1)
for k, v in out.items(): print(k, "ratio B/A per bin", [round(x, 3) for x in v["ratio_B_over_A"]], "total", round(v["total_ratio"], 3), "relspread %.3f" % v["ratio_rel_spread"], "band needed ±%.3f" % v["band_halfwidth_needed_rel"], "regional err", [round(x, 3) for x in v["regional_rel_error"]])
