"""Investigation of multiple roots of the entrance equation R(Q) = H_req(Q) - H_a at small discharge.

1. scan R(Q) on a fine logarithmic Q grid for several geometries/heads and classify sign changes (upward = stable branch);
2. check physical admissibility of each root (depths, Froude, control points) and whether the structure depends on the substep count;
3. branch continuity along a slow head ramp for the production (largest-root) and reference (bisection) solvers;
4. event timing: dry-to-wet activation of a slowly filling basin with both solvers, several activation thresholds.
"""
import json, os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src")); sys.path.insert(0, os.path.join(HERE, "..", "tests"))
import ice_network.hydraulics.channel as C
import ice_network.hydraulics.channel_ref as CR
import ice_network.evolution.network as NW
from helpers import chan, const, lake

G = 9.81; NM, KE = 0.04, 0.5
geoms = {"mild_100m": (100.4985, 1.0, 0.5, 100.0, 9), "flat_100m": (100.0, 0.02, 0.5, 100.0, 9), "short_20m": (100.0, 0.05, 0.3, 20.0, 5)}


def R_scan(z, w, dx, h_a, h_j, msub, qgrid):
    n = z.size; h = np.zeros(n); hf = np.zeros(n - 1); ch = np.zeros(n, np.int64)
    return np.array([C._req(q, h_j, z, w, NM, dx, KE, msub, h, hf, ch) - h_a for q in qgrid]), h, ch


out = {"scan": [], "ramp": {}, "timing": {}}
# ---------------------------------------------------------------- 1-2 root structure
for name, (z0, drop, w0, L, n) in geoms.items():
    z = z0 - drop * np.linspace(0, 1, n); w = np.full(n, w0); dx = L / (n - 1)
    for head in (1e-9, 1e-7, 1e-5, 1e-4, 1e-3, 1e-2, 1e-1):
        h_a = z0 + head; h_j = z[-1] - 5.0
        q_in = w0 * np.sqrt(G) * (head / (1.5 + 0.5 * KE)) ** 1.5
        qgrid = q_in * np.geomspace(1e-14, 1.0, 600)
        for msub in (4, 16, 64):
            R, h, ch = R_scan(z, w, dx, h_a, h_j, msub, qgrid)
            s = np.sign(R); idx = np.flatnonzero(s[:-1] * s[1:] < 0)
            roots = []
            for i in idx:
                lo, hi = qgrid[i], qgrid[i + 1]
                for _ in range(80):
                    mid = 0.5 * (lo + hi)
                    r_mid = R_scan(z, w, dx, h_a, h_j, msub, [mid])[0][0]
                    if np.sign(r_mid) == np.sign(R[i]): lo = mid
                    else: hi = mid
                q_root = 0.5 * (lo + hi); Rr, hh, cc = R_scan(z, w, dx, h_a, h_j, msub, [q_root])
                roots.append(dict(Q=float(q_root), direction="up" if R[i] < 0 else "down", Q_over_Qin=float(q_root / q_in), h_min_m=float(hh.min()), h_max_m=float(hh.max()),
                                  critical_control_stations=int(cc.sum()), froude_max=float((q_root / (w0 * hh) / np.sqrt(G * hh)).max())))
            out["scan"].append(dict(geometry=name, head_m=head, msub=msub, Q_in=float(q_in), R_at_tiny_Q=float(R[0]), n_roots=len(roots), roots=roots))
# ---------------------------------------------------------------- 3 branch continuity along a head ramp
z0, drop, w0, L, n = geoms["flat_100m"]
z = z0 - drop * np.linspace(0, 1, n); w = np.full(n, w0); dx = L / (n - 1)
heads = np.geomspace(1e-7, 1e-1, 400)
h = np.zeros(n); hf = np.zeros(n - 1); ch = np.zeros(n, np.int64); dg = np.zeros(4)
qf = np.array([C.edge_flow_diag(z0 + x, z[-1] - 5.0, z, w, NM, dx, KE, 16, h, hf, ch, dg)[0] for x in heads])
qr = np.array([CR.edge_flow_ref(z0 + x, z[-1] - 5.0, z, w, NM, dx, KE, 16, h, hf, ch)[0] for x in heads])
def max_jump(q):
    m = q > 1e-14
    r = np.abs(np.diff(q[m])) / q[m][:-1]
    return float(r.max()) if r.size else 0.0
out["ramp"] = dict(geometry="flat_100m", n_heads=len(heads), production_max_relative_step=max_jump(qf), reference_max_relative_step=max_jump(qr),
                   production_first_nonzero_head=float(heads[np.argmax(qf > 0)]), reference_first_nonzero_head=float(heads[np.argmax(qr > 0)]),
                   max_abs_diff=float(np.max(np.abs(qf - qr))), heads_where_roots_differ_gt_1em6_rel=[float(x) for x, a, b in zip(heads, qf, qr) if b > 1e-12 and abs(a - b) / b > 1e-6][:20])
# ---------------------------------------------------------------- 4 event timing (slowly filling basin crossing a sill)
def activation(fn, inflow, area=500.0):
    NW.edge_flow = fn
    e = chan(0, 1, 100.0, 0.02, 0.5, L=100.0, n=9)
    net = NW.Network([lake(99.0, area)], [dict(h_fixed=95.0)], [e], melting=False)
    t_cross = area * 0.02 / inflow
    t_end = 3 * t_cross
    r = net.run([area * 0.98], [0.0], t_end, t_end / 1500.0, const(inflow), const(0.0), evolve=False, rtol=1e-12, atol=1e-14)
    t = np.array(r["t"]); Q = np.array([abs(s["Q"][0]) for s in r["hist"]]); st = np.array([s["stage"][0] for s in r["hist"]])
    res = {}
    for thr in (1e-12, 1e-9, 1e-6, 1e-4, 1e-3):
        k = np.argmax(Q > thr); res[f"t_Q_gt_{thr:g}"] = float(t[k]) if Q.max() > thr else None
    res["t_stage_reaches_sill"] = float(t[np.argmax(st >= 100.0)])
    res["delivered_m3"] = float(r["budgets"]["delivered_m3"][0]); res["Q_end"] = float(Q[-1])
    return res, t, Q
tim = {}
for inflow in (5e-2, 5e-3, 5e-4):
    a, ta, Qa = activation(C.edge_flow, inflow); b, tb, Qb = activation(CR.edge_flow_ref, inflow)
    tim[f"inflow_{inflow:g}"] = dict(production=a, reference=b,
                                    max_abs_Q_diff=float(np.max(np.abs(Qa - Qb))), max_rel_Q_diff_where_Q_gt_1em6=float(np.max(np.abs(Qa - Qb)[Qb > 1e-6] / Qb[Qb > 1e-6])) if (Qb > 1e-6).any() else None,
                    delivered_rel_diff=abs(a["delivered_m3"] - b["delivered_m3"]) / b["delivered_m3"])
out["timing"] = tim
json.dump(out, open(os.path.join(HERE, "..", "records", "root_multiplicity_study.json"), "w"), indent=1)
# console summary
multi = [r for r in out["scan"] if r["n_roots"] > 1]
print("scan cases", len(out["scan"]), "with >1 root:", len(multi))
for r in out["scan"]:
    if r["msub"] == 16 and r["n_roots"] > 0:
        print(r["geometry"], "head", r["head_m"], "roots", [(x["direction"], "%.3g" % x["Q"], "hmin %.2e" % x["h_min_m"]) for x in r["roots"]])
print("ramp", {k: v for k, v in out["ramp"].items() if k != "heads_where_roots_differ_gt_1em6_rel"})
print(json.dumps(out["timing"], indent=1))
