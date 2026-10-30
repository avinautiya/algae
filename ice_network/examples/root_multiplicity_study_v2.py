"""Bounded, section-wise version (v2): each section writes its own record in records/root_study/ and prints progress; timing variants run one per process under a wall limit.
Investigation of multiple roots of the entrance equation R(Q) = H_req(Q) - H_a at small discharge.

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


OUTDIR = os.path.join(HERE, "..", "records", "root_study"); os.makedirs(OUTDIR, exist_ok=True)


def save(name, obj):
    p = os.path.join(OUTDIR, name + ".json"); json.dump(obj, open(p + ".tmp", "w"), indent=1); os.replace(p + ".tmp", p)


def section_scan():
    out = {"scan": []}
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
    save("scan", out)
    multi = [r for r in out["scan"] if r["n_roots"] > 1]
    print("scan cases", len(out["scan"]), "with >1 root:", len(multi), flush=True)


def section_ramp():
    out = {}
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
    save("ramp", out["ramp"]); print("ramp", {k: v for k, v in out["ramp"].items() if k != "heads_where_roots_differ_gt_1em6_rel"}, flush=True)


# ---------------------------------------------------------------- 4 event timing (slowly filling basin crossing a sill)
EVAL_CAP = int(os.environ.get("EVAL_CAP", 150000)); RTOL = float(os.environ.get("RTOL", 1e-8))
class Stalled(Exception): pass

def counting(fn, limit):
    cnt = [0]
    def wrapped(*args):
        cnt[0] += 1
        if cnt[0] > limit: raise Stalled
        return fn(*args)
    return wrapped, cnt

def activation(fn, inflow, h_floor, area=500.0, limit=EVAL_CAP):
    wrapped, cnt = counting(fn, limit)
    NW.edge_flow = wrapped
    e = chan(0, 1, 100.0, 0.02, 0.5, L=100.0, n=9)
    net = NW.Network([lake(99.0, area)], [dict(h_fixed=95.0)], [e], melting=False, h_floor=h_floor)
    t_cross = area * 0.02 / inflow
    t_end = 3 * t_cross
    t0 = time.perf_counter()
    try:
        r = net.run([area * 0.98], [0.0], t_end, t_end / 1500.0, const(inflow), const(0.0), evolve=False, rtol=RTOL, atol=1e-10)
    except Stalled:
        return dict(stalled=True, edge_flow_calls=cnt[0], t_cross_nominal=t_cross), None, None
    t = np.array(r["t"]); Q = np.array([abs(s["Q"][0]) for s in r["hist"]]); st = np.array([s["stage"][0] for s in r["hist"]])
    res = dict(stalled=False, edge_flow_calls=cnt[0], seconds=time.perf_counter() - t0)
    for thr in (1e-9, 1e-6, 1e-4, 1e-3 * inflow, 1e-2 * inflow):
        k = np.argmax(Q > thr); res[f"t_Q_gt_{thr:.3g}"] = float(t[k]) if Q.max() > thr else None
    res["t_stage_reaches_sill"] = float(t[np.argmax(st >= 100.0)])
    res["delivered_m3"] = float(r["budgets"]["delivered_m3"][0]); res["Q_end"] = float(Q[-1]); res["t_cross_nominal"] = t_cross
    return res, t, Q


VARIANTS = {"raw_fast": (C.edge_flow, 0.0), "raw_reference_bisection": (CR.edge_flow_ref, 0.0), "floor_1e-4": (C.edge_flow, 1e-4),
            "floor_1e-3": (C.edge_flow, 1e-3), "floor_1e-2": (C.edge_flow, 1e-2)}


def section_timing(inflow, name):
    fn, hf_ = VARIANTS[name]; t0 = time.perf_counter()
    res, t, Q = activation(fn, float(inflow), hf_)
    res["wall_s"] = time.perf_counter() - t0; res["eval_cap"] = EVAL_CAP; res["rtol"] = RTOL
    if t is not None: np.savez_compressed(os.path.join(OUTDIR, f"timing_{inflow}_{name}.npz"), t=t, Q=Q)
    save(f"timing_{inflow}_{name}", res); print(inflow, name, {k: v for k, v in res.items() if k.startswith("t_") or k in ("stalled", "edge_flow_calls", "wall_s")}, flush=True)


if __name__ == "__main__":
    sec = sys.argv[1]
    if sec == "scan": section_scan()
    elif sec == "ramp": section_ramp()
    elif sec == "timing": section_timing(sys.argv[2], sys.argv[3])
