"""Statistics of the fast discharge solver over a sweep of levels, tailwaters and geometries (reference = unchanged bisection)."""
import json, os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src"))
from ice_network.hydraulics.channel import edge_flow_diag
from ice_network.hydraulics.channel_ref import edge_flow_ref
NM, KE = 0.04, 0.5
rows = []
t_fast = t_ref = 0.0
for (z0, drop, w0, L, n) in [(100.4985, 1.0, 0.5, 100.0, 9), (100.0, 0.02, 0.5, 100.0, 9), (100.0, 0.05, 0.3, 20.0, 5)]:
    z = z0 - drop * np.linspace(0, 1, n); w = np.full(n, w0); dx = L / (n - 1)
    h = np.zeros(n); hf = np.zeros(n - 1); ch = np.zeros(n, np.int64); dg = np.zeros(4)
    for head in np.concatenate([[1e-8, 1e-6, 1e-4], np.geomspace(1e-3, 2.0, 40)]):
        for td in (-5.0, 0.0, 0.3 * head, 0.8 * head, 0.95 * head):
            h_i, h_j = z0 + head, z[-1] + td
            if h_j >= h_i: continue
            t0 = time.perf_counter(); qf = edge_flow_diag(h_i, h_j, z, w, NM, dx, KE, 16, h, hf, ch, dg)[0]; t_fast += time.perf_counter() - t0
            d = dg.copy()
            t0 = time.perf_counter(); qr = edge_flow_ref(h_i, h_j, z, w, NM, dx, KE, 16, h, hf, ch)[0]; t_ref += time.perf_counter() - t0
            rows.append((qf, qr, d))
qf = np.array([r[0] for r in rows]); qr = np.array([r[1] for r in rows]); d = np.array([r[2] for r in rows])
m = qr > 1e-8
out = dict(n_cases=len(rows), n_above_1em8=int(m.sum()), max_rel_diff=float(np.max(np.abs(qf[m] - qr[m]) / qr[m])), median_rel_diff=float(np.median(np.abs(qf[m] - qr[m]) / qr[m])),
           max_abs_diff_below_1em8=float(np.max(np.abs(qf[~m] - qr[~m]))) if (~m).any() else 0.0,
           outlet_control_cases=int((d[:, 3] == 2).sum()), entrance_control_cases=int((d[:, 3] == 1).sum()),
           brent_iterations_mean=float(d[d[:, 3] == 2, 0].mean()), brent_iterations_max=float(d[d[:, 3] == 2, 0].max()),
           max_entrance_residual_m=float(d[d[:, 3] == 2, 1].max()), speedup_per_call=t_ref / t_fast)
json.dump(out, open(os.path.join(HERE, "..", "records", "solver_sweep_stats.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
