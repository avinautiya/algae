"""Profile of the reduced solver's discharge-law cost: per-call Brent iterations and wall time, by stations. Usage: profile_reduced.py SCN STATIONS DT_GEOM"""
import sys, time, json, os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numpy as np, run_eval as R
import ice_network.evolution.network as NW
from ice_network.hydraulics.channel import edge_flow_diag
scn = next(s for s in R.CFG["scenarios"] if s["id"] == sys.argv[1]); st = int(sys.argv[2]); dt = float(sys.argv[3])
rows = []; dg = np.zeros(4)
def wrapped(h_i, h_j, z, w, n_man, dx, ke, msub, h, hf, ch):
    t0 = time.perf_counter(); o = edge_flow_diag(h_i, h_j, z, w, n_man, dx, ke, msub, h, hf, ch, dg); dt_ = time.perf_counter() - t0
    rows.append((dt_, dg[0], dg[3], abs(o[0]), max(h_i, h_j) - max(z[0], z[-1]) if True else 0)); return o
NW.edge_flow = wrapped
orig = R.Network.run
def run(self, *a, **k): k["sample_dt"] = None; return orig(self, *a, **k)
R.Network.run = run
t0 = time.perf_counter(); R.run_reduced(scn, st, dt, sample_dt=None); wall = time.perf_counter() - t0
a = np.array(rows); tot = a[:, 0].sum()
out = dict(scenario=sys.argv[1], stations=st, dt_geom=dt, wall_two_runs_s=wall, calls=len(a) / 2, time_in_edge_flow_s=tot, frac_time_in_edge_flow=tot / wall,
           mean_call_ms=1e3 * a[:, 0].mean(), brent_iter_mean=a[:, 1].mean(), brent_iter_max=a[:, 1].max(), frac_outlet_control=float((a[:, 2] == 2).mean()),
           frac_entrance_control=float((a[:, 2] == 1).mean()), frac_zero_Q=float((a[:, 3] == 0).mean()), ms_per_call_by_control={str(int(c)): 1e3 * a[a[:, 2] == c, 0].mean() for c in np.unique(a[:, 2])},
           ms_per_brent_iter=1e3 * tot / max(a[:, 1].sum(), 1))
print(json.dumps(out)); os.makedirs("records/profile", exist_ok=True); json.dump(out, open(f"records/profile/{sys.argv[1]}_r{st}.json", "w"), indent=1)
