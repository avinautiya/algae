"""Profile of the archived Dark Zone target over the algal concentration (competitor observation model), on a fixed subset.

For each grid value of algae (cells/mL) the other six unconstrained coordinates are re-optimised (multi-start L-BFGS-B).
Two objectives: 'target' = log posterior density as archived (Student-t nu=3, sigma=0.003, Beta priors, Jacobians) and
'likelihood' = the Student-t term only (no priors, no Jacobians), restricted to the same box (algae <= prior cap).
Interval rule (approximate, the t-likelihood is not chi2): drop of 1.92 nats from the maximum. A side is bound-limited when
the drop is not reached before the grid end (1e2 cells/mL or the prior cap 5e5).
"""
import json, os, sys, numpy as np, h5py
from scipy.optimize import minimize
from scipy.special import logit
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import reproduce_darkzone as RD
AR = os.environ["DARKZONE_ARCHIVE"]
layers = RD.load_mlp(os.path.join(AR, "data/emulator/mlp_icelaps_101025.h5"))
f = h5py.File(os.path.join(AR, "data/spectra/HCRF_ASDFieldspec4_spectra_Chevrollier_et_al_2026.nc"), "r")
obs = f["reflectance"][:][0:851:5, :]; tags = [t.decode() if isinstance(t, bytes) else t for t in f["measurement_tag"][:]]
per = {r["tag"]: r for r in json.load(open(os.path.join(HERE, "..", "records", "inversion_m1", "darkzone_repro", "per_spectrum.json")))}
sel = list(range(0, 105, 9))
grid = np.geomspace(1e2, 4.9e5, 11)
rng = np.random.default_rng(7)

def ll_only(z, y):
    th, _ = RD.z_to_theta(z); pred = RD.mlp(layers, th)[0, 10:181]; r = (pred - y) / RD.SIGMA
    return -0.5 * (RD.NU + 1) * np.sum(np.log1p(r * r / RD.NU))

def prof(y, with_prior, z_best):
    out = []
    for g in grid:
        u = (g - RD.DATA_MIN[6]) / RD.RANGE[6]; zg = float(logit(u)); best = -np.inf
        def obj(zr):
            z = np.concatenate([zr, [zg]])
            v = RD.log_target(z, layers, y)[0] if with_prior else ll_only(z, y)
            return -v if np.isfinite(v) else 1e12
        starts = [np.delete(z_best, 6)] + [np.delete(RD.init_z(rng), 6) for _ in range(3)]
        for s in starts:
            r = minimize(obj, s, method="L-BFGS-B", options={"maxiter": 200})
            best = max(best, -r.fun)
        out.append(best)
    return np.array(out)

res = []
for j in sel:
    y = obs[:, j]; p = per[tags[j]]
    # best z from the stored MAP fit is not saved; re-fit quickly (6 starts)
    (lt, zb, _), _ = RD.fit_one(layers, y, 6, rng)
    row = dict(tag=tags[j], map_algae=float(RD.DATA_MIN[6] + RD.z_to_theta(zb)[0][6] * RD.RANGE[6]))
    for name, wp in (("target", True), ("likelihood", False)):
        v = prof(y, wp, zb); d = v.max() - v
        ins = np.flatnonzero(d <= 1.92)
        row[name] = dict(drop_nats=[float(x) for x in d], lower_bound_limited=bool(ins[0] == 0), upper_bound_limited=bool(ins[-1] == len(grid) - 1),
                         interval_dex_width=float(np.log10(grid[ins[-1]] / grid[ins[0]])))
    res.append(row); print(row["tag"], {k: (row[k]["lower_bound_limited"], row[k]["upper_bound_limited"], round(row[k]["interval_dex_width"], 2)) for k in ("target", "likelihood")}, flush=True)
json.dump(dict(grid=grid.tolist(), rows=res), open(os.path.join(HERE, "..", "records", "inversion_m1", "darkzone_repro", "profile_algae.json"), "w"), indent=1)
