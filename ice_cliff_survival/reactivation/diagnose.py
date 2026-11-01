"""Fixed-bed, isothermal diagnostics of the partition discrepancy (M3 step 1). Original failed checks are preserved in results/verify_swe_v4.json; this file only adds isolation tests.
`python -m ice_cliff_survival.reactivation.diagnose OUT.json [test ...]`"""
import json, math, sys, time
import numpy as np
from . import cases, verify_swe as V
from .swe2d import Sim, G
def stage(s): return float(s.h[s.i0 - 3, s.nw + 1] + s.zb[s.i0 - 3, s.nw + 1])      # pool stage at the pool/lane junction row
def steady(s, tmax=4000.0): r = s.hydraulics(tmax=tmax, tol=1e-4, cool=False); return r
def single_lane(dx=0.05, pool_floor=-1.5, n=0.02, ztoe=0.0, zbyp=0.0, **kw):
    s = Sim(dx=dx, z_toe=ztoe, z_byp=zbyp, pool_floor=pool_floor, manning=n, **kw); steady(s); d = s.diagnostics(); eta = stage(s); Qt, Qb = d["Q_toe"], d["Q_byp"]
    out = dict(dx=dx, pool_floor=pool_floor, n=n, eta_model=eta, f_toe=d["f_toe"], Q_toe=Qt, Q_byp=Qb)
    # Cd of each lane against K (eta - z)^1.5 (frictionless broad-crested weir) and the GVF reference with friction
    K = 0.5 * (2 / 3) ** 1.5 * math.sqrt(G)
    out["Cd_toe"] = Qt / (K * max(eta - ztoe, 1e-9) ** 1.5); out["Cd_byp"] = Qb / (K * max(eta - zbyp, 1e-9) ** 1.5)
    eta_ref, f = V.gvf_ref(0.07, [(ztoe, 0.5), (zbyp, 0.5)], n=max(n, 1e-6)); out["eta_gvf"] = eta_ref; out["f_gvf"] = float(f[0]); return out
def profile(s, ys):
    j = ys.start + 3; rows = []
    for i in range(s.i0 - 2, s.nx + 1, 2):
        h = s.h[i, j]; u = s.hu[i, j] / h if h > 1e-4 else 0.0; E = s.zb[i, j] + s.solid[i, j] * s.zd[i, j] + h + u * u / (2 * G); rows.append(dict(x=(i - 0.5) * s.dx, eta=float(h + s.zb[i, j] + s.solid[i, j] * s.zd[i, j]), h=float(h), u=float(u), Fr=float(u / math.sqrt(G * h)) if h > 1e-4 else 0.0, E=float(E)))
    return rows
TESTS = {}
def test(f): TESTS[f.__name__] = f; return f
@test
def t1_symmetric_lanes():
    """identical lanes (no junction asymmetry, no deposit): stage and Cd vs weir/GVF, three resolutions, two friction values"""
    return [single_lane(dx=dx, n=n) for dx in (0.1, 0.05, 0.025) for n in (0.02, 0.001)]
@test
def t2_entrance_step_height():
    """does the pool->lane step height (entrance geometry) change the discharge coefficient?"""
    return [single_lane(dx=0.05, pool_floor=pf, n=0.001) for pf in (-0.1, -0.3, -0.6, -1.5)]
@test
def t3_tail_and_brink():
    """tail drop / slope: is the brink control or the entrance the control?"""
    return [single_lane(dx=0.05, n=0.001, tail_drop=td, tail_slope=ts) for td, ts in ((0.0, 0.5), (0.5, 0.5), (1.5, 0.5), (1.5, 1.0))]
@test
def t4_two_floors_no_deposit():
    return [single_lane(dx=0.05, ztoe=-0.08, zbyp=0.0, n=n) for n in (0.02, 0.001)]
@test
def t5_strip_alone():
    """solid strip, frictionless: toe flow vs the weir at the strip crest using the measured pool stage and the measured upstream-of-strip lane stage"""
    out = []
    for dx in (0.1, 0.05, 0.025):
        s = cases.build("strip", dx=dx, manning=0.001); steady(s); d = s.diagnostics(); eta = stage(s); K = 0.5 * (2 / 3) ** 1.5 * math.sqrt(G); pr = profile(s, s.toe_cols())
        mid = [r for r in pr if 1.4 < r["x"] < 1.5][0] if any(1.4 < r["x"] < 1.5 for r in pr) else pr[len(pr) // 2]
        out.append(dict(dx=dx, f_toe=d["f_toe"], eta_pool=eta, Q_toe=d["Q_toe"], Q_weir_pool=K * max(eta - 0.12, 0) ** 1.5, Q_weir_laneH=K * max(mid["E"] - 0.12, 0) ** 1.5, lane_E_before_strip=mid["E"], profile=pr))
    return out
@test
def t6_energy_loss_along_branches():
    s = cases.build("strip", dx=0.025); steady(s); return dict(toe=profile(s, s.toe_cols()), byp=profile(s, s.byp_cols()), eta_pool=stage(s))
if __name__ == "__main__":
    names = sys.argv[2:] or list(TESTS); res = {}
    for n in names:
        t = time.time(); res[n] = TESTS[n](); print(n, "done", round(time.time() - t), "s", flush=True)
        json.dump(res, open(sys.argv[1], "w"), indent=1, default=float)
