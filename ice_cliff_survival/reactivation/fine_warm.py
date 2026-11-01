"""Finest-grid boundary point with a WARM START (the cold-start dx=0.0125 settle needs >2 h per removal day). `python -m ice_cliff_survival.reactivation.fine_warm ZB_NPZ OUT.json removal_day`
1) settle the dx=0.025 steady state for the removal-day bed; 2) transfer bed AND free surface / unit discharges bilinearly to dx=0.0125; 3) march short chunks until f_toe drift < 1e-4 and |storage| < 2e-4 Q. Records f_toe(t_r+) at dx = 0.025 and 0.0125 plus the convergence history."""
import json, sys, time
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from . import cases
from .boundary_check import transfer, settle
def interp_field(F, dx_src, sim):
    ns = (F.shape[0] - 2, F.shape[1] - 2); xs = (np.arange(ns[0]) + 0.5) * dx_src; ys = (np.arange(ns[1]) + 0.5) * dx_src; it = RegularGridInterpolator((xs, ys), F[1:-1, 1:-1], bounds_error=False, fill_value=None)
    X, Y = np.meshgrid((np.arange(sim.nx) + 0.5) * sim.dx, (np.arange(sim.ny) + 0.5) * sim.dx, indexing="ij"); return it(np.stack([X, Y], -1))
def main(npz, out, day, dx_c=0.025, dx_f=0.0125, dtm=6 * 3600.0):
    zb = np.load(npz)["zb"][int(round(day * 86400 / dtm))].astype(float); hist = dict(removal_day=day); t0 = time.time()
    c = cases.build("bare", dx=dx_c); transfer(zb, 0.025, c); c.reset_state(-0.13); rc, dc = settle(c); hist["coarse"] = dict(dx=dx_c, f_toe=float(dc["f_toe"]), converged=bool(rc["converged"]), hyd_time=rc["hyd_time"], wall_s=time.time() - t0); print(hist["coarse"], flush=True)
    f = cases.build("bare", dx=dx_f); transfer(zb, 0.025, f); eta = interp_field(c.h + c.zb, dx_c, f); hu = interp_field(c.hu, dx_c, f); hv = interp_field(c.hv, dx_c, f)
    f.h[1:-1, 1:-1] = np.maximum(eta - f.zb[1:-1, 1:-1], 0.0) * (~f.wall[1:-1, 1:-1]); wet = f.h[1:-1, 1:-1] > 1e-4; f.hu[1:-1, 1:-1] = hu * wet; f.hv[1:-1, 1:-1] = hv * wet; f.hT = f.h * f.Tw
    steps = []; last = None
    for k in range(6):
        r = f.hydraulics(tmax=30.0, tol=0, cool=False); d = f.diagnostics(); steps.append(dict(chunk=k, hyd_time=30.0 * (k + 1), f_toe=float(d["f_toe"]), storage=float(d["water_storage_rate"]), wall_s=time.time() - t0)); print(steps[-1], flush=True)
        hist["fine_steps"] = steps; json.dump(hist, open(out, "w"), indent=1)
        if last is not None and abs(d["f_toe"] - last) < 1e-4 and abs(d["water_storage_rate"]) < 2e-4 * f.Q: break
        last = d["f_toe"]
    hist["fine"] = dict(dx=dx_f, f_toe=steps[-1]["f_toe"], converged=bool(len(steps) < 6), wall_s=time.time() - t0); json.dump(hist, open(out, "w"), indent=1); print(hist["coarse"], hist["fine"])
if __name__ == "__main__": main(sys.argv[1], sys.argv[2], float(sys.argv[3]))
