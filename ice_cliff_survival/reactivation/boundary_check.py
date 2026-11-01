"""Finest-grid removal-boundary check (fixed bed). `python -m ice_cliff_survival.reactivation.boundary_check ZB_NPZ OUT.json dx1,dx2,... day1,day2,...`
Uses bed snapshots of the medium no-removal strip run (V2M_strip_none_zb.npz, one per 6-h update). For each removal day the deposit is removed from the snapshot bed, the bed is transferred bilinearly to the target grid (cell-centre interpolation of the full padded bed), and the hydraulics are advanced to steady state (warm-started across days; heat sink off) -> toe flow fraction f_toe(t_r+). LIMITATION: the bed HISTORY is that of the medium run; only the post-removal hydraulics are resolved at the finer grid."""
import json, sys, time
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from . import cases
from .swe2d import Sim
def transfer(zb_src, dx_src, sim):
    nxs, nys = zb_src.shape[0] - 2, zb_src.shape[1] - 2; xs = (np.arange(nxs) + 0.5) * dx_src; ys = (np.arange(nys) + 0.5) * dx_src; it = RegularGridInterpolator((xs, ys), zb_src[1:-1, 1:-1], bounds_error=False, fill_value=None, method="linear")
    xt = (np.arange(sim.nx) + 0.5) * sim.dx; yt = (np.arange(sim.ny) + 0.5) * sim.dx; X, Y = np.meshgrid(xt, yt, indexing="ij"); z = it(np.stack([X, Y], -1)); sim.zb[1:-1, 1:-1] = z; sim.zb[:, 0] = sim.zb[:, 1]; sim.zb[:, -1] = sim.zb[:, -2]; sim.zb[0] = sim.zb[1]; sim.zb[-1] = sim.zb[-2]
def main(npz, out, dxs, days, dx_src=0.025, dtm=6 * 3600.0):
    zb = np.load(npz)["zb"]; res = []
    for dx in dxs:
        s = cases.build("bare", dx=dx); first = True; t0 = time.time()
        for d in days:
            i = int(round(d * 86400 / dtm)); transfer(zb[i].astype(float), dx_src, s)
            if first: s.reset_state(max(0.15, float(s.zb.max()) + 0.2)); first = False
            r = s.hydraulics(tmax=3000.0, tol=2e-4, cool=False); dd = s.diagnostics(); res.append(dict(dx=dx, removal_day=d, f_toe=float(dd["f_toe"]), Q_toe=float(dd["Q_toe"]), converged=r["converged"], hyd_time=r["hyd_time"], water_resid=float(dd["water_resid"]), wall_s=time.time() - t0))
            print(res[-1], flush=True); json.dump(res, open(out, "w"), indent=1)
if __name__ == "__main__": main(sys.argv[1], sys.argv[2], [float(x) for x in sys.argv[3].split(",")], [float(x) for x in sys.argv[4].split(",")])
