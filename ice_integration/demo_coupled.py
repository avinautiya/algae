"""Coupled demonstration run: plan-view swe2d (solver v4) driven by the field-forcing interface. `python -m ice_integration.demo_coupled OUT.json dx n_intervals [strict|reconstructed] [interval_h]`.
Constant Q is an IMPOSED SCENARIO (no discharge measurement is used); no observed target exists for this plan-view geometry: timing/resource and behaviour demonstration only. Atomic checkpoints every interval; resumes from OUT.ckpt.npz when present."""
import datetime as dt, json, os, resource, sys, time
from . import forcing, coupling
from ice_cliff_survival.reactivation import cases
def main(out, dx, n, mode="reconstructed", hours=6.0, start=dt.datetime(2019, 7, 20), Q=0.07):
    f = forcing.FieldForcing("P5"); sim = cases.build("strip", dx=dx); prov = dict(forcing=f.provenance(), window_start=start.isoformat(), interval_h=hours, Q_m3_s=Q, Q_basis="imposed_scenario", mode=mode, layout="strip", dx=dx, solver="reactivation.swe2d v4")
    drv = coupling.CoupledDriver(sim, prov, include_reconstructed=(mode == "reconstructed"), hyd_max=60.0); ck = out + ".ckpt.npz"; done = 0; t0 = time.time(); rows = []
    if os.path.exists(ck): meta = drv.restore(ck); done = meta["n_log"]; rows = json.load(open(out + ".rows.json")) if os.path.exists(out + ".rows.json") else []
    if done == 0: sim.hydraulics(tmax=3000.0, cool=False)
    for k in range(done, n):
        a = f.integrate(start + dt.timedelta(hours=hours * k), start + dt.timedelta(hours=hours * (k + 1))); t1 = time.time(); d = drv.advance(a, Q)
        rows.append(dict(k=k, observed_s=a.observed_s, reconstructed_s=a.reconstructed_s, unsupported_s=a.unsupported_s, Tw=None if d is None else sim.Tw, f_toe=None if d is None else d["f_toe"], wall_s=time.time() - t1, water_resid=None if d is None else d["water_resid"], energy_resid=None if d is None else d["energy_resid_rel"]))
        drv.checkpoint(ck); json.dump(rows, open(out + ".rows.json.tmp", "w")); os.replace(out + ".rows.json.tmp", out + ".rows.json")
    res = dict(provenance=prov, time=drv.time, melt_m3=drv.melt, rows=rows, wall_s=time.time() - t0, maxrss_MB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, checkpoint_bytes=os.path.getsize(ck), cells=sim.nx * sim.ny)
    json.dump(res, open(out + ".tmp", "w"), indent=1, default=float); os.replace(out + ".tmp", out); return res
if __name__ == "__main__":
    r = main(sys.argv[1], float(sys.argv[2]), int(sys.argv[3]), sys.argv[4] if len(sys.argv) > 4 else "reconstructed", float(sys.argv[5]) if len(sys.argv) > 5 else 6.0)
    print(json.dumps({k: r[k] for k in ("time", "melt_m3", "wall_s", "maxrss_MB", "checkpoint_bytes", "cells")}, default=float))
