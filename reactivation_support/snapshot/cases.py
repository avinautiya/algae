"""Layouts and the quasi-steady morphology runner (M3). Layout names: bare, corridor, strip (solid transverse), porous_strip (Ergun drag; same area/volume), plus the same with prescribed removal."""
import json, math, os, time, resource
import numpy as np
from .swe2d import Sim, RHO_I, LF
def build(layout, dx=0.05, depth=0.2, cover=0.5, Q=0.07, Tw=0.5, **kw):
    if layout == 'porous_strip': kw.setdefault('porous', (1.0, 0.4, 0.2))
    s = Sim(dx=dx, Q=Q, Tw=Tw, **kw); nw = s.nw; i0, i1 = s.i0, s.i1; nxr = i1 - i0 + 1; M = np.zeros((s.nx, s.ny), bool)          # M indexed [i-1, j-1]
    if layout == "bare": return s
    ncor = int(round(cover * nw))            # corridor: deposit covers `cover` of the toe-lane width over the whole reach
    if layout == "corridor": M[i0 - 1:i1, 0:ncor] = True; s.set_deposit(M, depth, "solid")
    elif layout in ("strip", "porous_strip"):
        L = int(round(cover * nxr)); a = i0 - 1 + (nxr - L) // 2; M[a:a + L, 0:nw] = True
        s.set_deposit(M, depth, "solid" if layout == "strip" else "porous")
    else: raise ValueError(layout)
    s.deposit_cells = int(M.sum()); s.reset_state(0.15); return s
def crest(s, lane):
    zs = s.zs(); sl = s.toe_cols() if lane == "toe" else s.byp_cols(); r = zs[s.i0:s.i1 + 1, sl]; return float(r.max(axis=0).min())      # lowest lane-wise controlling crest (a bare corridor lane gives its floor max)
def run(layout, dx=0.05, dtm=6 * 3600.0, days=45.0, removal_day=None, spin=400.0, log=None, hyd_max=60.0, **kw):
    """one development/evaluation run. Hydraulics are warm-started and advanced to quasi-steady (partition drift < 2e-4 Q per 4 s, storage < 1e-3 Q) between bed updates."""
    s = build(layout, dx=dx, **kw); t0w = time.time(); rec = []; cum = dict(toe=0.0, byp=0.0, heat_in=0.0, heat_out=0.0, latent=0.0, storage=0.0); t = 0.0; prev_E = None
    s.hydraulics(tmax=3000.0, tol=2e-4); nsteps = int(round(days * 86400 / dtm)); hyd_steps = 0; hyd_time = 0.0
    for n in range(nsteps + 1):
        if removal_day is not None and not s.removed and t >= removal_day * 86400 - 1e-6:
            s.remove_debris(); s.export_log[-1]["t"] = t; s.hydraulics(tmax=3000.0, tol=2e-4)
        d = s.diagnostics(); h_ = dict(day=t / 86400, **{k: float(v) for k, v in d.items()}); h_.update(crest_toe=crest(s, "toe"), crest_byp=crest(s, "byp"), floor_toe_mean=float(s.zb[s.i0:s.i1 + 1, s.toe_cols()].mean()), floor_toe_min=float(s.zb[s.i0:s.i1 + 1, s.toe_cols()].min()), floor_byp_mean=float(s.zb[s.i0:s.i1 + 1, s.byp_cols()].mean()),
            debris_volume=float((np.where(s.solid, s.zd, 0.0) + s.pdep).sum() * dx ** 2), cum_melt_toe=cum["toe"], cum_melt_byp=cum["byp"], hyd_time=s.hyd_time, hyd_steps=s.hyd_steps, dt_limit=s.cfl * dx / s.smax, wall_s=time.time() - t0w)
        rec.append(h_)
        if n == nsteps: break
        vol, vt, vb = s.lower_bed(dtm); cum["toe"] += vt; cum["byp"] += vb; cum["heat_in"] += d["heat_in_W"] * dtm; cum["heat_out"] += d["heat_out_W"] * dtm; cum["latent"] += RHO_I * LF * vol; cum["storage"] += 0.0
        t += dtm; s.t = t; r = s.hydraulics(tmax=hyd_max, tol=2e-4); hyd_steps += r["steps"]; hyd_time += r["hyd_time"]
    return dict(layout=layout, dx=dx, dtm=dtm, days=days, removal_day=removal_day, records=rec, cum=cum, export=s.export_log, cells=s.nx * s.ny, wall_s=time.time() - t0w, maxrss_MB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, deposit_volume=getattr(s, "deposit_volume", 0.0), hyd_steps=hyd_steps)
def contact_times(records, theta=1e-3, key="f_toe"):
    """first time f_toe falls below theta (linear interpolation of the continuous series) and later reconnection (first rise above theta after removal)"""
    t = np.array([r["day"] for r in records]); f = np.array([r[key] for r in records]); lost = None; rec = None
    for i in range(1, len(f)):
        if lost is None and f[i] < theta <= f[i - 1]: lost = float(t[i - 1] + (t[i] - t[i - 1]) * (f[i - 1] - theta) / (f[i - 1] - f[i]))
        elif lost is not None and f[i] >= theta > f[i - 1]: rec = float(t[i - 1] + (t[i] - t[i - 1]) * (theta - f[i - 1]) / (f[i] - f[i - 1])); break
    return lost, rec
