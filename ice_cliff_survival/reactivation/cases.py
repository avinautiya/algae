"""Layouts and the quasi-steady morphology runner (M3). Layout names: bare, corridor, strip (solid transverse), porous_strip (Ergun drag; same area/volume), plus the same with prescribed removal."""
import json, math, os, time, resource
import numpy as np
from .swe2d import Sim, RHO_I, LF
V_TARGET = 0.1      # equal solid volume [m3] for every deposit layout (= 0.2 m x 0.5 m x 1.0 m rectangle); equal footprint 0.5 m2
def deposit_shape(layout, s, taper=0.3):
    """unit-height shape (nx, ny): trapezoidal (side slope taper/height); strip: spans the lane width, centred, footprint 1.0 m x 0.5 m; corridor: full reach length (tapered at both ends), width 0.25 m at the lane side wall, tapered on the free edge."""
    nw, dx = s.nw, s.dx; i0, i1 = s.i0, s.i1; nxr = i1 - i0 + 1; H = np.zeros((s.nx, s.ny)); ramp = lambda d: np.clip(d / taper, 0.0, 1.0)
    if layout == "corridor":
        ncor = int(round(0.5 * nw)); xs = (np.arange(nxr) + 0.5) * dx; fx = np.minimum(ramp(xs), ramp(nxr * dx - xs)); ys = (np.arange(ncor) + 0.5) * dx; fy = ramp(ncor * dx - ys)
        H[i0 - 1:i1, 0:ncor] = np.outer(fx, fy)
    else:
        L = int(round(0.5 * nxr)); a = i0 - 1 + (nxr - L) // 2; xs = (np.arange(L) + 0.5) * dx; fx = np.minimum(ramp(xs), ramp(L * dx - xs)); H[a:a + L, 0:nw] = fx[:, None]
    return H
def set_regions(s, fp):
    """region masks (full padded arrays): toe_fp = deposit footprint cells in the toe reach (identical footprint to the reduced model's covered control), toe_open = other toe-reach cells, byp = bypass reach, ramp_toe/ramp_byp = entrance ramps (ice, not in any reduced control). For 'bare' toe_fp is empty."""
    rm = np.zeros_like(s.ice); rm[s.i0:s.i1 + 1] = True; F = np.zeros_like(s.ice)
    if fp is not None: F[1:-1, 1:-1] = fp
    s.regions = dict(toe_fp=s.toe & rm & F, toe_open=s.toe & rm & ~F, byp=s.byp & rm, ramp_toe=s.toe & ~rm, ramp_byp=s.byp & ~rm); s.region_area = {k: float(m.sum() * s.dx ** 2) for k, m in s.regions.items()}
def build(layout, dx=0.05, cover=0.5, Q=0.07, Tw=0.5, taper=0.1, **kw):
    if layout == 'porous_strip': kw.setdefault('porous', (1.0, 0.4, 0.2))
    s = Sim(dx=dx, Q=Q, Tw=Tw, **kw)
    if layout == "bare": set_regions(s, None); return s
    if layout not in ("corridor", "strip", "porous_strip"): raise ValueError(layout)
    sh = deposit_shape("corridor" if layout == "corridor" else "strip", s, taper); scale = V_TARGET / (sh.sum() * dx ** 2)
    s.set_deposit(sh * scale, "solid" if layout != "porous_strip" else "porous"); set_regions(s, sh > 0); s.deposit_cells = int((sh > 0).sum()); s.deposit_scale = float(scale); s.reset_state(0.15); return s
def crest(s, lane):
    zs = s.zs(); sl = s.toe_cols() if lane == "toe" else s.byp_cols(); r = zs[s.i0:s.i1 + 1, sl]; return float(r.max(axis=0).min())      # lowest lane-wise controlling crest (a bare corridor lane gives its floor max)
def run(layout, dx=0.05, dtm=6 * 3600.0, days=45.0, removal_day=None, spin=400.0, log=None, hyd_max=60.0, theta_dump=None, **kw):
    """one development/evaluation run. Hydraulics are warm-started and advanced to quasi-steady (partition drift < 2e-4 Q per 4 s, storage < 1e-3 Q) between bed updates."""
    s = build(layout, dx=dx, **kw); t0w = time.time(); rec = []; snaps = []; cum = dict(toe=0.0, byp=0.0, toe_reach=0.0, byp_reach=0.0, heat_in=0.0, heat_out=0.0, latent=0.0, storage=0.0); t = 0.0; prev_E = None
    s.hydraulics(tmax=3000.0, tol=2e-4); nsteps = int(round(days * 86400 / dtm)); hyd_steps = 0; hyd_time = 0.0
    for n in range(nsteps + 1):
        if removal_day is not None and not s.removed and t >= removal_day * 86400 - 1e-6:
            s.remove_debris(); s.export_log[-1]["t"] = t; s.hydraulics(tmax=3000.0, tol=2e-4)
        d = s.diagnostics(); h_ = dict(day=t / 86400, **{k: float(v) for k, v in d.items()}); h_.update(crest_toe=crest(s, "toe"), crest_byp=crest(s, "byp"), floor_toe_mean=float(s.zb[s.i0:s.i1 + 1, s.toe_cols()].mean()), floor_toe_min=float(s.zb[s.i0:s.i1 + 1, s.toe_cols()].min()), floor_byp_mean=float(s.zb[s.i0:s.i1 + 1, s.byp_cols()].mean()),
            region_heat_W=s.region_heat(), region_cum_melt_m3=dict(s.region_vol), debris_volume=float((np.where(s.solid, s.zd, 0.0) + s.pdep).sum() * dx ** 2), cum_melt_toe=cum["toe"], cum_melt_byp=cum["byp"], cum_melt_toe_reach=cum["toe_reach"], cum_melt_byp_reach=cum["byp_reach"], cum_latent_J=cum["latent"], f_toe_cum_water_m3=cum.get("water_toe", 0.0), hyd_time=s.hyd_time, hyd_steps=s.hyd_steps, dt_limit=s.cfl * dx / s.smax, wall_s=time.time() - t0w, toe_floor_x=[round(float(v), 5) for v in s.zb[s.i0:s.i1 + 1, s.toe_cols()].mean(axis=1)], byp_floor_x=[round(float(v), 5) for v in s.zb[s.i0:s.i1 + 1, s.byp_cols()].mean(axis=1)], toe_crest_x=[round(float(v), 5) for v in s.zs()[s.i0:s.i1 + 1, s.toe_cols()].max(axis=1)], heat_in=float(d['heat_in_W']), heat_out=float(d['heat_out_W']), storage_W=float(d['storage_rate_W']))
        rec.append(h_); snaps.append(s.zb.astype(np.float32).copy())
        if n == nsteps: break
        vol, vt, vb = s.lower_bed(dtm); cum["toe"] += vt; cum["byp"] += vb; cum["toe_reach"] += s.last_reach[0]; cum["byp_reach"] += s.last_reach[1]; cum["water_toe"] = cum.get("water_toe", 0.0) + d["Q_toe"] * dtm; cum["heat_in"] += d["heat_in_W"] * dtm; cum["heat_out"] += d["heat_out_W"] * dtm; cum["latent"] += RHO_I * LF * vol; cum["storage"] += 0.0
        t += dtm; s.t = t; r = s.hydraulics(tmax=hyd_max, tol=2e-4); hyd_steps += r["steps"]; hyd_time += r["hyd_time"]
    return dict(_snaps=np.array(snaps), region_area=s.region_area, incision=s.incision, layout=layout, dx=dx, dtm=dtm, days=days, removal_day=removal_day, records=rec, cum=cum, export=s.export_log, cells=s.nx * s.ny, wall_s=time.time() - t0w, maxrss_MB=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, deposit_volume=getattr(s, "deposit_volume", 0.0), hyd_steps=hyd_steps)
def contact_times(records, theta=1e-3, key="f_toe"):
    """first time f_toe falls below theta (linear interpolation of the continuous series) and later reconnection (first rise above theta after removal)"""
    t = np.array([r["day"] for r in records]); f = np.array([r[key] for r in records]); lost = None; rec = None
    for i in range(1, len(f)):
        if lost is None and f[i] < theta <= f[i - 1]: lost = float(t[i - 1] + (t[i] - t[i - 1]) * (f[i - 1] - theta) / (f[i - 1] - f[i]))
        elif lost is not None and f[i] >= theta > f[i - 1]: rec = float(t[i - 1] + (t[i] - t[i - 1]) * (theta - f[i - 1]) / (f[i] - f[i - 1])); break
    return lost, rec
