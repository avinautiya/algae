"""Baseline forward prediction of thermal incision in a clear-water, temperate-ice open reach (v2).   python -m ice_network_v2.predict case.json outdir [--restart ckpt.npz] [--stop-at SECONDS]
Model (baseline, NO new closure): 1-D shallow-water finite volume (HLL, hydrostatic reconstruction, semi-implicit Manning friction), upwind temperature transport, Dittus-Boelter wall heat transfer
(times nu_factor), the SAME interface power used for water cooling and ice melt, melt water added where produced, geometry updated every step (see reach.py, derived from ice_network_reference/solver.py).
Depth-averaged: it does not resolve vertical structure, free-surface/jump interior, or boundary-layer development. Not CFD.
Inputs: only initial geometry/water state, forcing time series and constants. No later geometry is read. Discharge at the inlet is FORCING (prescribed), not a prediction.
Outputs (outdir): timeseries.csv, snapshots.npz, summary.json (ledgers, flags, provenance), checkpoint.npz (atomic). Width is output only if evolve_width is true.
Case JSON keys: id, length_m, n_cells, bed{z_inlet_m, slope | z_m[]}, width_m (scalar|list), bank_height_m, n_manning, ice{condition:"temperate"}, inflow{Q:[[t,Q]..], T_C:[[t,T]..], T_reference:"flux-weighted bulk temperature"},
 outlet{type:"free_outfall"|"tailwater", stage_m}, initial{depth_m}, spinup_s, t_end_s, record_every_s, snapshot_every_s, checkpoint_every_s, nu_factor, heat_dissipation, evolve_width, cfl."""
import hashlib, json, os, subprocess, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, ".."))
from ice_network_v2.reach import channel_step_q, G, RHO_W, RHO_I, CP_W, LF, NU_W
VERSION = "v2-baseline-0.1"


def interp(tab, t):
    a = np.asarray(tab, float); return float(np.interp(t, a[:, 0], a[:, 1]))


def build(case):
    n = int(case["n_cells"]); L = float(case["length_m"]); dx = L / n; x = (np.arange(n) + 0.5) * dx; b = case["bed"]
    z = np.array(b["z_m"], float) if "z_m" in b else b["z_inlet_m"] - b["slope"] * x
    w = np.full(n, float(case["width_m"])) if np.isscalar(case["width_m"]) else np.array(case["width_m"], float)
    zs = z + float(case.get("bank_height_m", 1.0))
    return dict(n=n, L=L, dx=dx, x=x, z=z.copy(), w=w.copy(), zs=zs, z0=z.copy(), w0=w.copy())


def refuse(case):
    if case.get("ice", {}).get("condition") != "temperate": raise SystemExit("REFUSED: only ice.condition == 'temperate' is supported (cold-ice conduction is not implemented)")
    if case["outlet"]["type"] not in ("free_outfall", "tailwater"): raise SystemExit("REFUSED: outlet.type must be free_outfall or tailwater")
    for k in ("inflow", "bed", "width_m", "t_end_s", "n_cells", "length_m"):
        if k not in case: raise SystemExit(f"REFUSED: missing '{k}'")
    if "T_reference" not in case["inflow"]: raise SystemExit("REFUSED: inflow.T_reference must state the observation definition of the inlet temperature (e.g. 'flux-weighted bulk')")
    if any(k in case for k in ("later_geometry", "observed_geometry", "z_final")): raise SystemExit("REFUSED: case contains later-geometry fields; the predictor does not consume future geometry")


class Run:
    def __init__(self, case):
        self.case = case; g = build(case); self.g = g; n = g["n"]; d0 = float(case.get("initial", {}).get("depth_m", 0.05))
        self.h = np.full(n, d0); self.Q = np.full(n, interp(case["inflow"]["Q"], 0.0)); self.T = np.full(n, interp(case["inflow"]["T_C"], 0.0)); self.t = 0.0
        self.tot = dict(V_in=0.0, V_out=0.0, sens_in=0.0, sens_out=0.0, heat_ice=0.0, melt_w=0.0, diss=0.0, melted_ice=0.0)
        self.flags = dict(min_depth_m=1e9, max_Fr=0.0, min_Fr=1e9, min_Re_Dh=1e12, max_Re_Dh=0.0, dry_cell_seen=False, negative_inlet_T=False)
        self.rec = []; self.next_rec = 0.0; self.nsteps = 0; self.E0 = None; self.V0 = None

    def state_E(self): return RHO_W * CP_W * float(np.sum(self.g["w"] * self.h * self.T) * self.g["dx"])
    def state_V(self): return float(np.sum(self.g["w"] * self.h) * self.g["dx"])

    def step(self, melting, tend):
        c = self.case; g = self.g; out = np.zeros(8)
        qin = interp(c["inflow"]["Q"], self.t); tin = interp(c["inflow"]["T_C"], self.t)
        wet = self.h > 1e-6; u = np.abs(self.Q / (g["w"] * np.maximum(self.h, 1e-6))); cmax = float(np.max(np.where(wet, u + np.sqrt(G * self.h), 0.0)))
        cmax = max(cmax, qin / (g["w"][0] * max(self.h[0], 1e-3)) + np.sqrt(G * max(self.h[0], 1e-3)))
        dt = min(float(c.get("cfl", 0.4)) * g["dx"] / cmax, tend - self.t, float(c.get("dt_max_s", 5.0)))
        pool = c["outlet"].get("stage_m", -1e9) if c["outlet"]["type"] == "tailwater" else -1e9
        channel_step_q(self.h, self.Q, g["w"], g["z"], g["zs"], self.T, g["dx"], dt, float(c.get("n_manning", 0.04)), 0.0, qin, tin, pool, float(c.get("nu_factor", 1.0)),
                       melting, bool(c.get("heat_dissipation", True)), bool(c.get("evolve_width", True)), out)
        if melting:
            for k, i in (("V_in", 0), ("V_out", 1), ("sens_in", 2), ("heat_ice", 3), ("sens_out", 4), ("melt_w", 5), ("diss", 6), ("melted_ice", 7)): self.tot[k] += out[i]
            h = np.maximum(self.h, 1e-9); w = g["w"]; v = np.abs(self.Q) / (w * h); dh = 4 * w * h / (w + 2 * h); fr = v / np.sqrt(G * h); re = v * dh / NU_W
            wetm = self.h > 1e-4
            if wetm.any():
                f = self.flags; f["min_depth_m"] = min(f["min_depth_m"], float(self.h[wetm].min())); f["max_Fr"] = max(f["max_Fr"], float(fr[wetm].max())); f["min_Fr"] = min(f["min_Fr"], float(fr[wetm].min()))
                f["min_Re_Dh"] = min(f["min_Re_Dh"], float(re[wetm].min())); f["max_Re_Dh"] = max(f["max_Re_Dh"], float(re[wetm].max()))
            if (~wetm).any(): self.flags["dry_cell_seen"] = True
            if tin < 0: self.flags["negative_inlet_T"] = True
        self.t += dt; self.nsteps += 1
        return qin, tin, dt

    def record(self):
        g = self.g; self.rec.append(dict(t=self.t, stage_in=float(g["z"][0] + self.h[0]), stage_out=float(g["z"][-1] + self.h[-1]), Q_in=float(interp(self.case["inflow"]["Q"], self.t)), Q_out=float(self.Q[-1]), T_out_bulk_C=float(self.T[-1]),
                         P_out_W=float(RHO_W * CP_W * self.Q[-1] * self.T[-1]), melted_ice_m3=self.tot["melted_ice"], heat_to_ice_J=self.tot["heat_ice"], bed_lowering_max_m=float(np.max(g["z0"] - g["z"])), bed_lowering_mean_m=float(np.mean(g["z0"] - g["z"]))))

    def residuals(self):
        t = self.tot; dV = self.state_V() - self.V0; dE = self.state_E() - self.E0
        rw = (t["V_in"] + t["melt_w"] - t["V_out"] - dV); re = (t["sens_in"] + t["diss"] - t["sens_out"] - t["heat_ice"] - dE)
        g = self.g; geom = float(np.sum(((g["zs"] - g["z"]) * g["w"] - (g["zs"] - g["z0"]) * g["w0"])) * g["dx"])
        return dict(water=dict(residual_m3=rw, denominator_m3=t["V_in"], rel=rw / max(abs(t["V_in"]), 1e-30)), energy=dict(residual_J=re, denominator_J=max(abs(t["sens_in"]), abs(t["heat_ice"]), 1e-30), rel=re / max(abs(t["sens_in"]), abs(t["heat_ice"]), 1e-30)),
                    mass_ice_vs_meltwater=dict(rho_i_x_melted_ice_kg=RHO_I * t["melted_ice"], rho_w_x_meltwater_kg=RHO_W * t["melt_w"]), geometry_vs_latent_rel=(geom - t["melted_ice"]) / max(t["melted_ice"], 1e-30),
                    latent_energy_J=RHO_I * LF * t["melted_ice"], heat_to_ice_J=t["heat_ice"], power_identity_rel=(RHO_I * LF * t["melted_ice"] - t["heat_ice"]) / max(t["heat_ice"], 1e-30))

    def save(self, path, snaps):
        g = self.g; tmp = path + f".tmp{os.getpid()}.npz"
        np.savez(tmp, h=self.h, Q=self.Q, T=self.T, z=g["z"], w=g["w"], t=self.t, tot=json.dumps(self.tot), flags=json.dumps(self.flags), E0=self.E0, V0=self.V0, next_rec=self.next_rec, nsteps=self.nsteps, rec=json.dumps(self.rec), snaps=json.dumps(snaps))
        os.replace(tmp, path)

    def load(self, path):
        d = np.load(path, allow_pickle=False); g = self.g
        self.h = d["h"].copy(); self.Q = d["Q"].copy(); self.T = d["T"].copy(); g["z"] = d["z"].copy(); g["w"] = d["w"].copy(); self.t = float(d["t"]); self.tot = json.loads(str(d["tot"])); self.flags = json.loads(str(d["flags"]))
        self.E0 = float(d["E0"]); self.V0 = float(d["V0"]); self.next_rec = float(d["next_rec"]); self.nsteps = int(d["nsteps"]); self.rec = json.loads(str(d["rec"])); return json.loads(str(d["snaps"]))


def main(argv):
    case = json.load(open(argv[1])); out = argv[2]; refuse(case); os.makedirs(out, exist_ok=True); R = Run(case); snaps = []
    restart = argv[argv.index("--restart") + 1] if "--restart" in argv else None; stop_at = float(argv[argv.index("--stop-at") + 1]) if "--stop-at" in argv else None
    wall0 = time.time(); tend = float(case["t_end_s"]); rec_dt = float(case.get("record_every_s", 60.0)); snap_dt = float(case.get("snapshot_every_s", 3600.0)); ck_dt = float(case.get("checkpoint_every_s", 1800.0))
    if restart: snaps = R.load(restart)
    else:
        sp = float(case.get("spinup_s", 600.0)); ts = 0.0                    # spin-up with FROZEN geometry and no melting: documented, not part of the forecast clock
        while ts < sp:
            R.t = 0.0; ts += R.step(False, 1e30)[2]
        R.nsteps = 0
        R.t = 0.0; R.E0 = R.state_E(); R.V0 = R.state_V(); R.next_rec = 0.0
        # record geometry at the forecast origin
        snaps.append(dict(t=0.0, z=R.g["z"].tolist(), w=R.g["w"].tolist(), h=R.h.tolist(), T=R.T.tolist()))
    next_snap = (len(snaps)) * snap_dt if not restart else (int(R.t // snap_dt) + 1) * snap_dt; next_ck = R.t + ck_dt
    while R.t < tend - 1e-9:
        if R.t >= R.next_rec - 1e-9: R.record(); R.next_rec += rec_dt
        lim = min(tend, next_snap if next_snap > R.t else tend, R.next_rec if R.next_rec > R.t else tend)
        if stop_at is not None: lim = min(lim, stop_at)
        R.step(True, lim)
        if R.t >= next_snap - 1e-9: snaps.append(dict(t=R.t, z=R.g["z"].tolist(), w=R.g["w"].tolist(), h=R.h.tolist(), T=R.T.tolist())); next_snap += snap_dt
        if R.t >= next_ck: R.save(os.path.join(out, "checkpoint.npz"), snaps); next_ck += ck_dt
        if stop_at is not None and R.t >= stop_at - 1e-9: R.save(os.path.join(out, "checkpoint.npz"), snaps); print("stopped at", R.t); return
    R.record(); snaps.append(dict(t=R.t, z=R.g["z"].tolist(), w=R.g["w"].tolist(), h=R.h.tolist(), T=R.T.tolist())); R.save(os.path.join(out, "checkpoint.npz"), snaps)
    keys = list(R.rec[0].keys()); open(os.path.join(out, "timeseries.csv"), "w").write(",".join(keys) + "\n" + "\n".join(",".join("%.10g" % r[k] for k in keys) for r in R.rec) + "\n")
    g = R.g; ev = bool(case.get("evolve_width", True)); sn = dict(t=np.array([s["t"] for s in snaps]), x=g["x"], bed_lowering_m=np.array([np.array(snaps[0]["z"]) - np.array(s["z"]) for s in snaps]), depth_m=np.array([s["h"] for s in snaps]), T_bulk_C=np.array([s["T"] for s in snaps]))
    if ev: sn["width_m"] = np.array([s["w"] for s in snaps])
    np.savez(os.path.join(out, "snapshots.npz"), **sn)
    f = R.flags; flags = dict(f); flags["laminar_or_transitional_heat_transfer_closure_invalid"] = bool(f["min_Re_Dh"] < 1e4); flags["dittus_boelter_nominal_range_Re_ge_1e4"] = not flags["laminar_or_transitional_heat_transfer_closure_invalid"]
    flags["supercritical_throughout"] = bool(f["min_Fr"] > 1); flags["jump_possible_unresolved_interior"] = bool(f["max_Fr"] > 1 and f["min_Fr"] < 1); flags["depth_averaged_hydraulics"] = True
    flags["width_is_output"] = ev; flags["inlet_discharge_is_prescribed_forcing"] = True; flags["unsupported_processes"] = ["cold ice conduction", "debris", "roof closure", "3-D turbulence", "air entrainment", "atmospheric heat flux"]
    try: gh = subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=HERE, text=True).strip()
    except Exception: gh = "unknown"
    summ = dict(version=VERSION, git=gh, case_sha256=hashlib.sha256(json.dumps(case, sort_keys=True).encode()).hexdigest(), case_id=case.get("id"), t_end_s=R.t, steps=R.nsteps, wall_s=time.time() - wall0,
                outputs=dict(bed_lowering_max_m=R.rec[-1]["bed_lowering_max_m"], bed_lowering_mean_m=R.rec[-1]["bed_lowering_mean_m"], melted_ice_volume_m3=R.tot["melted_ice"], heat_to_ice_J=R.tot["heat_ice"], outlet_sensible_export_J=R.tot["sens_out"], delivered_volume_m3=R.tot["V_out"], inlet_volume_m3=R.tot["V_in"]),
                ledgers=R.residuals(), flags=flags)
    json.dump(summ, open(os.path.join(out, "summary.json"), "w"), indent=1, default=float); print(json.dumps(summ["outputs"]), json.dumps(summ["ledgers"]["energy"]), json.dumps(summ["ledgers"]["water"]))


if __name__ == "__main__": main(sys.argv)
