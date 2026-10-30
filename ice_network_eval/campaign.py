"""Observable, recoverable campaign runner. One run = one subprocess + one immutable directory.

  python campaign.py cases frozen_s1_s5            # write cases/<name>.json (list of case specs)
  python campaign.py drive cases/<name>.json [--workers 2]
  python campaign.py manifest                      # rebuild records/manifest.json from per-run dirs (never edited by workers)

A run directory runs/<id>/ holds spec.json (immutable), progress.json (atomic, overwritten by that run only),
record.json (written once at termination, status in {completed, failed, timeout, stalled}) and series.npz.
Failed/timed-out cases are NOT retried automatically; a retry needs a new spec (e.g. a version tag) and so a new id."""
import hashlib, json, os, resource, subprocess, sys, time
HERE = os.path.dirname(os.path.abspath(__file__)); RUNS = os.path.join(HERE, "records", "runs")
LIMITS = dict(reference_wall_s=1200, reduced_wall_s=900, stall_s=240)      # documented limits (stall = no progress update)


def commit():
    try: return subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=HERE, text=True).strip()
    except Exception: return "unknown"


def spec_id(spec):
    h = hashlib.sha256(json.dumps(spec, sort_keys=True).encode()).hexdigest()[:10]
    return f"{spec['label']}_{h}"


def atomic_write(path, obj):
    tmp = path + f".tmp{os.getpid()}"
    with open(tmp, "w") as f: json.dump(obj, f, indent=1, default=float)
    os.replace(tmp, path)


def applic(m, scn):
    """compact applicability summary (diagnostics only, from the sampled run); percentiles over active-edge snapshots with Q > 1 % of the mean inflow"""
    import numpy as np
    from ice_network.diagnostics.applicability import diagnose
    try:
        d = diagnose(m["net"], m["result"]); rows = [r for r in d["rows"] if r["Q"] > 0.01 * scn["inflow_m3_s"]]
        out = dict(summary=d.get("summary"), n_rows_active=len(rows))
        for k in ("r_qs", "froude_max", "channel_over_node", "melt_over_flow", "Re_min", "Re_max"):
            v = np.array([r[k] for r in rows if np.isfinite(r[k])]); out[k + "_p50_p90_max"] = [float(np.percentile(v, 50)), float(np.percentile(v, 90)), float(v.max())] if v.size else None
        out["frac_active_snapshots_with_froude_ge_1"] = float(np.mean([r["n_stations_froude_ge_1"] > 0 for r in rows])) if rows else None
        return out
    except Exception as e:
        return dict(error=f"{type(e).__name__}: {e}")


def worker(rd):
    sys.path.insert(0, HERE)
    spec = json.load(open(os.path.join(rd, "spec.json"))); t_start = time.time(); t0 = time.perf_counter()
    prog_path = os.path.join(rd, "progress.json"); last = [0.0]; state = dict(n_updates=0)
    def progress(t, n, dt=None, resid=None):
        now = time.perf_counter()
        if now - last[0] < 3.0: return
        last[0] = now; state["n_updates"] += 1
        atomic_write(prog_path, dict(status="running", sim_t_s=t, count=n, count_kind="rhs_evals_cumulative_per_macrostep" if spec["kind"] == "reduced" else "steps",
                                     dt=dt, min_depth_m=resid, wall_s=now - t0, last_progress_unix=time.time(), start_unix=t_start,
                                     peak_rss_mb=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024))
    import numpy as np
    t_imp = time.perf_counter(); import run_eval as R; t_imp = time.perf_counter() - t_imp
    if spec.get("cfg_patch"): R.set_config(spec["cfg_patch"])
    scn = next(s for s in R.CFG["scenarios"] if s["id"] == spec["scn"]) if "scn_override" not in spec else spec["scn_override"]
    ov = spec.get("overrides", {}); rec = dict(spec=spec, id=os.path.basename(rd), commit=spec.get("commit"), start_unix=t_start)
    timing = dict(import_s=t_imp)
    # warm-up (numba cache load / dispatcher compile): a 120 s simulation of the same problem, timed separately and discarded
    tw = time.perf_counter(); T_full = R.T_END; R.T_END = 120.0
    try:
        if spec["kind"] == "reduced": R.run_reduced(scn, spec["stations"], 60.0, sample_dt=None, **{k: v for k, v in ov.items() if k in ("n_man", "nu_factor", "msub", "rtol")})
        else: R.run_reference(scn, min(spec["ncell"], 50), **{k: v for k, v in ov.items() if k in ("n_man", "nu_factor", "cfl")})
    except BaseException as e:
        timing["warmup_error"] = f"{type(e).__name__}: {e}"
    R.T_END = T_full; timing["warmup_s"] = time.perf_counter() - tw
    try:
        if spec["kind"] == "reduced":
            m = R.run_reduced(scn, spec["stations"], spec["dt_geom"], timing=timing, progress=progress, **ov)
            timing["solver_s"] = m["wall"]; extra = dict(applicability=applic(m, scn), discharge_law_evaluations=m["nfev"], conservation={k: m["budgets"][k] for k in ("water_residual_rel", "mass_residual_rel", "energy_residual_rel", "geometry_vs_latent_rel_max", "geometry_vs_latent_kg_abs_max")},
                                                       n_macrosteps=int(len(m["result"]["t"]) and round(R.T_END / spec["dt_geom"])))
        else:
            m = R.run_reference(scn, spec["ncell"], progress=progress, **ov)
            timing["solver_s"] = m["wall"]; extra = dict(nsteps=m["nsteps"], diss_J=float(m["diss"]))
        tm = time.perf_counter(); met = R.metrics(m, scn)
        np.savez_compressed(os.path.join(rd, "series.npz"), t=m["t"], Q1=m["Q1"], Q2=m["Q2"], stage=m["stage"], sill1=m["sill1"], sill2=m["sill2"])
        timing["metrics_serialization_s"] = time.perf_counter() - tm
        rec.update(status="completed", metrics=met, timing=timing, **extra)
    except BaseException as e:
        rec.update(status="failed", error=f"{type(e).__name__}: {e}", timing=timing)
    rec.update(wall_total_s=time.perf_counter() - t0, peak_rss_mb=resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024, end_unix=time.time())
    if os.path.exists(os.path.join(rd, "record.json")): sys.exit("record exists; refusing to overwrite")
    atomic_write(os.path.join(rd, "record.json"), rec)


def drive(cases_path, workers=2):
    cases = json.load(open(cases_path)); os.makedirs(RUNS, exist_ok=True); sha = commit(); pending = []
    for c in cases:
        c = dict(c); c["commit"] = sha_for(c, sha); rid = spec_id({k: v for k, v in c.items() if k != "commit"}); rd = os.path.join(RUNS, rid)
        if os.path.exists(os.path.join(rd, "record.json")): continue
        pending.append((rid, rd, c))
    print(f"{len(cases)} cases, {len(pending)} pending, workers={workers}", flush=True); active = {}
    while pending or active:
        while pending and len(active) < workers:
            rid, rd, c = pending.pop(0); os.makedirs(rd, exist_ok=True)
            if os.path.exists(os.path.join(rd, "progress.json")): os.remove(os.path.join(rd, "progress.json"))
            atomic_write(os.path.join(rd, "spec.json"), c)
            p = subprocess.Popen([sys.executable, "-u", os.path.abspath(__file__), "worker", rd], stdout=open(os.path.join(rd, "log.txt"), "w"), stderr=subprocess.STDOUT)
            active[rid] = (p, rd, c, time.time()); print("start", rid, flush=True)
        time.sleep(2)
        for rid in list(active):
            p, rd, c, t0 = active[rid]; rc = p.poll(); now = time.time()
            lim = LIMITS["reference_wall_s" if c["kind"] == "reference" else "reduced_wall_s"]
            reason = None
            if rc is None:
                pp = os.path.join(rd, "progress.json")
                if now - t0 > lim: reason = "timeout"
                elif os.path.exists(pp) and now - json.load(open(pp)).get("last_progress_unix", now) > LIMITS["stall_s"]: reason = "stalled"
                elif not os.path.exists(pp) and now - t0 > LIMITS["stall_s"]: reason = "stalled"
                if reason:
                    p.kill(); p.wait()
            if rc is None and not reason: continue
            if not os.path.exists(os.path.join(rd, "record.json")):
                pr = json.load(open(os.path.join(rd, "progress.json"))) if os.path.exists(os.path.join(rd, "progress.json")) else None
                atomic_write(os.path.join(rd, "record.json"), dict(spec=c, id=rid, status=reason or "failed", exit_code=rc, limits=LIMITS, wall_total_s=now - t0, last_progress=pr,
                             error=None if reason else open(os.path.join(rd, "log.txt")).read()[-1500:]))
            r = json.load(open(os.path.join(rd, "record.json"))); print("end", rid, r["status"], round(r.get("wall_total_s", 0), 1), flush=True); del active[rid]
    manifest()


def sha_for(c, sha): return c.get("commit", sha)


def manifest():
    rows = []
    for rid in sorted(os.listdir(RUNS)) if os.path.isdir(RUNS) else []:
        rd = os.path.join(RUNS, rid); rp = os.path.join(rd, "record.json"); pp = os.path.join(rd, "progress.json")
        if os.path.exists(rp):
            r = json.load(open(rp)); rows.append(dict(id=rid, state=r["status"], wall_s=r.get("wall_total_s"), peak_rss_mb=r.get("peak_rss_mb"), label=r["spec"]["label"]))
        elif os.path.exists(pp): rows.append(dict(id=rid, state="running_or_orphaned", progress=json.load(open(pp))))
        else: rows.append(dict(id=rid, state="pending"))
    from collections import Counter
    out = dict(generated_unix=time.time(), counts=dict(Counter(r["state"] for r in rows)), total_wall_s=sum(r.get("wall_s") or 0 for r in rows), runs=rows)
    atomic_write(os.path.join(HERE, "records", "manifest.json"), out); print(out["counts"], "total wall s", round(out["total_wall_s"]))


PARAM_RANGES = {   # name: (low, high, scale, status).  Idealised bare-ice ranges bracketing order-of-magnitude supraglacial values; NOT tied to any glacier and no literature value was verified in this work.
    "inflow_m3_s": (0.01, 0.10, "log"), "pulse_amp": (0.0, 0.6, "lin"), "pulse_period_h": (1.0, 4.0, "lin"), "inflow_T_C": (0.0, 2.0, "lin"),
    "deficit_m": (0.01, 0.30, "log"), "area_m2": (1000.0, 6000.0, "log"), "sill_offset_m": (0.0, 0.15, "lin"),
    "w1_m": (0.2, 1.0, "lin"), "L1_m": (50.0, 200.0, "lin"), "slope1": (0.002, 0.02, "log"), "slope2": (0.005, 0.05, "log"), "L2_m": (30.0, 120.0, "lin"),
    "n_man": (0.03, 0.05, "lin"), "nu_factor": (0.5, 2.0, "log")}


CHALLENGE_RANGES = {   # deliberately wider than PARAM_RANGES (outside the development envelope on every axis); negative deficit = lake initially ABOVE sill 1 (initially wet route 1); negative sill_offset = route 2 is the LOWER sill
    "inflow_m3_s": (0.005, 0.30, "log"), "pulse_amp": (0.0, 0.9, "lin"), "pulse_period_h": (0.25, 8.0, "log"), "inflow_T_C": (0.0, 4.0, "lin"),
    "deficit_m": (-0.10, 0.50, "lin"), "area_m2": (300.0, 12000.0, "log"), "sill_offset_m": (-0.15, 0.20, "lin"),
    "w1_m": (0.1, 2.0, "log"), "L1_m": (20.0, 400.0, "log"), "slope1": (0.0005, 0.08, "log"), "slope2": (0.002, 0.15, "log"), "L2_m": (20.0, 300.0, "log"),
    "n_man": (0.02, 0.07, "lin"), "nu_factor": (0.33, 3.0, "log")}


def param_point(u, ranges=None):
    import math
    p = {}
    for x, (k, (lo, hi, sc)) in zip(u, (ranges or PARAM_RANGES).items()):
        p[k] = float(math.exp(math.log(lo) + x * (math.log(hi) - math.log(lo)))) if sc == "log" else float(lo + x * (hi - lo))
    return p


def param_spec(p, idx, set_name, kind, **kw):
    patch = dict(lake=dict(area_m2=p["area_m2"], initial_stage_below_sill1_m=p["deficit_m"]),
                 channels=[dict(width_m=p["w1_m"], length_m=p["L1_m"], slope=p["slope1"]), dict(width_m=0.4, length_m=p["L2_m"], slope=p["slope2"])])
    scn = dict(id=f"{set_name}{idx:03d}", sill_offset_m=p["sill_offset_m"], inflow_m3_s=p["inflow_m3_s"], inflow_T_C=p["inflow_T_C"],
               pulse_amp=p["pulse_amp"], pulse_period_s=p["pulse_period_h"] * 3600.0)
    return dict(kind=kind, scn=scn["id"], scn_override=scn, cfg_patch=patch, overrides=dict(n_man=p["n_man"], nu_factor=p["nu_factor"]), params=p, label=f"{kind[:3]}_{set_name}{idx:03d}" + kw.pop("suffix", ""), version="v1.0", **kw)


def param_cases(set_name, seed, n, ranges=None):
    from scipy.stats import qmc
    U = qmc.Sobol(d=len(PARAM_RANGES), scramble=True, seed=seed).random(n)
    out = []
    for i, u in enumerate(U):
        p = param_point(u, ranges)
        for st, dtg in ((9, 900), (17, 450)):
            out.append(param_spec(p, i, set_name, "reduced", stations=st, dt_geom=dtg, suffix=f"_r{st}"))
        out.append(param_spec(p, i, set_name, "reference", ncell=200, suffix="_200"))
        if i % 4 == 0:
            out.append(param_spec(p, i, set_name, "reference", ncell=100, suffix="_100"))
            out.append(param_spec(p, i, set_name, "reference", ncell=400, suffix="_400"))
    return out


def cases(name):
    sys.path.insert(0, HERE); import run_eval as R
    out = []
    if name == "frozen_s1_s5":      # same scenarios/resolutions/tolerances as the original frozen evaluation, re-executed under the observable runner (reduced_v1.0 / reference_v1.0 physics)
        for s in R.CFG["scenarios"]:
            for n in (100, 200, 400): out.append(dict(kind="reference", scn=s["id"], ncell=n, label=f"ref_{s['id']}_{n}", version="v1.0"))
            for rc in R.CFG["reduced_model_configurations"]:
                out.append(dict(kind="reduced", scn=s["id"], stations=rc["stations"], dt_geom=rc["dt_geom_s"], label=f"red_{s['id']}_{rc['id']}", version="v1.0"))
    if name.startswith("param_"):    # param_dev (seed 11, 32 pts), param_eval (seed 22, 64 pts): independent Sobol designs, deterministic
        out = param_cases({"param_dev": "D", "param_eval": "E"}[name], {"param_dev": 11, "param_eval": 22}[name], {"param_dev": 32, "param_eval": 64}[name])
    if name == "convergence":   # one factor at a time about the baseline r17 / dt 450 / msub 16 / rtol 1e-7; reference mesh and CFL separately
        for sid in ("S1", "S2", "S3", "S4"):    # S1 competing outlets, S2 near activation, S3 strongly incising (largest inflow), S4 weak melting / smooth
            base = dict(stations=17, dt_geom=450, msub=16, rtol=1e-7)
            def red(tag, **kw):
                d = dict(base); d.update(kw); ov = {k: d[k] for k in ("msub", "rtol")}
                out.append(dict(kind="reduced", scn=sid, stations=d["stations"], dt_geom=d["dt_geom"], overrides=ov, label=f"red_{sid}_conv_{tag}", version="v1.0"))
            for st in (5, 9, 17, 33, 65): red(f"st{st}", stations=st)
            for dtg in (1800, 900, 450, 225, 112.5): red(f"dt{dtg}", dt_geom=dtg)
            for ms in (2, 4, 8, 16, 32): red(f"msub{ms}", msub=ms)
            for rt in (1e-5, 1e-6, 1e-7, 1e-8, 1e-9): red(f"rtol{rt:g}", rtol=rt)
            for n in (100, 200, 400, 800): out.append(dict(kind="reference", scn=sid, ncell=n, label=f"ref_{sid}_conv_{n}", version="v1.0"))
            for cfl in (0.2, 0.1): out.append(dict(kind="reference", scn=sid, ncell=200, overrides=dict(cfl=cfl), label=f"ref_{sid}_conv_200_cfl{cfl}", version="v1.0"))
        seen = set(); uniq = []
        for o in out:
            k = json.dumps({a: b for a, b in o.items() if a != "label"}, sort_keys=True)
            if k not in seen: seen.add(k); uniq.append(o)
        out = uniq
    if name == "challenge":   # frozen BEFORE any challenge run (commit precedes execution): Sobol seed 33, 32 points, outside the development ranges
        out = param_cases("C", 33, 32, CHALLENGE_RANGES)
    if name == "closure_sensitivity":   # closure variants applied identically to both models (matched) so closure sensitivity and numerical error are separable
        variants = {"nu0.5": dict(nu_factor=0.5), "nu2": dict(nu_factor=2.0), "n0.03": dict(n_man=0.03), "n0.05": dict(n_man=0.05)}
        for s in R.CFG["scenarios"]:
            for vn, ov in variants.items():
                out.append(dict(kind="reference", scn=s["id"], ncell=200, overrides=ov, label=f"ref_{s['id']}_200_{vn}", version="v1.0"))
                for st, dtg in ((9, 900), (17, 450)):
                    out.append(dict(kind="reduced", scn=s["id"], stations=st, dt_geom=dtg, overrides=ov, label=f"red_{s['id']}_r{st}_{vn}", version="v1.0"))
    if name == "timing_serial":   # repeated timings, run with --workers 1 on an otherwise idle machine
        for s in ("S1", "S2", "S5"):
            for rep in range(3):
                for n in (100, 200): out.append(dict(kind="reference", scn=s, ncell=n, label=f"ref_{s}_{n}_rep{rep}", version="v1.0"))
                for st, dtg in ((9, 900), (17, 450), (33, 225)): out.append(dict(kind="reduced", scn=s, stations=st, dt_geom=dtg, label=f"red_{s}_r{st}_rep{rep}", version="v1.0"))
    os.makedirs(os.path.join(HERE, "cases"), exist_ok=True); atomic_write(os.path.join(HERE, "cases", name + ".json"), out); print(len(out), "cases")


if __name__ == "__main__":
    c = sys.argv[1]
    if c == "worker": worker(sys.argv[2])
    elif c == "drive": drive(sys.argv[2], int(sys.argv[sys.argv.index("--workers") + 1]) if "--workers" in sys.argv else 2)
    elif c == "manifest": manifest()
    elif c == "cases": cases(sys.argv[2])
