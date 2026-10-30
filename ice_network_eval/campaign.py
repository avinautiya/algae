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
    scn = next(s for s in R.CFG["scenarios"] if s["id"] == spec["scn"]) if "scn_override" not in spec else spec["scn_override"]
    ov = spec.get("overrides", {}); rec = dict(spec=spec, id=os.path.basename(rd), commit=spec.get("commit"), start_unix=t_start)
    timing = dict(import_s=t_imp)
    try:
        if spec["kind"] == "reduced":
            tw = time.perf_counter(); R.run_reduced(scn, 3, 86400.0, msub=2, sample_dt=None, **{k: v for k, v in ov.items() if k in ("n_man", "nu_factor")}) if False else None
            m = R.run_reduced(scn, spec["stations"], spec["dt_geom"], timing=timing, progress=progress, **ov)
            timing["solver_s"] = m["wall"]; extra = dict(discharge_law_evaluations=m["nfev"], conservation={k: m["budgets"][k] for k in ("water_residual_rel", "mass_residual_rel", "energy_residual_rel", "geometry_vs_latent_rel_max", "geometry_vs_latent_kg_abs_max")},
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


def cases(name):
    sys.path.insert(0, HERE); import run_eval as R
    out = []
    if name == "frozen_s1_s5":      # same scenarios/resolutions/tolerances as the original frozen evaluation, re-executed under the observable runner (reduced_v1.0 / reference_v1.0 physics)
        for s in R.CFG["scenarios"]:
            for n in (100, 200, 400): out.append(dict(kind="reference", scn=s["id"], ncell=n, label=f"ref_{s['id']}_{n}", version="v1.0"))
            for rc in R.CFG["reduced_model_configurations"]:
                out.append(dict(kind="reduced", scn=s["id"], stations=rc["stations"], dt_geom=rc["dt_geom_s"], label=f"red_{s['id']}_{rc['id']}", version="v1.0"))
    if name == "closure_sensitivity":   # closure variants applied identically to both models (matched) so closure sensitivity and numerical error are separable
        variants = {"nu0.5": dict(nu_factor=0.5), "nu2": dict(nu_factor=2.0), "n0.03": dict(n_man=0.03), "n0.05": dict(n_man=0.05)}
        for s in R.CFG["scenarios"]:
            for vn, ov in variants.items():
                out.append(dict(kind="reference", scn=s["id"], ncell=200, overrides=ov, label=f"ref_{s['id']}_200_{vn}", version="v1.0"))
                for st, dtg in ((9, 900), (17, 450)):
                    out.append(dict(kind="reduced", scn=s["id"], stations=st, dt_geom=dtg, overrides=ov, label=f"red_{s['id']}_r{st}_{vn}", version="v1.0"))
    os.makedirs(os.path.join(HERE, "cases"), exist_ok=True); atomic_write(os.path.join(HERE, "cases", name + ".json"), out); print(len(out), "cases")


if __name__ == "__main__":
    c = sys.argv[1]
    if c == "worker": worker(sys.argv[2])
    elif c == "drive": drive(sys.argv[2], int(sys.argv[sys.argv.index("--workers") + 1]) if "--workers" in sys.argv else 2)
    elif c == "manifest": manifest()
    elif c == "cases": cases(sys.argv[2])
