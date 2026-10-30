"""Summaries generated ONLY from immutable per-run records (records/runs/*/record.json, series.npz).
  python analyze_campaign.py <prefix-or-all> <out.json>
Groups runs by (scenario, overrides); targets are the finest reference (400 cells) with the 200->400 change as its uncertainty estimate."""
import glob, json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); RUNS = os.path.join(HERE, "records", "runs")
KEYS = ("t_act1_s", "t_act2_s", "t_dom2_s", "sill_lowering1_m", "sill_lowering2_m", "delivered1_m3", "delivered2_m3", "heat_to_ice_J", "stage_max_m")


def load(selector=None):
    runs = []
    for f in sorted(glob.glob(os.path.join(RUNS, "*", "record.json"))):
        r = json.load(open(f))
        if selector and not selector(r): continue
        r["_dir"] = os.path.dirname(f); runs.append(r)
    return runs


def series(r):
    z = np.load(os.path.join(r["_dir"], "series.npz")); return {k: z[k] for k in z.files}


def route(m):   # outlet selection: route with the larger delivered volume over the run; "none" if neither delivered anything
    a, b = m["delivered1_m3"], m["delivered2_m3"]
    return "none" if a + b <= 0 else ("R1" if a >= b else "R2")


def ser_err(red, ref, I):
    t = ref["t"]; out = {}
    for k in ("Q1", "Q2"):
        d = np.interp(t, red["t"], red[k]) - ref[k]
        out[k] = dict(rmse_over_I=float(np.sqrt(np.mean(d ** 2)) / I), max_over_I=float(np.abs(d).max() / I), cumvol_err_m3=float(np.trapezoid(d, t)))
    tot = float(np.trapezoid(ref["Q1"] + ref["Q2"], t))
    out["cumvol_err_rel_total"] = float((abs(out["Q1"]["cumvol_err_m3"]) + abs(out["Q2"]["cumvol_err_m3"])) / max(tot, 1e-30))
    out["hydrograph_rmse_over_I"] = float(np.hypot(out["Q1"]["rmse_over_I"], out["Q2"]["rmse_over_I"]))
    return out


def summarize(runs):
    groups = {}
    for r in runs:
        s = r["spec"]; groups.setdefault((s["scn"], json.dumps(s.get("overrides", {}), sort_keys=True)), []).append(r)
    res = {}
    for (scn, ov), rs in groups.items():
        refs = {r["spec"]["ncell"]: r for r in rs if r["spec"]["kind"] == "reference"}
        reds = {r["spec"]["label"]: r for r in rs if r["spec"]["kind"] == "reduced"}
        entry = dict(overrides=json.loads(ov), failed=[dict(id=r["id"], status=r["status"], last_progress=r.get("last_progress"), error=r.get("error")) for r in rs if r["status"] != "completed"])
        ok = {n: r for n, r in refs.items() if r["status"] == "completed"}
        if not ok: res[f"{scn}|{ov}"] = entry; continue
        top = max(ok); tgt = ok[top]; I = None
        sc = r["spec"].get("scn_override"); I = sc["inflow_m3_s"] if sc else None
        if I is None:
            import run_eval as R; I = next(s for s in R.CFG["scenarios"] if s["id"] == scn)["inflow_m3_s"]
        tser = series(tgt); tm = tgt["metrics"]
        conv = {}
        ns = sorted(ok)
        for k in KEYS:
            vals = {n: ok[n]["metrics"][k] for n in ns}
            row = dict(values=vals)
            if len(ns) >= 3 and all(v is not None for v in vals.values()):
                d1 = vals[ns[1]] - vals[ns[0]]; d2 = vals[ns[2]] - vals[ns[1]]
                row.update(change_fine=d2, change_coarse=d1, observed_order=float(np.log2(abs(d1) / abs(d2))) if d1 * d2 > 0 and d2 != 0 else None, richardson_first_order=vals[ns[2]] + d2)
            conv[k] = row
        entry["reference_convergence"] = conv
        entry["reference_hydrograph_change"] = {f"{ns[i]}->{ns[i+1]}": ser_err(series(ok[ns[i]]), series(ok[ns[i + 1]]), I)["hydrograph_rmse_over_I"] for i in range(len(ns) - 1)}
        entry["target"] = dict(ncell=top, route=route(tm), metrics=tm, wall_s={n: ok[n]["timing"]["solver_s"] for n in ns}, steps={n: ok[n].get("nsteps") for n in ns})
        entry["reference_cost_vs_error"] = {n: dict(hydrograph_rmse_over_I=ser_err(series(ok[n]), tser, I)["hydrograph_rmse_over_I"], wall_s=ok[n]["timing"]["solver_s"]) for n in ns if n != top}
        entry["reduced"] = {}
        for lab, r in sorted(reds.items()):
            if r["status"] != "completed": continue
            m = r["metrics"]; e = ser_err(series(r), tser, I)
            ev = {}
            for k in ("t_act2_s", "t_dom2_s"):
                a, b = m[k], tm[k]
                ev[k] = None if (a is None and b is None) else ("missing_one_side" if (a is None or b is None) else a - b)
            entry["reduced"][r["spec"]["label"].split("_", 2)[-1] if r["spec"].get("scn_override") is None else r["spec"]["label"]] = dict(
                route=route(m), route_match=route(m) == route(tm), event_time_error_s=ev, hydrograph=e,
                sill_lowering_err_m=[m["sill_lowering1_m"] - tm["sill_lowering1_m"], m["sill_lowering2_m"] - tm["sill_lowering2_m"]],
                sill_lowering_rel_err=[(m["sill_lowering1_m"] - tm["sill_lowering1_m"]) / max(abs(tm["sill_lowering1_m"]), 1e-12), (m["sill_lowering2_m"] - tm["sill_lowering2_m"]) / max(abs(tm["sill_lowering2_m"]), 1e-12)],
                delivered_err_m3=[m["delivered1_m3"] - tm["delivered1_m3"], m["delivered2_m3"] - tm["delivered2_m3"]],
                conservation=r.get("conservation"), solver_s=r["timing"]["solver_s"], sampling_run_s=r["timing"].get("sampling_run_s"), peak_rss_mb=r["peak_rss_mb"], discharge_law_evals=r.get("discharge_law_evaluations"))
        res[f"{scn}|{ov}"] = entry
    return res


if __name__ == "__main__":
    pre = sys.argv[1]; sel = None if pre == "all" else (lambda r: r["spec"]["label"].split("_", 1)[1].startswith(pre) or r["spec"]["scn"].startswith(pre))
    out = summarize(load(sel)); json.dump(out, open(sys.argv[2], "w"), indent=1, default=float); print(len(out), "groups ->", sys.argv[2])
