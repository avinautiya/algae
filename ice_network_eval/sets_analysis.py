"""Per-point analysis of the parameter/challenge sets from per-run records, with the declared applicability flags and failure classes (cases/CHALLENGE_PROTOCOL.md).
python sets_analysis.py param_dev|param_eval|challenge"""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import analyze_campaign as A
TOL_HYD = 0.10; TIE = 0.10


def flags(app):
    if not app or "error" in app or app.get("r_qs_p50_p90_max") is None: return None
    return dict(quasi_steady=bool(app["r_qs_p50_p90_max"][1] <= 1.0), subcritical=bool((app.get("frac_active_snapshots_with_froude_ge_1") or 0.0) <= 0.05),
                storage_neglected=bool((app["channel_over_node_p50_p90_max"] or [0, 0, 0])[1] <= 0.05))


def analyze(setname):
    runs = A.load(lambda r: A.case_set(r["spec"]) == setname)
    pts = {}
    for r in runs:
        s = r["spec"]; pts.setdefault(s["scn"], dict(params=s["params"], scn=s["scn_override"], ref={}, red={}))
        (pts[s["scn"]]["ref"].__setitem__(s["ncell"], r) if s["kind"] == "reference" else pts[s["scn"]]["red"].__setitem__(s["stations"], r))
    out = {}
    for pid, p in sorted(pts.items()):
        I = p["scn"]["inflow_m3_s"]; refs = {n: r for n, r in p["ref"].items() if r["status"] == "completed"}
        row = dict(params=p["params"], ref_status={n: r["status"] for n, r in p["ref"].items()}, red_status={n: r["status"] for n, r in p["red"].items()})
        if not refs: row["class"] = "reference_failed"; out[pid] = row; continue
        top = max(refs); tgt = refs[top]; ts = A.series(tgt); tm = tgt["metrics"]
        v1, v2 = tm["delivered1_m3"], tm["delivered2_m3"]; margin = abs(v1 - v2) / max(v1 + v2, 1e-12)
        row.update(target_ncell=top, route_ref=A.route(tm), route_margin=margin, t_act2_ref=tm["t_act2_s"], t_dom2_ref=tm["t_dom2_s"], sill_lowering_ref=[tm["sill_lowering1_m"], tm["sill_lowering2_m"]], delivered_ref=[v1, v2],
                   ref_hydro_change=(A.ser_err(A.series(refs[sorted(refs)[-2]]), ts, I)["hydrograph_rmse_over_I"] if len(refs) > 1 else None), reds={})
        for st, r in sorted(p["red"].items()):
            if r["status"] != "completed":
                fl = None; row["reds"][st] = dict(status=r["status"], last_progress=r.get("last_progress"), error=r.get("error"), class_="failure_numerical")
                continue
            m = r["metrics"]; e = A.ser_err(A.series(r), ts, I); fl = flags(r.get("applicability")); inside = bool(fl and all(fl.values()))
            rm = A.route(m) == A.route(tm); near = margin <= TIE
            ev = lambda k: None if (m[k] is None or tm[k] is None) else m[k] - tm[k]
            bad = (e["hydrograph_rmse_over_I"] > TOL_HYD) or ((not rm) and not near)
            cls = "ok" if not bad else ("failure_unexplained" if inside else "failure_outside_regime")
            row["reds"][st] = dict(status="completed", route=A.route(m), route_match=rm, near_tie=near, hydro_rmse_over_I=e["hydrograph_rmse_over_I"], cumvol_err_rel_total=e["cumvol_err_rel_total"], t_act2_err_s=ev("t_act2_s"),
                                  switch_ref=tm["t_dom2_s"] is not None, switch_red=m["t_dom2_s"] is not None, t_dom2_err_s=ev("t_dom2_s"),
                                  sill_err_m=[m["sill_lowering1_m"] - tm["sill_lowering1_m"], m["sill_lowering2_m"] - tm["sill_lowering2_m"]], flags=fl, inside_envelope=inside,
                                  conservation=r.get("conservation"), solver_s=r["timing"]["solver_s"], class_=cls, applicability=r.get("applicability", {}).get("summary") if r.get("applicability") else None)
        out[pid] = row
    return out


def aggregate(res):
    agg = {}
    for st in (9, 17):
        rows = [v["reds"][st] for v in res.values() if "reds" in v and st in v["reds"]]
        comp = [x for x in rows if x["status"] == "completed"]; cls = {}
        for x in rows: cls[x["class_"]] = cls.get(x["class_"], 0) + 1
        clear = [x for x in comp if not x["near_tie"]]
        sw = [x for x in comp if x["switch_ref"] or x["switch_red"]]
        agg[f"r{st}"] = dict(n=len(rows), classes=cls, hydro_rmse_median=float(np.median([x["hydro_rmse_over_I"] for x in comp])) if comp else None, hydro_rmse_p90=float(np.percentile([x["hydro_rmse_over_I"] for x in comp], 90)) if comp else None,
                             route_accuracy_clear_cases=(float(np.mean([x["route_match"] for x in clear])), len(clear)) if clear else None, n_near_tie=len(comp) - len(clear),
                             switching_agree=(int(sum(x["switch_ref"] == x["switch_red"] for x in sw)), len(sw)), inside_envelope=int(sum(x["inside_envelope"] for x in comp)), outside_envelope=int(sum(not x["inside_envelope"] for x in comp)))
    agg["reference_failed"] = int(sum(v.get("class") == "reference_failed" for v in res.values()))
    return agg


if __name__ == "__main__":
    res = analyze(sys.argv[1]); a = aggregate(res)
    json.dump(dict(points=res, aggregate=a), open(os.path.join(HERE, "records", f"set_summary_{sys.argv[1]}.json"), "w"), indent=1, default=float)
    print(len(res), "points"); print(json.dumps(a, indent=1, default=float))
