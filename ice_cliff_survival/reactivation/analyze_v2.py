"""v2 analysis: identical-footprint melt/heat vs reduced model, corrected recovery classes (v1 + v2), removal-time bracket, counterfactual. `python -m ice_cliff_survival.reactivation.analyze_v2 V1_EVAL V2_EVAL SPEC_DIR OUT_DIR`"""
import glob, json, os, sys
import numpy as np
from . import metrics as M
def jl(p): return json.load(open(p))
def runs(d, pat):
    o = {}
    for f in sorted(glob.glob(os.path.join(d, pat))):
        r = jl(f)
        if r.get("status") == "complete": o[os.path.basename(f)[:-5]] = r
    return o
def main(v1, v2, specdir, out, tag="v2"):
    os.makedirs(out, exist_ok=True); spec = jl(os.path.join(specdir, f"production_spec_{tag}.json")); T_ = tag.upper(); P = spec["predictions"]; res = dict(); T = []
    V2 = runs(v2, T_ + "*.json"); V1m = runs(v1, "E*_medium.json"); V1c = runs(v1, "E*_coarse.json")
    fref = M.arrays(V1m["E01_bare_none_medium"])["f"].mean(); fref_c = M.arrays(V1c["E01_bare_none_coarse"])["f"].mean()
    # ---- 1. v2 coarse repeats vs v1 coarse: identical flow? and identical-footprint melt/heat vs reduced
    fp = []; same = {}
    for n in [e["name"] for e in spec["runs"] if e["role"] == "v2_coarse_repeat_with_regions"]:
        k = n; e = n[4:]
        if k not in V2: fp.append(dict(run=e, status="MISSING")); continue
        r = V2[k]; a = M.arrays(r); b = M.arrays(V1c[e + "_coarse"]); same[e] = float(np.max(np.abs(a["f"] - b["f"]))); last = r["records"][-1]; days = r["records"][-1]["day"]
        cfd = last["region_cum_melt_m3"]; heat = {kk: float(np.mean([x["region_heat_W"][kk] for x in r["records"]])) for kk in last["region_heat_W"]}; pred = P["footprint"].get(e)
        row = dict(run=e, loss=M.events(a["day"], a["f"], 1e-3, r.get("removal_day"))["loss"], f_end=float(a["f"][-1]), cfd_melt_m3={kk: cfd[kk] for kk in ("toe_fp", "toe_open", "byp")}, cfd_mean_heat_W={kk: heat[kk] for kk in ("toe_fp", "toe_open", "byp")}, ramp_melt_m3=dict(toe=cfd["ramp_toe"], byp=cfd["ramp_byp"]), region_area=r["region_area"])
        if pred: row["reduced_melt_m3"] = pred["melt_m3"]; row["reduced_mean_heat_W"] = pred["mean_heat_W"]; row["rel_err_melt"] = {kk: (pred["melt_m3"][kk] - cfd[kk]) / cfd[kk] if cfd[kk] > 1e-9 else None for kk in pred["melt_m3"]}
        fp.append(row)
    res["identical_footprint"] = fp; res["v2_vs_v1_max_f_difference_coarse"] = same
    # ---- 2. corrected classification of v1 runs (preserving the original label) + v2
    cl = []
    for n, r in V1m.items():
        e = n[:-7]; rem = r.get("removal_day"); a = M.arrays(r); c = M.classify(a["day"], a["f"], rem, f_ref=fref); s = M.summarize(dict(r, name=e)); cl.append(dict(run=e, removal_day=rem, original_label=("sustained reconnection at %.2f" % s["sustained_reconnection"]) if s["sustained_reconnection"] is not None else "no sustained reconnection", **c))
    res["corrected_classes_v1_medium"] = cl
    cl3 = []; f3 = M.arrays(V2[T_ + "C_E01_bare_none"])["f"].mean()
    for e in [x["name"] for x in spec["runs"] if x["role"] == "v2_coarse_repeat_with_regions"]:
        r = V2.get(e)
        if r: a = M.arrays(r); c = M.classify(a["day"], a["f"], r.get("removal_day"), f_ref=f3); cl3.append(dict(run=e[4:], removal_day=r.get("removal_day"), f_end=float(a["f"][-1]), **c))
    res["classes_solver_v3_coarse"] = cl3
    # ---- 3. bracket
    br = []
    for r_ in (8.5, 9.0, 9.5, 10.0, 10.5):
        k = f"{T_}M_strip_r{r_:g}"
        if k not in V2: br.append(dict(removal_day=r_, status="MISSING")); continue
        r = V2[k]; a = M.arrays(r); c = M.classify(a["day"], a["f"], r_, f_ref=fref); pr = P["bracket"][f"{r_:g}"]; br.append(dict(removal_day=r_, f_plus1d=c["f_post1"], cls=c["cls"], R=c["R"], f_end=float(a["f"][-1]), pred_f_plus1d=pr["f_at_removal_plus_1d"], pred_cls=pr["cls"], class_match=(pr["cls"] == c["cls"]), loss=M.events(a["day"], a["f"], 1e-3, r_)["loss"], region_heat_toe_fp_mean_W=float(np.mean([x["region_heat_W"]["toe_fp"] for x in r["records"]])), melt_toe_fp_m3=r["records"][-1]["region_cum_melt_m3"]["toe_fp"], melt_toe_open_m3=r["records"][-1]["region_cum_melt_m3"]["toe_open"]))
    pts = [(8.0, 0.104)] + [(b["removal_day"], b["f_plus1d"]) for b in br if "f_plus1d" in b] + [(12.0, 0.0)]; pts.sort(); tb = None
    for (x0, f0), (x1, f1) in zip(pts[:-1], pts[1:]):
        if f0 >= 1e-3 > f1: tb = x0 + (x1 - x0) * (f0 - 1e-3) / (f0 - f1); break
    res["bracket"] = br; res["cfd_boundary_day"] = tb; res["reduced_boundary_day"] = P["reduced_boundary_day"]; res["boundary_pass"] = (tb is not None and abs(tb - P["reduced_boundary_day"]) <= 1.0)
    # ---- 4. counterfactual
    cf = []
    for k, r in V2.items():
        if "_CF_" not in k: continue
        a = M.arrays(r); e = M.events(a["day"], a["f"], 1e-3, r.get("removal_day")); ref_name = {}.get(k)
        cf.append(dict(run=k, removal_day=r.get("removal_day"), incision=r.get("incision"), contact_lost=e["loss"] is not None, loss=e["loss"], f0=float(a["f"][0]), f_end=float(a["f"][-1]), min_f=float(a["f"].min()), toe_fp_melt_m3=r["records"][-1]["region_cum_melt_m3"]["toe_fp"], byp_melt_m3=r["records"][-1]["region_cum_melt_m3"]["byp"], crest_toe_end=float(a["crest_toe"][-1]), crest_byp_end=float(a["crest_byp"][-1])))
    res["counterfactual"] = cf; res["reduced_counterfactual"] = P["counterfactual"]
    json.dump(res, open(os.path.join(out, "analysis_v2.json"), "w"), indent=1, default=float)
    fm = lambda x, p=3: "—" if x is None else f"{x:.{p}f}"
    L = ["### Identical-footprint melt and heat (21 d, coarse dx=0.05, v2 runs; flow identical to v1 to %.1e)" % max(same.values() or [0]), "| run | loss day | f_end | CFD melt fp / open / byp [m³] | reduced melt fp / open / byp | rel. err fp / open / byp | CFD mean heat fp / open / byp [W] | reduced heat |", "|---|---|---|---|---|---|---|---|"]
    for r in fp:
        if r.get("status"): L.append(f"| {r['run']} | MISSING |"); continue
        cm = r["cfd_melt_m3"]; hm = r["cfd_mean_heat_W"]; rm = r.get("reduced_melt_m3"); rh = r.get("reduced_mean_heat_W"); er = r.get("rel_err_melt")
        L.append(f"| {r['run']} | {fm(r['loss'],2)} | {r['f_end']:.3f} | {cm['toe_fp']:.4f} / {cm['toe_open']:.4f} / {cm['byp']:.4f} | " + (f"{rm.get('toe_fp',0):.4f} / {rm.get('toe_open',0):.4f} / {rm.get('byp',0):.4f} | " + " / ".join(fm(None if er[k] is None else 100 * er[k], 1) + "%" if k in er else "n/a" for k in ("toe_fp", "toe_open", "byp")) + f" | {hm['toe_fp']:.1f} / {hm['toe_open']:.1f} / {hm['byp']:.1f} | " + f"{rh.get('toe_fp',0):.1f} / {rh.get('toe_open',0):.1f} / {rh.get('byp',0):.1f} |" if rm else f"n/a (no reduced prediction) | | {hm['toe_fp']:.1f} / {hm['toe_open']:.1f} / {hm['byp']:.1f} | |"))
    L += ["", "### Corrected recovery classes of the v1 evaluation runs (medium)", "| run | removal | original v1 label | corrected class | R = post-removal share / bare share | decline flag |", "|---|---|---|---|---|---|"] + [f"| {c['run']} | {fm(c['removal_day'],1)} | {c['original_label']} | **{c['cls']}** | {fm(c['R'])} | {c['decline']} |" for c in cl]
    L += ["", "### Recovery classes with the corrected solver v3 (coarse dx=0.05, 21 d); solver-v1 labels above are superseded where they differ", "| run | removal | class | R | decline | f_end |", "|---|---|---|---|---|---|"] + [f"| {c['run']} | {fm(c['removal_day'],1)} | **{c['cls']}** | {fm(c['R'])} | {c['decline']} | {c['f_end']:.3f} |" for c in cl3]
    L += ["", f"### Removal-time bracket (NEW runs, medium). Reduced boundary {P['reduced_boundary_day']:.2f} d; CFD boundary {fm(tb,2)} d; pass(±1.0 d)={res['boundary_pass']}", "| removal day | f at +1 d (CFD / pred) | class CFD / pred | match | R | f_end | loss day | mean heat on former footprint W | footprint melt m³ |", "|---|---|---|---|---|---|---|---|---|"] + [(f"| {b['removal_day']:g} | MISSING |" if b.get("status") else f"| {b['removal_day']:g} | {b['f_plus1d']:.4f} / {b['pred_f_plus1d']:.4f} | {b['cls']} / {b['pred_cls']} | {b['class_match']} | {fm(b['R'])} | {b['f_end']:.4f} | {fm(b['loss'],2)} | {b['region_heat_toe_fp_mean_W']:.1f} | {b['melt_toe_fp_m3']:.4f} |") for b in br]
    L += ["", "### Counterfactual: differential incision disabled (uniform lowering, same total melt)", "| run | removal | contact lost | min f | f0 → f_end | toe-fp melt | bypass melt | crest toe → byp (end) |", "|---|---|---|---|---|---|---|---|"] + [f"| {c['run']} | {fm(c['removal_day'],1)} | {c['contact_lost']} | {c['min_f']:.3f} | {c['f0']:.3f} → {c['f_end']:.3f} | {c['toe_fp_melt_m3']:.4f} | {c['byp_melt_m3']:.4f} | {c['crest_toe_end']:.3f} / {c['crest_byp_end']:.3f} |" for c in cf]
    open(os.path.join(out, "results_v2.md"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
if __name__ == "__main__": main(*sys.argv[1:6])
