"""Tables and accuracy-cost summary from eval_results.json (reduced network vs resolved reference)."""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
R = json.load(open(os.path.join(HERE, "records", "eval_results.json")))
SCALARS = ["sill_lowering1_m", "sill_lowering2_m", "delivered1_m3", "delivered2_m3", "heat_to_ice_J", "stage_end_m"]
EVENTS = ["t_act1_s", "t_act2_s", "t_dom2_s"]
out = {}
for sid, d in R.items():
    ref = d["reference_metrics"]; m200, m400 = ref["200"], ref["400"]
    ext, unc = {}, {}
    for k in SCALARS + EVENTS:
        a, b = m200.get(k), m400.get(k)
        if a is None or b is None:
            ext[k], unc[k] = (b, None); continue
        ext[k] = b + (b - a); unc[k] = abs(b - a)              # first-order Richardson; uncertainty taken as |400 - 200| (the correction size)
    rows = {}
    for cid, v in d["reduced"].items():
        mt = v["metrics"]; e = {}
        for k in SCALARS + EVENTS:
            x, r_ = mt.get(k), ext[k]
            if x is None or r_ is None:
                e[k] = dict(reduced=x, reference=r_, err=None); continue
            e[k] = dict(reduced=x, reference=r_, err=x - r_, rel=(x - r_) / abs(r_) if r_ else None, ref_unc=unc[k], resolved=bool(abs(x - r_) > (unc[k] or 0)))
        rows[cid] = dict(errors=e, wall_s=v["wall_s"], nfev=v["discharge_law_evaluations"], series=v["series_error_vs_finest_reference"], conservation=v["conservation"])
    out[sid] = dict(reference_extrapolated=ext, reference_uncertainty=unc, reference_wall_s=d["reference_wall_s"], reduced=rows, scenario=d["scenario"])
json.dump(out, open(os.path.join(HERE, "records", "eval_summary.json"), "w"), indent=1, default=float)
# console
for sid, d in out.items():
    print("==", sid, d["scenario"])
    print("  ref (extrapolated):", {k: (None if v is None else round(v, 3)) for k, v in d["reference_extrapolated"].items() if k in EVENTS + ["sill_lowering1_m", "sill_lowering2_m", "delivered1_m3", "delivered2_m3"]})
    print("  ref uncertainty   :", {k: (None if v is None else round(v, 3)) for k, v in d["reference_uncertainty"].items() if k in EVENTS + ["sill_lowering1_m", "sill_lowering2_m", "delivered1_m3", "delivered2_m3"]})
    for cid, r in d["reduced"].items():
        e = r["errors"]
        f = lambda k: "None" if e[k].get("err") is None else "%+.1f%%" % (100 * e[k]["rel"]) if abs(e[k]["reference"]) > 0 else "n/a"
        print("  %-10s wall %5.1fs nfev %6d | sill1 %s sill2 %s deliv1 %s deliv2 %s | t_act2 %s t_dom2 %s | Q1 rmse/I %.3f Q2 %.3f" % (
            cid, r["wall_s"], r["nfev"], f("sill_lowering1_m"), f("sill_lowering2_m"), f("delivered1_m3"), f("delivered2_m3"), f("t_act2_s"), f("t_dom2_s"),
            r["series"]["Q1"]["rmse_over_I"], r["series"]["Q2"]["rmse_over_I"]))
