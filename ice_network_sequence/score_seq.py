"""Score the frozen sequence-effect predictions against simulations with the frozen thresholds (FREEZE.md).  python score_seq.py -> scores_test.json, SCORES_TEST.md"""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import evaluate_seq as EV, sdesign as design
P = json.load(open(os.path.join(HERE, "TEST_PREDICTIONS.json"))); TH = P["thresholds"]; PRED = P["predictions"]
E = EV.load("test"); json.dump(E, open(os.path.join(HERE, "eval_test.json"), "w"), default=float)


def eff(fl, key): return fl["AB"][key] - fl["BA"][key] if fl and "AB" in fl and "BA" in fl else None
def heat_rel(fl): return fl["AB"]["heat_to_ice_J"] / fl["BA"]["heat_to_ice_J"] - 1 if fl and "AB" in fl and "BA" in fl and fl["BA"]["heat_to_ice_J"] > 0 else None
def inc_rel(fl): return max(abs(fl["AB"]["sill1_mm"] - fl["BA"]["sill1_mm"]) / max(0.5 * (fl["AB"]["sill1_mm"] + fl["BA"]["sill1_mm"]), 1.0), abs(fl["AB"]["sill2_mm"] - fl["BA"]["sill2_mm"]) / max(0.5 * (fl["AB"]["sill2_mm"] + fl["BA"]["sill2_mm"]), 1.0)) if fl and "AB" in fl and "BA" in fl else None


def observed(fl):
    if not fl or "AB" not in fl or "BA" not in fl: return None
    ab, ba = fl["AB"], fl["BA"]
    flip = ab["dominant_end"] != ba["dominant_end"] and not ab["tie_end"] and not ba["tie_end"]
    o = dict(share_end=ab["share_end"] - ba["share_end"], share_deliv=ab["share_deliv"] - ba["share_deliv"], cum2=ab["cum2"] - ba["cum2"], heat_rel=heat_rel(fl), incision_rel=inc_rel(fl), stage_end=ab["stage_end"] - ba["stage_end"],
             persistence=None if abs(ab["share_cont_start"] - ba["share_cont_start"]) < 1e-4 else (ab["share_end"] - ba["share_end"]) / (ab["share_cont_start"] - ba["share_cont_start"]), flip=bool(flip), switch_differs=(ab["t_dom2"] is None) != (ba["t_dom2"] is None),
             dominant=[ab["dominant_end"], ba["dominant_end"], fl["AVG"]["dominant_end"] if "AVG" in fl else None], tie=[ab["tie_end"], ba["tie_end"]])
    o["meaningful_routing"] = bool(abs(o["share_end"]) >= TH["share_end"] or abs(o["share_deliv"]) >= TH["share_deliv"] or o["flip"] or o["switch_differs"])
    o["meaningful_incision"] = bool((o["incision_rel"] or 0) >= TH["incision_rel"] or abs(o["heat_rel"] or 0) >= TH["heat_rel"])
    if "AVG" in fl:
        o["avg_error_share_end"] = [fl["AVG"]["share_end"] - ab["share_end"], fl["AVG"]["share_end"] - ba["share_end"]]; o["avg_error_share_deliv"] = [fl["AVG"]["share_deliv"] - ab["share_deliv"], fl["AVG"]["share_deliv"] - ba["share_deliv"]]
        o["avg_error_heat_rel"] = [fl["AVG"]["heat_to_ice_J"] / ab["heat_to_ice_J"] - 1, fl["AVG"]["heat_to_ice_J"] / ba["heat_to_ice_J"] - 1] if ab["heat_to_ice_J"] > 0 and ba["heat_to_ice_J"] > 0 else [0.0, 0.0]
        o["avg_adequate"] = bool(max(abs(x) for x in o["avg_error_share_end"] + o["avg_error_share_deliv"]) < TH["share_end"] and max(abs(x) for x in o["avg_error_heat_rel"]) < TH["heat_rel"])
    return o


res = {}
for cid in design.TEST:
    r = E[cid]; pr = PRED[cid]; row = dict(description=r["description"], predicted=dict(share_end=pr["rtm_effect"]["share_end"], heat_rel=pr["rtm_effect"]["heat_rel"], eta=pr["eta_analytic"], product_rule=pr["product_rule_share_effect"], meaningful_routing=pr["predicted_meaningful_routing"], meaningful_incision=pr["predicted_incision_meaningful"], avg_error_share_end=[pr["rtm_share"]["AVG"] - pr["rtm_share"]["AB"], pr["rtm_share"]["AVG"] - pr["rtm_share"]["BA"]]))
    row["net17"] = observed(r.get("thermal_r17")); row["net17_fixed"] = observed(r.get("fixed_r17")); row["net33"] = observed(r.get("thermal_r33")); row["ref200"] = observed(r.get("thermal_ref200")); row["ref400"] = observed(r.get("thermal_ref400"))
    row["closure"] = {k: observed(r.get(f"thermal{k}_r17")) for k in ("nu0.5", "nu2")}; row["tol"] = observed(r.get("thermalTOL_r17"))
    # numerical uncertainty of the effect: net r17 vs r33, ref200 vs ref400
    u = {}
    if row["net17"] and row["net33"]: u["net_share_end"] = abs(row["net17"]["share_end"] - row["net33"]["share_end"]); u["net_heat_rel"] = abs(row["net17"]["heat_rel"] - row["net33"]["heat_rel"])
    if row["ref200"] and row["ref400"]: u["ref_share_end"] = abs(row["ref200"]["share_end"] - row["ref400"]["share_end"]); u["ref_heat_rel"] = abs(row["ref200"]["heat_rel"] - row["ref400"]["heat_rel"])
    row["numerical_uncertainty"] = u; res[cid] = row
json.dump(res, open(os.path.join(HERE, "scores_test.json"), "w"), indent=1, default=float)
print("id | pred avg error | net17 avg error (AVG-AB, AVG-BA)")
for cid, r in res.items():
    o = r["net17"]
    if o and "avg_error_share_end" in o: print(cid, [round(x, 3) for x in r["predicted"]["avg_error_share_end"]], [round(x, 3) for x in o["avg_error_share_end"]], "adequate" if o["avg_adequate"] else "NOT adequate")
print("id | pred share eff | net17 | net33 | ref200 | ref400 | fixed net17 | heat rel pred / net17 / ref | meaningful routing pred/net17/ref | flips")
for cid, r in res.items():
    f = lambda o, k="share_end": "-" if (not o or o.get(k) is None) else "%+.4f" % o[k]
    ref = r["ref400"] or r["ref200"]
    print(cid, "| %+.4f | %s | %s | %s | %s | %s | %+.3f / %s / %s | %s/%s/%s | %s" % (r["predicted"]["share_end"], f(r["net17"]), f(r["net33"]), f(r["ref200"]), f(r["ref400"]), f(r["net17_fixed"]), r["predicted"]["heat_rel"], f(r["net17"], "heat_rel"), f(ref, "heat_rel"), r["predicted"]["meaningful_routing"], None if not r["net17"] else r["net17"]["meaningful_routing"], None if not ref else ref["meaningful_routing"], None if not r["net17"] else r["net17"]["flip"]))
