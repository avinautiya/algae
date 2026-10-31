"""Frozen test predictions: RTM_seq effects, closed-form heat-difference eta, product-rule share effect eta * dphi_inc, and the predicted classification under the frozen thresholds."""
import json, os, sys, hashlib
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import sdesign as design, theory_eta, rtm
TH = dict(share_end=0.05, share_deliv=0.03, tie=0.03, incision_rel=0.15, heat_rel=0.10)
out = {}
for cid, (desc, c) in design.TEST.items():
    pr = json.load(open(os.path.join(HERE, "pred_test", cid + ".json"))); p = dict(c["geom"], deficit=-0.10); pp = rtm.full(p); V = pp["area"] * (pp["sill1"] - pp["deficit"] - pp["bottom"])
    e = theory_eta.eta(p, c["pair"], c["tau_h"] * 3600, c["tail_h"] * 3600, V); t, f = pr["thermal"], pr["fixed"]
    dphi = t["AVG"]["share2_end"] - f["AVG"]["share2_end"]
    eff = dict(share_end=t["AB"]["share2_end"] - t["BA"]["share2_end"], share_deliv=t["AB"]["share2_deliv"] - t["BA"]["share2_deliv"], cum2=t["AB"]["cum2"] - t["BA"]["cum2"], sill1_mm=t["AB"]["sill1_mm"] - t["BA"]["sill1_mm"], sill2_mm=t["AB"]["sill2_mm"] - t["BA"]["sill2_mm"],
               heat_rel=(t["AB"]["heat1"] + t["AB"]["heat2"]) / (t["BA"]["heat1"] + t["BA"]["heat2"]) - 1, fixed_share_end=f["AB"]["share2_end"] - f["BA"]["share2_end"])
    dom = {k: t[k]["dominant_end"] for k in ("AB", "BA", "AVG")}; sh = {k: t[k]["share2_end"] for k in ("AB", "BA", "AVG")}
    flip = dom["AB"] != dom["BA"] and all(abs(sh[k] - 0.5) >= TH["tie"] for k in ("AB", "BA"))
    meaningful_routing = bool(abs(eff["share_end"]) >= TH["share_end"] or abs(eff["share_deliv"]) >= TH["share_deliv"] or flip)
    incis = max(abs(eff["sill1_mm"]) / max(0.5 * (t["AB"]["sill1_mm"] + t["BA"]["sill1_mm"]), 1.0), abs(eff["sill2_mm"]) / max(0.5 * (t["AB"]["sill2_mm"] + t["BA"]["sill2_mm"]), 1.0))
    out[cid] = dict(description=desc, config=c, rtm_effect=eff, rtm_dominant=dom, rtm_share=sh, eta_analytic=e["eta"], dphi_inc=dphi, product_rule_share_effect=e["eta"] * dphi, predicted_flip=bool(flip), predicted_meaningful_routing=meaningful_routing,
                    predicted_incision_effect_rel=incis, predicted_incision_meaningful=bool(incis >= TH["incision_rel"] or abs(e["eta"]) >= TH["heat_rel"]), predicted_cum=dict(AB=[t["AB"]["cum1"], t["AB"]["cum2"]], BA=[t["BA"]["cum1"], t["BA"]["cum2"]], AVG=[t["AVG"]["cum1"], t["AVG"]["cum2"]]))
json.dump(dict(thresholds=TH, predictions=out), open(os.path.join(HERE, "TEST_PREDICTIONS.json"), "w"), indent=1, default=float)
L = ["| id | configuration | RTM share effect (AB-BA) | product rule eta*dphi | eta | heat rel | incision rel | dominant AB/BA/AVG | meaningful routing? | meaningful incision? |", "|---|---|---|---|---|---|---|---|---|---|"]
for cid, r in out.items(): L.append(f"| {cid} | {r['description']} | {r['rtm_effect']['share_end']:+.4f} | {r['product_rule_share_effect']:+.4f} | {r['eta_analytic']:+.3f} | {r['rtm_effect']['heat_rel']:+.3f} | {r['predicted_incision_effect_rel']:.3f} | {'/'.join(r['rtm_dominant'].values())} | {r['predicted_meaningful_routing']} | {r['predicted_incision_meaningful']} |")
open(os.path.join(HERE, "TEST_PREDICTIONS.md"), "w").write("\n".join(L)); print("\n".join(L))
