"""Score the frozen test predictions against simulations with the frozen acceptance rules (TEST_FREEZE.md). python score_test.py -> scores_test.json, SCORES_TEST.md"""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "ice_network_eval"))
import analyze_campaign as A
import evaluate as EV, design
HZ = 12 * 3600.0
P = json.load(open(os.path.join(HERE, "TEST_PREDICTIONS.json")))
runs = {r["spec"]["label"]: r for r in A.load(lambda r: A.case_set(r["spec"]) == "capture")}


def sim(label):
    r = runs.get(label)
    if r is None: return None
    if r["status"] != "completed": return dict(status=r["status"], error=(r.get("error") or "")[:160])
    s = EV.sim_summary(r); s.pop("series"); s["status"] = "completed"; return s


def score(pred, s):
    """pred: dict(switch, t_dom2_s, phi2_end, sill1_mm, sill2_mm, gap_drift_mm); s: simulated summary"""
    out = {}
    if s is None or s.get("status") != "completed": return dict(note="simulation missing/failed")
    nondraining = s["delivered"][0] + s["delivered"][1] <= 0
    sw_s = s["t_dom2"] is not None; sw_p = pred["switch"]
    near = (pred["t_dom2_s"] is not None and pred["t_dom2_s"] > 0.85 * HZ) or (s["t_dom2"] is not None and s["t_dom2"] > 0.85 * HZ) or abs(s["phi2_deliv"] - 0.5) < 0.03
    out["nondraining"] = bool(nondraining); out["switch_pred"] = sw_p; out["switch_sim"] = sw_s; out["near_tie_switch"] = bool(near)
    out["switch_ok"] = None if (near or nondraining) else (sw_p == sw_s)
    out["t_dom2_rel_err"] = None if (not sw_p or not sw_s or near) else (pred["t_dom2_s"] - s["t_dom2"]) / s["t_dom2"]
    out["t_dom2_ok"] = None if out["t_dom2_rel_err"] is None else abs(out["t_dom2_rel_err"]) <= 0.10
    d_s, d_p = s["gap_drift_mm"], pred["gap_drift_mm"]
    out["drift_sim_mm"], out["drift_pred_mm"] = d_s, d_p
    out["drift_ok"] = None if (abs(d_s) < 2 or abs(d_p) < 2) else (np.sign(d_s) == np.sign(d_p))
    out["phi2_err"] = pred["phi2_end"] - s["phi2_deliv"]; out["phi2_ok"] = abs(out["phi2_err"]) <= 0.03
    e1 = (pred["sill1_mm"] - s["d1"] * 1e3) / max(abs(s["d1"] * 1e3), 1.0); e2 = (pred["sill2_mm"] - s["d2"] * 1e3) / max(abs(s["d2"] * 1e3), 1.0)
    out["sill_rel_err"] = [e1, e2]; out["sill_ok"] = bool(abs(e1) <= 0.15 and abs(e2) <= 0.15)
    return out


res = {}
for cid in design.TEST:
    p = P[cid]; row = dict(description=p["description"], gap0=p["gap0"], qs_threshold=p["qs_capture_threshold_gap"], qs_kappa=p["qs_kappa"])
    row["sim"] = {k: sim(f"cap_{cid}_{k}") for k in ("thermal_r17", "dissipation_r17", "fixed_r17", "thermal_r33", "thermal_ref200", "thermal_ref400", "fixed_ref200", "fixed_ref400")}
    th = p["rtm_thermal"]
    row["rtm17_vs_net17"] = score(th, row["sim"]["thermal_r17"])
    if row["sim"]["thermal_ref200"] or row["sim"]["thermal_ref400"]:
        tgt = row["sim"]["thermal_ref400"] if row["sim"]["thermal_ref400"] and row["sim"]["thermal_ref400"].get("status") == "completed" else row["sim"]["thermal_ref200"]
        row["target_ref"] = "ref400" if tgt is row["sim"]["thermal_ref400"] else "ref200"
        row["rtm33_vs_ref"] = score(p["rtm_thermal33"], tgt); row["rtm17_vs_ref"] = score(th, tgt); row["net17_vs_ref"] = score(dict(switch=row["sim"]["thermal_r17"]["t_dom2"] is not None, t_dom2_s=row["sim"]["thermal_r17"]["t_dom2"], phi2_end=row["sim"]["thermal_r17"]["phi2_deliv"], sill1_mm=row["sim"]["thermal_r17"]["d1"] * 1e3, sill2_mm=row["sim"]["thermal_r17"]["d2"] * 1e3, gap_drift_mm=row["sim"]["thermal_r17"]["gap_drift_mm"]), tgt) if row["sim"]["thermal_r17"].get("status") == "completed" else None
        r33 = row["sim"]["thermal_r33"]
        if r33 and r33.get("status") == "completed": row["net33_vs_ref"] = score(dict(switch=r33["t_dom2"] is not None, t_dom2_s=r33["t_dom2"], phi2_end=r33["phi2_deliv"], sill1_mm=r33["d1"] * 1e3, sill2_mm=r33["d2"] * 1e3, gap_drift_mm=r33["gap_drift_mm"]), tgt)
        a, b = row["sim"]["thermal_ref200"], row["sim"]["thermal_ref400"]
        if a and b and a.get("status") == "completed" and b.get("status") == "completed":
            row["ref_200_vs_400"] = dict(phi2=b["phi2_deliv"] - a["phi2_deliv"], sill1_mm=(b["d1"] - a["d1"]) * 1e3, sill2_mm=(b["d2"] - a["d2"]) * 1e3, drift_mm=b["gap_drift_mm"] - a["gap_drift_mm"], t_dom2=None if a["t_dom2"] is None or b["t_dom2"] is None else b["t_dom2"] - a["t_dom2"], switch_200=a["t_dom2"] is not None, switch_400=b["t_dom2"] is not None)
    # comparators against the thermal network r17 (all cases) : fixed geometry and dissipation-only incision
    t17 = row["sim"]["thermal_r17"]
    for name, key in (("fixed_geometry", "fixed_r17"), ("dissipation_only", "dissipation_r17")):
        c = row["sim"][key]
        if t17 and t17.get("status") == "completed" and c and c.get("status") == "completed":
            row[f"{name}_vs_thermal_net"] = dict(switch_agree=(c["t_dom2"] is not None) == (t17["t_dom2"] is not None), drift_sim=c["gap_drift_mm"], drift_thermal=t17["gap_drift_mm"], drift_sign_agree=bool((abs(c["gap_drift_mm"]) < 2 and abs(t17["gap_drift_mm"]) < 2) or (abs(c["gap_drift_mm"]) >= 2 and abs(t17["gap_drift_mm"]) >= 2 and np.sign(c["gap_drift_mm"]) == np.sign(t17["gap_drift_mm"]))),
                                                 phi2_diff=c["phi2_deliv"] - t17["phi2_deliv"], sill_rel_err=[(c["d1"] - t17["d1"]) / max(t17["d1"], 1e-6), (c["d2"] - t17["d2"]) / max(t17["d2"], 1e-6)])
    res[cid] = row
json.dump(res, open(os.path.join(HERE, "scores_test.json"), "w"), indent=1, default=float)
def tally(key, fields):
    out = {}
    for f in fields:
        v = [r[key][f] for r in res.values() if r.get(key) and r[key].get(f) is not None]; out[f] = (int(sum(bool(x) for x in v)), len(v))
    return out
summary = dict(rtm17_vs_net17=tally("rtm17_vs_net17", ("switch_ok", "t_dom2_ok", "drift_ok", "phi2_ok", "sill_ok")), rtm33_vs_ref=tally("rtm33_vs_ref", ("switch_ok", "t_dom2_ok", "drift_ok", "phi2_ok", "sill_ok")), rtm17_vs_ref=tally("rtm17_vs_ref", ("switch_ok", "t_dom2_ok", "drift_ok", "phi2_ok", "sill_ok")),
               net17_vs_ref=tally("net17_vs_ref", ("switch_ok", "t_dom2_ok", "drift_ok", "phi2_ok", "sill_ok")), net33_vs_ref=tally("net33_vs_ref", ("switch_ok", "t_dom2_ok", "drift_ok", "phi2_ok", "sill_ok")),
               fixed_geometry_vs_thermal_net=dict(switch_agree=(sum(r["fixed_geometry_vs_thermal_net"]["switch_agree"] for r in res.values() if "fixed_geometry_vs_thermal_net" in r), sum("fixed_geometry_vs_thermal_net" in r for r in res.values())), drift_sign_agree=(sum(r["fixed_geometry_vs_thermal_net"]["drift_sign_agree"] for r in res.values() if "fixed_geometry_vs_thermal_net" in r), sum("fixed_geometry_vs_thermal_net" in r for r in res.values()))),
               dissipation_only_vs_thermal_net=dict(switch_agree=(sum(r["dissipation_only_vs_thermal_net"]["switch_agree"] for r in res.values() if "dissipation_only_vs_thermal_net" in r), sum("dissipation_only_vs_thermal_net" in r for r in res.values())), drift_sign_agree=(sum(r["dissipation_only_vs_thermal_net"]["drift_sign_agree"] for r in res.values() if "dissipation_only_vs_thermal_net" in r), sum("dissipation_only_vs_thermal_net" in r for r in res.values()))))
json.dump(summary, open(os.path.join(HERE, "scores_test_summary.json"), "w"), indent=1); print(json.dumps(summary, indent=1))
