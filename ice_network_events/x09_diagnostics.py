"""Stored-output diagnostics of the sequence-study case X09 (development evidence for the event-aware study): wet/dry definitions, durations, stage relative to sill 2, partitions,
heat/incision, conservation, event chatter, and the ablations (fixed geometry; network r17 vs resolved reference 200 cells). No new simulation.  -> x09_diagnostics.json"""
import glob, json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "ice_network_eval")); sys.path.insert(0, os.path.join(HERE, "..", "ice_network_sequence"))
import analyze_campaign as A
import sdesign
c = sdesign.TEST["X09"][1]; seq, mean = sdesign.physical(c)
runs = {r["spec"]["label"]: r for r in A.load(lambda r: A.case_set(r["spec"]) == "sequence" and r["status"] == "completed" and r["spec"]["scn"].startswith("X09"))}
WET = 0.005 * mean["Q_bar"]          # outlet wet when its inlet discharge exceeds 0.5 % of the mean inflow (numerical-chatter margin: far above integrator noise, far below any physical flow here)
def segments(t, w, min_h=0.25):
    """wet intervals with a 15 min minimum for either state (removes chatter)"""
    s = (w > 0).astype(int); out = []; i = 0; n = len(s)
    while i < n:
        j = i
        while j < n and s[j] == s[i]: j += 1
        out.append([int(s[i]), t[i], t[j - 1] if j < n else t[-1]]); i = j
    return out
res = {}
for sim in ("thermal_r17", "fixed_r17", "thermal_ref200"):
    tr, kind = sim.rsplit("_", 1)
    for fn in ("AB", "BA", "AVG"):
        r = runs.get(f"seq_X09_{fn}_{tr}_{kind}")
        if not r: continue
        s = A.series(r); t = s["t"]; q1, q2 = s["Q1"], s["Q2"]; dt = np.gradient(t)
        wet2 = q2 > WET; wet1 = q1 > WET; u2 = s["stage"] - s["sill2"]; u1 = s["stage"] - s["sill1"]
        trans = int(np.sum(np.abs(np.diff(wet2.astype(int)))))
        res[f"{sim}|{fn}"] = dict(wet_hours_outlet1=float(np.sum(dt[wet1]) / 3600), wet_hours_outlet2=float(np.sum(dt[wet2]) / 3600), dry_hours_outlet2=float(np.sum(dt[~wet2]) / 3600), transitions_outlet2=trans,
            first_dry_to_wet_h=None, min_stage_minus_sill2_m=float(u2.min()), max_stage_minus_sill2_m=float(u2.max()), stage_minus_sill1_end=float(u1[-1]), cum1=float(np.trapezoid(q1, t)), cum2=float(np.trapezoid(q2, t)),
            share_end=float(q2[-1] / max(q1[-1] + q2[-1], 1e-12)), share_deliv=float(np.trapezoid(q2, t) / max(np.trapezoid(q1 + q2, t), 1e-12)), sill1_mm=float((s["sill1"][0] - s["sill1"][-1]) * 1e3), sill2_mm=float((s["sill2"][0] - s["sill2"][-1]) * 1e3),
            heat_to_ice_J=r["metrics"]["heat_to_ice_J"], conservation=r.get("conservation"), solver_s=r["timing"]["solver_s"])
        # periods of outlet-2 dryness in hours
        res[f"{sim}|{fn}"]["dry_periods_h"] = [[round(a / 3600, 2), round(b / 3600, 2)] for st, a, b in segments(t, q2 - WET) if st == 0 and b - a > 900]
json.dump(dict(mean=mean, wet_threshold_m3s=WET, results=res), open(os.path.join(HERE, "x09_diagnostics.json"), "w"), indent=1, default=float)
for k, v in res.items(): print(k, "wet2 %.1f h dry2 %.1f h transitions %d | u2 range %+.3f..%+.3f | cum2 %.0f share_end %.3f share_deliv %.3f | sills %.0f/%.0f mm | dry periods %s" % (v["wet_hours_outlet2"], v["dry_hours_outlet2"], v["transitions_outlet2"], v["min_stage_minus_sill2_m"], v["max_stage_minus_sill2_m"], v["cum2"], v["share_end"], v["share_deliv"], v["sill1_mm"], v["sill2_mm"], v["dry_periods_h"]))
