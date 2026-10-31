"""Sequence-effect metrics from per-run records (network r17/r33, resolved reference) and RTM_seq predictions.  python evaluate_seq.py dev|test -> eval_<set>.json
Per configuration/treatment/simulator: metrics for AB, BA, AVG; sequence effect = AB - BA; averaging error = AVG - (AB or BA)."""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "ice_network_eval"))
import analyze_campaign as A
import sdesign as design
TIE = 0.03


def metrics_of(r, p):
    s = A.series(r); m = r["metrics"]; t = s["t"]; tail = p["tail_s"]; i_c = np.searchsorted(t, t[-1] - tail)
    q1, q2 = s["Q1"], s["Q2"]; sh = q2 / np.maximum(q1 + q2, 1e-12)
    d1, d2 = m["delivered1_m3"], m["delivered2_m3"]
    return dict(cum1=d1, cum2=d2, share_deliv=d2 / max(d1 + d2, 1e-12), share_end=float(sh[-1]), share_cont_start=float(sh[i_c]), Q1_end=float(q1[-1]), Q2_end=float(q2[-1]),
                dominant_end="R2" if q2[-1] > q1[-1] else "R1", tie_end=bool(abs(sh[-1] - 0.5) < TIE), dominant_cum="R2" if d2 > d1 else "R1", tie_cum=bool(abs(d2 / max(d1 + d2, 1e-12) - 0.5) < TIE),
                sill1_mm=m["sill_lowering1_m"] * 1e3, sill2_mm=m["sill_lowering2_m"] * 1e3, gap_end_mm=float(1e3 * (s["sill2"][-1] - s["sill1"][-1])), heat_to_ice_J=m["heat_to_ice_J"], stage_end=float(s["stage"][-1]),
                t_dom2=m["t_dom2_s"], t_act2=m["t_act2_s"], share_series=sh[::10].tolist(), t_series=t[::10].tolist(), conservation=r.get("conservation"), solver_s=r["timing"]["solver_s"])


def load(which):
    table = {"dev": design.DEV, "test": design.TEST}[which]
    runs = {r["spec"]["label"]: r for r in A.load(lambda r: A.case_set(r["spec"]) == "sequence" and r["status"] == "completed")}
    out = {}
    for cid, (desc, c) in table.items():
        seq, mean = design.physical(c); p = dict(tail_s=c["tail_h"] * 3600.0); row = dict(description=desc, config=c, mean=mean)
        for sim in ("thermal_r17", "fixed_r17", "dissipation_r17", "thermal_r33", "thermal_ref200", "thermal_ref400", "thermalnu0.5_r17", "thermalnu2_r17", "thermalTOL_r17"):
            tr, kind = sim.rsplit("_", 1)[0], sim.rsplit("_", 1)[1]; fl = {}
            for fn in ("AB", "BA", "AVG"):
                r = runs.get(f"seq_{cid}_{fn}_{tr}_{kind}")
                if r: fl[fn] = metrics_of(r, p)
            if fl: row[sim] = fl
        pf = os.path.join(HERE, f"pred_{which}", cid + ".json")
        if os.path.exists(pf): row["rtm"] = json.load(open(pf))
        out[cid] = row
    return out


def effect(fl, key):
    return None if not ("AB" in fl and "BA" in fl) else fl["AB"][key] - fl["BA"][key]


if __name__ == "__main__":
    out = load(sys.argv[1]); json.dump(out, open(os.path.join(HERE, f"eval_{sys.argv[1]}.json"), "w"), default=float)
    print("id | net r17 AB-BA: share_end share_deliv cum2(m3) sill1 sill2(mm) heat(rel %) | dominant AB/BA/AVG | RTM_seq AB-BA share_end | fixed net share_end | AVG-AB share_end, AVG-BA")
    for cid, r in out.items():
        n = r.get("thermal_r17"); f = r.get("fixed_r17"); rt = r.get("rtm")
        if not n or "AB" not in n or "BA" not in n or "AVG" not in n: continue
        ab, ba, av = n["AB"], n["BA"], n["AVG"]
        print(cid, "| %+.3f %+.3f %+.0f %+.1f %+.1f %+.1f%% | %s/%s/%s | %s | %s | %+.3f %+.3f" % (ab["share_end"] - ba["share_end"], ab["share_deliv"] - ba["share_deliv"], ab["cum2"] - ba["cum2"], ab["sill1_mm"] - ba["sill1_mm"], ab["sill2_mm"] - ba["sill2_mm"], 100 * (ab["heat_to_ice_J"] / ba["heat_to_ice_J"] - 1), ab["dominant_end"], ba["dominant_end"], av["dominant_end"],
              None if not rt else "%+.3f" % (rt["thermal"]["AB"]["share2_end"] - rt["thermal"]["BA"]["share2_end"]), None if not f or "AB" not in f else "%+.3f" % (f["AB"]["share_end"] - f["BA"]["share_end"]), av["share_end"] - ab["share_end"], av["share_end"] - ba["share_end"]))
