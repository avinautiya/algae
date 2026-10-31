"""Compare frozen RTM predictions with simulator outputs (network r17 for all treatments; resolved reference where run) from per-run records. python evaluate.py dev|test"""
import glob, json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "ice_network_eval"))
import analyze_campaign as A
import design
HORIZON = 12 * 3600.0


def gap_series(ser):   # sill gap s2 - s1 (m) on the stored grid
    return ser["sill2"] - ser["sill1"]


def sim_summary(r):
    s = A.series(r); m = r["metrics"]; t = s["t"]; I = r["spec"]["params"]["I"]
    tot1, tot2 = m["delivered1_m3"], m["delivered2_m3"]
    # gap drift over the second half of the run (after both routes are normally active)
    i0 = np.searchsorted(t, 0.5 * HORIZON); gap = gap_series(s)
    return dict(t_dom2=m["t_dom2_s"], t_act2=m["t_act2_s"], phi2_deliv=tot2 / max(tot1 + tot2, 1e-12), d1=m["sill_lowering1_m"], d2=m["sill_lowering2_m"], gap0=float(gap[0]), gap_end=float(gap[-1]), gap_drift_mm=float(1e3 * (gap[-1] - gap[i0])),
                active_end=m["active_at_end"], delivered=[tot1, tot2], series=s)


def pred_summary(pr, key):
    e = pr[key]["events"]; tot = e["delivered1_m3"] + e["delivered2_m3"]
    return dict(t_dom2=e["t_dom2_s"], t_act2=e["t_act2_s"], phi2_deliv=e["delivered2_m3"] / max(tot, 1e-12), d1=e["sill1_m"], d2=e["sill2_m"], gap_drift_mm=e["gap_drift_mm"])


def compare(which):
    table = {"dev": design.DEV, "test": getattr(design, "TEST", {})}[which]; rows = []
    runs = A.load(lambda r: A.case_set(r["spec"]) == "capture" and r["status"] == "completed")
    idx = {r["spec"]["label"]: r for r in runs}
    for cid in table:
        pf = os.path.join(HERE, f"pred_{which}", cid + ".json")
        if not os.path.exists(pf): continue
        pr = json.load(open(pf)); row = dict(id=cid, params=pr["params"])
        for tr in ("thermal", "dissipation", "fixed"):
            r = idx.get(f"cap_{cid}_{tr}_r17")
            if not r: continue
            sim = sim_summary(r); prd = pred_summary(pr, f"{tr}_17" if tr != "thermal" else "thermal_17")
            row[tr] = dict(sim={k: v for k, v in sim.items() if k != "series"}, pred=prd)
        for n in (200,):
            for tr in ("thermal", "fixed"):
                r = idx.get(f"cap_{cid}_{tr}_ref{n}")
                if r: sim = sim_summary(r); row[f"{tr}_ref{n}"] = dict(sim={k: v for k, v in sim.items() if k != "series"}, pred=pred_summary(pr, f"{tr}_33" if tr == "thermal" else "fixed_17"))
        rows.append(row)
    return rows


def occ(x): return "none" if x is None else "yes"


if __name__ == "__main__":
    rows = compare(sys.argv[1]); json.dump(rows, open(os.path.join(HERE, f"eval_{sys.argv[1]}.json"), "w"), indent=1, default=float)
    print("id | treat | t_dom2 pred/sim (s) | phi2 pred/sim | d1,d2 pred/sim (mm) | gap drift sim (mm)")
    for r in rows:
        for tr in ("thermal", "dissipation", "fixed"):
            if tr in r:
                p, s = r[tr]["pred"], r[tr]["sim"]
                print(r["id"], tr[:4], "| %s / %s" % (None if p["t_dom2"] is None else round(p["t_dom2"]), None if s["t_dom2"] is None else round(s["t_dom2"])), "| %.2f / %.2f" % (p["phi2_deliv"], s["phi2_deliv"]), "| %.1f,%.1f / %.1f,%.1f" % (p["d1"] * 1e3, p["d2"] * 1e3, s["d1"] * 1e3, s["d2"] * 1e3), "| %.1f" % s["gap_drift_mm"])
