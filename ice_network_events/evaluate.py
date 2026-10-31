"""Event metrics for network runs (campaign records) and RTM methods.  python evaluate.py dev|test -> eval_<set>.json
Truth = network time-resolved run (TR, r17).  Per config and method: wet-hours of outlet 2, wet/dry class, delivered volume of outlet 2 (and its share of the total), final share; errors against truth."""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "ice_network_eval"))
import analyze_campaign as A, design, sim, methods as M, series as S


def net_metrics(r, qbar):
    s = A.series(r); m = sim.event_metrics(s["t"], s["Q1"], s["Q2"], qbar); m["solver_s"] = r["timing"]["solver_s"]; m["sill1_mm"] = r["metrics"]["sill_lowering1_m"] * 1e3; m["sill2_mm"] = r["metrics"]["sill_lowering2_m"] * 1e3; return m


def run(which):
    c, dev, test = design.selection(); tab = dev if which == "dev" else test
    runs = {r["spec"]["label"]: r for r in A.load(lambda r: A.case_set(r["spec"]) == "events" and r["status"] == "completed")}
    out = {}
    for cid, cfg in tab.items():
        p = dict(cfg["geom"], deficit=-0.10); ser = S.make_series(cfg["forcing"], design.DUR_H); mm = S.means(ser)
        rt, th, ft = sim.all_methods(p, ser); row = dict(category=cfg["category"], near=cfg["near"], crit=cfg["crit"], se_meta=cfg["se"], rtm=rt, net={}, mean=mm)
        for sfx, key in (("r17", "net"), ("ref200", "ref")):
            d = {}
            for m in design.METHODS:
                r = runs.get(f"evt_{cid}_{m}_{sfx}")
                if r: d[m] = net_metrics(r, mm["Q_bar"])
            if d: row[key] = d
        for tag, val in (("nu0.5", "0.5"), ("nu2", "2.0")):
            d = {m: net_metrics(runs[f"evt_{cid}_{m}_nu_factor{val}_r17"], mm["Q_bar"]) for m in ("TR", "AVG", "EA") if f"evt_{cid}_{m}_nu_factor{val}_r17" in runs}
            if d: row["net_" + tag] = d
        r = runs.get(f"evt_{cid}_TR__r17")
        if r: row["net_legacytol"] = net_metrics(r, mm["Q_bar"])
        out[cid] = row
    json.dump(out, open(os.path.join(HERE, f"eval_{which}.json"), "w"), indent=1, default=float); return out


def table(out, truth_key="net"):
    print("%-4s %-13s %5s %5s | %-22s | wet_h: TR  SE  AVG PW6 PWn EA | d2 vol: TR  AVG PW6 PWn EA | class TR/AVG/EA" % ("id", "category", "ratio", "R", "net TR wet/cls"))
    for cid, r in out.items():
        n = r.get(truth_key)
        if not n or "TR" not in n: continue
        rt = r["rtm"]; w = lambda d, m: d[m]["wet_h"] if m in d else float("nan"); v = lambda d, m: d[m]["cum2"] if m in d else float("nan")
        print("%-4s %-13s %5.2f %5.2f | %5.1f %-12s | %s | %s | %s" % (cid, r["category"] + ("*" if r["near"] else ""), r["crit"]["ratio_mean"], r["crit"]["R"], n["TR"]["wet_h"], n["TR"]["cls"],
              " ".join("%4.1f" % w(n, m) for m in ("TR", "AVG", "PW6", "PWn", "EA")), " ".join("%5.0f" % v(n, m) for m in ("TR", "AVG", "PW6", "PWn", "EA")), "/".join(n[m]["cls"][:3] if m in n else "-" for m in ("TR", "AVG", "EA"))))
        print("      rtm wet_h:", {m: round(rt[m]["wet_h"], 1) for m in rt}, "net-vs-rtm TR cum2 %.0f vs %.0f" % (n["TR"]["cum2"], rt["TR"]["cum2"]))


if __name__ == "__main__":
    o = run(sys.argv[1]); table(o)
