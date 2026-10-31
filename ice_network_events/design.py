"""Configuration draw and campaign-case builder for the event study.  python design.py draw   -> candidates.json (RTM-only screening, no network run)
                                                                    python design.py cases dev|test -> cases_<set>.json (campaign specs)
Configurations: two-outlet lake (initial stage 0.10 m above sill 1), 36 h hourly forcing series (diurnal x spell).  Seeded draw (seed 71); categories from the criterion and the water-storage model SE only."""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "ice_network_capture")); sys.path.insert(0, os.path.join(HERE, "..", "ice_network_sequence"))
import series as S, methods as M, sim, rtm
DUR_H = 36.0
METHODS = ("TR", "AVG", "PW6", "PW3", "EA", "PWn")      # network runs; SE is RTM-only


def draw(n, rng):
    out = []
    for _ in range(n):
        geom = dict(offset=float(rng.uniform(0.05, 0.30)), area=float(np.exp(rng.uniform(np.log(500), np.log(4000)))), L2=float(rng.choice([60.0, 100.0, 120.0])))
        f = dict(Qm=float(np.exp(rng.uniform(np.log(0.008), np.log(0.10)))), amp=float(rng.uniform(0.2, 0.9)), phase_h=float(rng.uniform(0, 24)), depth=float(rng.choice([0.0, rng.uniform(0.3, 0.9)])),
                 spell_start_h=float(rng.uniform(6, 18)), spell_len_h=float(rng.uniform(5, 14)), Tm=float(np.exp(rng.uniform(np.log(0.3), np.log(1.2)))), gamma=0.5)
        out.append(dict(geom=geom, forcing=f))
    return out


def classify(c):
    p = dict(c["geom"], deficit=-0.10); ser = S.make_series(c["forcing"], DUR_H); th = M.thresholds(p); crit = M.criterion(ser, th, p); ftabs = sim.fixed_tabs(p); r, se = sim.se_run(p, ser, ftabs)
    return dict(crit={k: (v if not isinstance(v, (np.floating, np.bool_)) else v.item()) for k, v in crit.items()}, se=se, near=bool(abs(crit["ratio_mean"] - 1) < 0.25), category=se["cls"])


def pick(cands, quota):
    """quota: list of (category, near, count); first candidates (in draw order) meeting each requirement"""
    chosen = []; used = set()
    for cat, near, k in quota:
        for i, c in enumerate(cands):
            if i in used or c["category"] != cat or c["near"] != near: continue
            chosen.append(i); used.add(i); k -= 1
            if k == 0: break
        if k > 0: print("WARNING quota unmet", cat, near, k)
    return chosen


def spec(cid, c, method, kind="reduced", ncell=200, stations=17, extra=None, wall=None):
    p = rtm.full(dict(c["geom"], deficit=-0.10)); ser = S.make_series(c["forcing"], DUR_H); m = S.means(ser); dur = m["duration"]
    p["duration_h"] = dur / 3600.0
    patch = dict(lake=dict(area_m2=p["area"], initial_stage_below_sill1_m=p["deficit"]), duration_h=p["duration_h"], channels=[dict(width_m=p["w1"], length_m=p["L1"], slope=p["S1"], crest_length_m=p["crest"]), dict(width_m=p["w2"], length_m=p["L2"], slope=p["S2"], crest_length_m=p["crest"])])
    scn = dict(id=f"{cid}{method}", sill_offset_m=p["offset"], inflow_m3_s=m["Q_bar"], inflow_T_C=m["T_bar"], lake_T0=ser["H"][0] / ser["Q"][0], ramp_s=900.0)
    if method == "TR": scn["series"] = ser
    else:
        th = M.thresholds(dict(c["geom"], deficit=-0.10)); ftabs = sim.fixed_tabs(dict(c["geom"], deficit=-0.10)); r, _ = sim.se_run(dict(c["geom"], deficit=-0.10), ser, ftabs)
        scn["segments"] = M.forcing_for(method, ser, th, dur, dict(t=r["t"], wet=r["Q2"] > 0.005 * m["Q_bar"]))
    ov = dict(n_man=p["n_man"], nu_factor=p["nu_factor"], atol_energy=1e-3); ov.update(extra or {})
    base = dict(kind=kind, scn=scn["id"], scn_override=scn, cfg_patch=patch, overrides=ov, params=dict(p, cid=cid, method=method, mean=m), version="events_v1")
    tag = "" if not extra else "_" + "_".join(f"{k}{v}" for k, v in extra.items() if k != "atol_energy")
    if kind == "reduced": o = dict(base, stations=stations, dt_geom=450 if stations == 17 else 225, label=f"evt_{cid}_{method}{tag}_r{stations}")
    else: o = dict(base, ncell=ncell, label=f"evt_{cid}_{method}{tag}_ref{ncell}")
    if wall: o["wall_limit_s"] = wall
    return o


if __name__ == "__main__":
    if sys.argv[1] == "draw":
        rng = np.random.default_rng(71); raw = draw(200, rng); cands = []
        for i, c in enumerate(raw):
            r = classify(c); cands.append(dict(idx=i, **c, **r)); 
        json.dump(cands, open(os.path.join(HERE, "candidates.json"), "w"), indent=1, default=float)
        import collections; print(collections.Counter((c["category"], c["near"]) for c in cands))


def selection():
    c = json.load(open(os.path.join(HERE, "candidates.json")))
    dev = pick(c[:100], [("wet", False, 1), ("intermittent", True, 2), ("intermittent", False, 2), ("dry", False, 1), ("dry", True, 1)])
    test = [100 + i for i in pick(c[100:], [("wet", False, 2), ("intermittent", True, 3), ("intermittent", False, 3), ("dry", False, 1), ("dry", True, 1)])]
    return c, {f"D{k+1:02d}": c[i] for k, i in enumerate(dev)}, {f"E{k+1:02d}": c[i] for k, i in enumerate(test)}


def build(which):
    c, dev, test = selection(); tab = dev if which == "dev" else test; out = []
    for cid, cfg in tab.items():
        for m in METHODS: out.append(spec(cid, cfg, m))
    if which == "test":
        mid = [k for k, v in test.items() if v["category"] == "intermittent"][:3]
        for cid in mid:
            for tag, ov in (("nu0.5", dict(nu_factor=0.5)), ("nu2", dict(nu_factor=2.0))):
                for m in ("TR", "AVG", "EA"): out.append(spec(cid, test[cid], m, extra=ov))
        for cid in mid[:2]: out.append(spec(cid, test[cid], "TR", extra=dict(atol_energy=None)))
        ref = [k for k, v in test.items() if v["category"] == "intermittent" and v["near"]][:2] + [k for k, v in test.items() if v["category"] == "dry" and v["near"]][:1] + [k for k, v in test.items() if v["category"] == "intermittent" and not v["near"]][:1]
        for cid in ref:
            for m in ("TR", "AVG", "EA"): out.append(spec(cid, test[cid], m, kind="reference", ncell=200, wall=2400))
    return out, tab


if __name__ == "__main__" and sys.argv[1] in ("dev", "test"):
    out, tab = build(sys.argv[1]); json.dump(out, open(os.path.join(HERE, f"cases_{sys.argv[1]}.json"), "w"), indent=1, default=float)
    json.dump({k: dict(idx=v["idx"], category=v["category"], near=v["near"], geom=v["geom"], forcing=v["forcing"], crit=v["crit"], se=v["se"]) for k, v in tab.items()}, open(os.path.join(HERE, f"configs_{sys.argv[1]}.json"), "w"), indent=1, default=float)
    print(len(out), "cases"); [print(k, v["category"], v["near"], round(v["crit"]["ratio_mean"], 2), "R=%.2f" % v["crit"]["R"]) for k, v in tab.items()]
