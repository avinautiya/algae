"""RTM_seq predictions (AB, BA, AVG; thermal and fixed geometry) -> pred_<set>/<id>.json (atomic, resumable).  python predict_seq.py dev|test"""
import json, os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import rtm_seq, sdesign as design, rtm


def summarize(r, p):
    t = r["t"]; Q1, Q2 = r["Q1"], r["Q2"]; I = Q1 + Q2
    i_sw = np.flatnonzero(Q2 - Q1 > 0)
    return dict(Q1_end=float(Q1[-1]), Q2_end=float(Q2[-1]), share2_end=float(Q2[-1] / max(Q1[-1] + Q2[-1], 1e-12)), cum1=float(r["cum1"][-1]), cum2=float(r["cum2"][-1]), share2_deliv=float(r["cum2"][-1] / max(r["cum1"][-1] + r["cum2"][-1], 1e-12)),
                sill1_mm=float(r["s1"][-1] * 1e3), sill2_mm=float(r["s2"][-1] * 1e3), gap_end_mm=float(1e3 * (p["offset"] - (r["s2"][-1] - r["s1"][-1]))), TL_end=float(r["TL"][-1]), V_end=float(r["V"][-1]),
                heat1=float(r["heat1"][-1]), heat2=float(r["heat2"][-1]), dominant_end="R2" if Q2[-1] > Q1[-1] else "R1", dominant_cum="R2" if r["cum2"][-1] > r["cum1"][-1] else "R1",
                first_R2_over_R1_s=float(t[i_sw[0]]) if i_sw.size else None, active2_end=bool(Q2[-1] > 0.01 * I[-1]))


def predict(cid, c):
    seq, mean = design.physical(c); p = rtm.full(dict(c["geom"], deficit=-0.10)); out = dict(id=cid, mean=mean, params=p, t0=time.time())
    q1 = dict(p, T=1.0)
    for mode in ("thermal", "fixed"):
        tabs = (rtm.tabulate(q1, 1, 17, "fixed" if mode == "fixed" else "thermal"), rtm.tabulate(q1, 2, 17, "fixed" if mode == "fixed" else "thermal"))
        out[mode] = {fn: summarize(rtm_seq.run(p, s, mode=mode, tabs=tabs, lake_T0=mean["T_bar"]), p) for fn, s in seq.items()}
    return out


if __name__ == "__main__":
    which = sys.argv[1]; table = {"dev": design.DEV, "test": design.TEST}[which]; os.makedirs(os.path.join(HERE, f"pred_{which}"), exist_ok=True)
    for cid, (_, c) in table.items():
        f = os.path.join(HERE, f"pred_{which}", cid + ".json")
        if os.path.exists(f): continue
        o = predict(cid, c); json.dump(o, open(f + ".tmp", "w")); os.replace(f + ".tmp", f); print(cid, "done", round(time.time() - o["t0"]), "s", flush=True)
