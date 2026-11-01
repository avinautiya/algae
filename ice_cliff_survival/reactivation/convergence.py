"""Event-time / coupling convergence (development runs). `python -m ice_cliff_survival.reactivation.convergence DIR OUT.json`: separate spatial (dx), hydrodynamic CFL and bed-update-interval (dtm) series for the strip layout without removal and with removal at 8 d; also the corridor (no removal)."""
import glob, json, math, os, sys
import numpy as np
from . import metrics as M
def load(d):
    out = {}
    for f in sorted(glob.glob(os.path.join(d, "conv_*.json"))):
        r = json.load(open(f)); r["name"] = os.path.basename(f)[:-5]
        if r.get("status") == "complete": out[r["name"]] = r
    return out
def qoi(r):
    s = M.summarize(r); a = M.arrays(r); days = np.array([1.0, 2.0, 3.0, 5.0, 8.0, 10.0, 15.0]); f = np.interp(days, a["day"], a["f"])
    return dict(loss=s["loss"], reconnect=s["sustained_reconnection"], f_end=s["f_end"], delivered=s["delivered_fraction"], toe_melt_all=s["toe_melt_all_m3"], f_series=f.tolist(), crest_toe_end=s["crest_toe_end"], crest_byp_end=s["crest_byp_end"], max_water_resid=s["max_water_resid"], max_energy_resid=s["max_energy_resid"], wall_s=s["wall_s"], cells=s["cells"])
def order(a, b, c, ratio=2.0):
    try: return float(math.log(abs(a - b) / abs(b - c)) / math.log(ratio))
    except Exception: return float("nan")
def study(d):
    R = load(d); out = {}
    for tag in ("n", "r8"):
        S = {}
        for key, vals, fmt in (("dx", (0.1, 0.05, 0.025), "conv_strip_%s_dx%g"), ("dtm", (12, 6, 3), "conv_strip_%s_dtm%g"), ("cfl", (0.8, 0.4, 0.2), "conv_strip_%s_cfl%g")):
            rows = []
            for v in vals:
                nm = fmt % (tag, v) if not ((key == "dtm" and v == 6) or (key == "cfl" and v == 0.4)) else "conv_strip_%s_dx0.05" % tag
                if nm in R: rows.append(dict(value=v, name=nm, **qoi(R[nm])))
            S[key] = rows
            if len(rows) == 3:
                for q in ("loss", "reconnect", "f_end", "delivered", "toe_melt_all"):
                    v = [x[q] for x in rows]
                    if all(t is not None for t in v): S.setdefault("diffs", {}).setdefault(key, {})[q] = dict(values=v, d12=v[1] - v[0], d23=v[2] - v[1], observed_order=order(*v))
        out["strip_" + tag] = S
    out["corridor_none"] = [dict(dx=v, **qoi(R["conv_corr_n_dx%g" % v])) for v in (0.1, 0.05, 0.025) if "conv_corr_n_dx%g" % v in R]
    return out
if __name__ == "__main__":
    o = study(sys.argv[1]); json.dump(o, open(sys.argv[2], "w"), indent=1, default=float)
    for k in ("strip_n", "strip_r8"):
        print(k)
        for s, d in o[k].get("diffs", {}).items():
            for q, v in d.items(): print("  ", s, q, [round(x, 5) if x is not None else None for x in v["values"]], "order", round(v["observed_order"], 2))
    print("corridor", [(c["dx"], round(c["f_end"], 4), round(c["delivered"], 4), c["loss"]) for c in o["corridor_none"]])
