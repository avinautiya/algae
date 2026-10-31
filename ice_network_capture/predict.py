"""RTM predictions for a configuration table -> predictions_<set>.json. Writes one atomic file per configuration (resumable).  python predict.py dev|test"""
import json, os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import rtm, design


def predict(cid, p):
    q = rtm.full(dict(design.H12, **p)); out = dict(id=cid, params=q, t0=time.time())
    for mode, st in (("thermal", 17), ("dissipation", 17), ("fixed", 17), ("thermal", 33)):
        t0 = time.time(); r = rtm.rtm_run(q, mode=mode, stations=st); ev = rtm.events(r, q["I"])
        # fractional capacity growth rates r_i = lam_i E_i at the end of the run's activation (t = 4 h) and at start of both-active phase
        i = int(4 * 3600 / 60); tabs = r["tabs"]; rr = []
        for k, tab in ((1, tabs[0]), (2, tabs[1])):
            u = (r["H"][i] - tab["sill"] + (r["d1"][i] if k == 1 else r["d2"][i])); g = np.interp(u, tab["u"], tab["g"]); Qk = np.interp(u, tab["u"], tab["Q"]); E = np.interp(u, tab["u"], tab["E"])
            rr.append(dict(u=float(u), Q=float(Qk), g=float(g), lam=float(g / max(Qk, 1e-12)), E=float(E), r=float(g / max(Qk, 1e-12) * E)))
        out[f"{mode}_{st}"] = dict(events=ev, rates_at_4h=rr, wall_s=time.time() - t0, success=r["success"],
                                   series=dict(t=r["t"][::10].tolist(), Q1=r["Q1"][::10].tolist(), Q2=r["Q2"][::10].tolist(), d1=r["d1"][::10].tolist(), d2=r["d2"][::10].tolist()))
    return out


if __name__ == "__main__":
    which = sys.argv[1]; table = {"dev": design.DEV, "test": getattr(design, "TEST", {})}[which]
    os.makedirs(os.path.join(HERE, f"pred_{which}"), exist_ok=True)
    for cid, (_, p) in table.items():
        f = os.path.join(HERE, f"pred_{which}", cid + ".json")
        if os.path.exists(f): continue
        o = predict(cid, p); json.dump(o, open(f + ".tmp", "w")); os.replace(f + ".tmp", f); print(cid, "done", round(time.time() - o["t0"]), "s", flush=True)
