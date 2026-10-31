"""Local exponents of the tabulated routes: gamma = dlnQ/dlnu (rating elasticity), alpha = dlnE/dlnQ (incision sensitivity), index = dln(lambda E)/dlnQ (replicator stability index, <0 stabilising at fixed lake level) at Q = 0.5 I.
python indices.py -> rtm_index.json (development configurations + the symmetric weir-controlled test geometry; thermal and dissipation-only)"""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import rtm, design
out = {}
cfgs = dict(design.DEV); cfgs["T14"] = design.TEST["T14"]; cfgs["T13"] = design.TEST["T13"]
for cid, (desc, p) in cfgs.items():
    q = rtm.full(dict(design.H12, **p)); row = dict(description=desc)
    for mode in ("thermal", "dissipation"):
        for k in (1, 2):
            t = rtm.tabulate(q, k, 17, mode); m = (t["Q"] > 1e-6) & (t["E"] > 0); Q = t["Q"][m]; u = t["u"][m]; E = t["E"][m]; lam = t["g"][m] / Q
            i = int(np.argmin(abs(Q - 0.5 * q["I"]))); lq = np.log(Q)
            row[f"{mode}_route{k}"] = dict(Q=float(Q[i]), u=float(u[i]), gamma=float(np.gradient(lq, np.log(u))[i]), alpha=float(np.gradient(np.log(E), lq)[i]), index=float(np.gradient(np.log(lam * E), lq)[i]),
                                          h_T_over_E_note="E in m/s", E=float(E[i]), incision_over_head_per_12h=float(E[i] * 43200 / u[i]))
    out[cid] = row
json.dump(out, open(os.path.join(HERE, "rtm_index.json"), "w"), indent=1)
for cid, r in out.items(): print(cid, "thermal r1: gamma %.2f alpha %.2f index %.2f | diss alpha %.2f index %.2f | incision/head (12 h) %.2f" % (r["thermal_route1"]["gamma"], r["thermal_route1"]["alpha"], r["thermal_route1"]["index"], r["dissipation_route1"]["alpha"], r["dissipation_route1"]["index"], r["thermal_route1"]["incision_over_head_per_12h"]))
