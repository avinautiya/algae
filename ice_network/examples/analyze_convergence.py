"""Observed convergence orders from successive differences of the refinement study (no reference solution needed), and
Richardson-extrapolated errors at the default settings."""
import json, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
c = json.load(open(os.path.join(HERE, "..", "records", "convergence_study.json")))
out = {}
def orders(rows, key, h_of, r):
    vals = np.array([row[key] for row in rows]); hs = np.array([h_of(row) for row in rows])
    d = np.abs(np.diff(vals))
    ords = [float(np.log(d[i] / d[i + 1]) / np.log(hs[i] / hs[i + 1])) for i in range(len(d) - 1) if d[i] > 0 and d[i + 1] > 0]
    p = ords[-1] if ords else float("nan")
    # Richardson estimate of the converged value from the two finest members
    h1, h2 = hs[-2], hs[-1]; v1, v2 = vals[-2], vals[-1]
    ext = v2 + (v2 - v1) / ((h1 / h2) ** p - 1) if np.isfinite(p) and p > 0 else float("nan")
    return dict(observed_orders=ords, richardson_value=float(ext), values=[float(v) for v in vals])
for fam, h_of in (("dt_geom", lambda r: r["dt"]), ("stations", lambda r: 100.0 / (r["n"] - 1)), ("msub", lambda r: 1.0 / r["msub"])):
    rows = c[fam]; out[fam] = {}
    for key in ("sill_lowering_m", "delivered_m3", "final_storage_m3"):
        out[fam][key] = orders(rows, key, h_of, 2)
    qk = [dict(outlet_Q_6h=r["outlet_Q_hourly"][-1]) for r in rows]
    for key in ("outlet_Q_6h",):
        tmp = [dict(**{key: r["outlet_Q_hourly"][-1]}, dt=r["dt"], n=r["n"], msub=r["msub"]) for r in rows]
        out[fam][key] = orders(tmp, key, h_of, 2)
# error of the default configuration (n=9, dt=900, msub=16, rtol=1e-9) relative to the Richardson value, per family
for fam, idx in (("dt_geom", 2), ("stations", 1), ("msub", 2)):
    out[fam]["default_error_vs_richardson"] = {k: float(abs(out[fam][k]["values"][idx] - out[fam][k]["richardson_value"]) / abs(out[fam][k]["richardson_value"]))
                                               for k in ("sill_lowering_m", "delivered_m3", "final_storage_m3", "outlet_Q_6h")}
json.dump(out, open(os.path.join(HERE, "..", "records", "convergence_orders.json"), "w"), indent=1)
for fam, d in out.items():
    print(fam)
    for k, v in d.items():
        if k == "default_error_vs_richardson":
            print("   default error vs Richardson:", {a: float("%.2e" % b) for a, b in v.items()})
        else:
            print("  ", k, "orders", [round(o, 2) for o in v["observed_orders"]], "extrapolated", float("%.7g" % v["richardson_value"]))
