"""Sensitivity of scenario S1 outcomes to the PHYSICAL closures (heat-transfer multiplier, Manning n), applied identically to both models, compared
with the reduced-vs-reference numerical error under the nominal closures; plus integrator-tolerance independence of the reduced model."""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, ".."))
import run_eval as E
scn = E.CFG["scenarios"][0]
variants = {"nominal": dict(nu=1.0, n=0.04), "nu_0.5": dict(nu=0.5, n=0.04), "nu_2": dict(nu=2.0, n=0.04), "manning_0.03": dict(nu=1.0, n=0.03), "manning_0.05": dict(nu=1.0, n=0.05)}
out = {}
for name, v in variants.items():
    E.G_COMMON["n_man"] = v["n"]
    ref = E.run_reference(scn, 200, nu_factor=v["nu"]); red = E.run_reduced(scn, 17, 450.0, nu_factor=v["nu"])
    out[name] = dict(settings=v, reference200=E.metrics(ref, scn), reduced_r17=E.metrics(red, scn), series_error=E.series_error(red, ref, scn["inflow_m3_s"]))
    print(name, {k: (round(x, 3) if isinstance(x, float) else x) for k, x in out[name]["reference200"].items() if k in ("t_act2_s", "t_dom2_s", "sill_lowering1_m", "sill_lowering2_m", "delivered1_m3", "delivered2_m3")}, flush=True)
E.G_COMMON["n_man"] = 0.04
# integrator tolerance independence of the reduced model (r9)
tol = {}
for rtol in (1e-5, 1e-7, 1e-9):
    m = E.run_reduced(scn, 9, 900.0, rtol=rtol); tol[str(rtol)] = dict(metrics=E.metrics(m, scn), wall_s=m["wall"], nfev=m["nfev"])
    print("rtol", rtol, round(m["wall"], 1), m["nfev"], flush=True)
out["rtol_independence_r9"] = tol
json.dump(out, open(os.path.join(HERE, "records", "sensitivity_S1.json"), "w"), indent=1, default=float)
