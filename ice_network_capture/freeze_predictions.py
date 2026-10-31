"""Collect the frozen test predictions (RTM + quasi-steady-lake capture quantities) into TEST_PREDICTIONS.json/.md. Run BEFORE any test case is executed."""
import json, os, sys, hashlib
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import rtm, design
HORIZON = 12 * 3600.0; rows = []; out = {}
for cid, (desc, p) in design.TEST.items():
    pr = json.load(open(os.path.join(HERE, "pred_test", cid + ".json"))); q = pr["params"]
    # quasi-steady-lake capture quantities from the frozen geometry (thermal)
    tabs = (rtm.tabulate(q, 1, 17, "thermal"), rtm.tabulate(q, 2, 17, "thermal")); pd = rtm.gap_dynamics(tabs, q["I"])
    th = pr["thermal_17"]["events"]; fx = pr["fixed_17"]["events"]; ds = pr["dissipation_17"]["events"]
    tdom = th["t_dom2_s"]; near_sw = (tdom is not None and tdom > 0.85 * HORIZON) or (tdom is None and False)
    drift = th["gap_drift_mm"]; sign = "closing" if drift < -2 else ("opening" if drift > 2 else "flat")
    rec = dict(id=cid, description=desc, qs_kappa=pd["kappa"], qs_capture_threshold_gap=[f["gap"] for f in pd["fixed_points"] if not f["stable"]], qs_sigma_per_s=[f["slope"] for f in pd["fixed_points"]],
               gap0=q["offset"], rtm_thermal=dict(switch=tdom is not None, t_dom2_s=tdom, near_tie_switch=bool(near_sw), phi2_end=th["phi2_end"], sill1_mm=th["sill1_m"] * 1e3, sill2_mm=th["sill2_m"] * 1e3, gap_drift_mm=drift, gap_drift_sign=sign, t_act2_s=th["t_act2_s"]),
               rtm_fixed=dict(switch=fx["t_dom2_s"] is not None, t_dom2_s=fx["t_dom2_s"], phi2_end=fx["phi2_end"]), rtm_dissipation=dict(switch=ds["t_dom2_s"] is not None, t_dom2_s=ds["t_dom2_s"], phi2_end=ds["phi2_end"], sill1_mm=ds["sill1_m"] * 1e3, sill2_mm=ds["sill2_m"] * 1e3),
               rtm_thermal33=dict(switch=pr["thermal_33"]["events"]["t_dom2_s"] is not None, t_dom2_s=pr["thermal_33"]["events"]["t_dom2_s"], phi2_end=pr["thermal_33"]["events"]["phi2_end"], sill1_mm=pr["thermal_33"]["events"]["sill1_m"] * 1e3, sill2_mm=pr["thermal_33"]["events"]["sill2_m"] * 1e3, gap_drift_mm=pr["thermal_33"]["events"]["gap_drift_mm"]))
    out[cid] = rec
json.dump(out, open(os.path.join(HERE, "TEST_PREDICTIONS.json"), "w"), indent=1, default=float)
L = ["| id | configuration | QS kappa / threshold gap (m) | gap0 | RTM17: switch (t_dom2 s) | fixed-geometry switch | phi2 end | sills mm (1/2) | gap drift (mm) |", "|---|---|---|---|---|---|---|---|---|"]
for cid, r in out.items():
    t = r["rtm_thermal"]; L.append(f"| {cid} | {r['description']} | {None if r['qs_kappa'] is None else round(r['qs_kappa'], 4)} / {[round(x, 4) for x in r['qs_capture_threshold_gap']]} | {r['gap0']} | {'yes' if t['switch'] else 'no'} ({None if t['t_dom2_s'] is None else round(t['t_dom2_s'])}) | {'yes' if r['rtm_fixed']['switch'] else 'no'} | {t['phi2_end']:.2f} | {t['sill1_mm']:.0f} / {t['sill2_mm']:.0f} | {t['gap_drift_mm']:+.1f} ({t['gap_drift_sign']}) |")
open(os.path.join(HERE, "TEST_PREDICTIONS.md"), "w").write("\n".join(L)); print("\n".join(L))
