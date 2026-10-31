"""Numerical blocker: cumulative energy states carried an absolute tolerance of 1e-9 J, so noise in heat rates at micro-depth discharge (raw root jumps +-30 % between heads 1e-7 m apart) forced
1e-5 s steps (red_C009_r9 stall). Correction under test: atol_energy = 1e-3 J on energy-type states (volumes keep 1e-9 m3). Default None = legacy behaviour.
Part 1: the three challenge failures that are not scope assertions, with the correction (wall limit 240 s each). Part 2: equivalence on S1, S2 against legacy (r17) and on C-cases that completed before."""
import json, os, signal, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "ice_network_eval"))
import run_eval as R
CH = {c["label"]: c for c in json.load(open(os.path.join(HERE, "..", "ice_network_eval", "cases", "challenge.json")))}
class Stop(Exception): pass
def alarm(*a): raise Stop
signal.signal(signal.SIGALRM, alarm)
out = dict(correction="atol_energy=1e-3 J on cumulative energy states", failures=[], equivalence=[])
def run_case(label, atol_energy, limit=240):
    R.set_config(CH[label]["cfg_patch"]); scn = CH[label]["scn_override"]; ov = CH[label]["overrides"]
    signal.alarm(limit); t0 = time.perf_counter()
    try:
        m = R.run_reduced(scn, CH[label]["stations"], CH[label]["dt_geom"], atol_energy=atol_energy, **ov); signal.alarm(0)
        return dict(status="completed", wall_s=m["wall"], evals=m["nfev"], metrics=R.metrics(m, scn), budgets={k: m["budgets"][k] for k in ("water_residual_rel", "energy_residual_rel")})
    except Stop: return dict(status="wall limit %d s" % limit, wall_s=time.perf_counter() - t0)
    except Exception as e: signal.alarm(0); return dict(status=f"failed: {type(e).__name__}: {e}", wall_s=time.perf_counter() - t0)
for lab in ("red_C009_r9", "red_C019_r9", "red_C031_r17"):
    r = run_case(lab, 1e-3); r["label"] = lab; out["failures"].append(r); print(lab, r["status"], round(r["wall_s"], 1), r.get("evals"), flush=True)
    json.dump(out, open(os.path.join(HERE, "records_blocker.json"), "w"), indent=1, default=float)
# equivalence: frozen S1, S2 (r17), legacy vs corrected
import importlib; R2 = importlib.reload(R)
for scn_id in ("S1", "S2"):
    R.CFG["scenarios"]; scn = next(s for s in R.CFG["scenarios"] if s["id"] == scn_id)
    a = R.run_reduced(scn, 17, 450, atol_energy=None); b = R.run_reduced(scn, 17, 450, atol_energy=1e-3)
    ma, mb = R.metrics(a, scn), R.metrics(b, scn)
    d = {k: (None if ma[k] is None or mb[k] is None else mb[k] - ma[k]) for k in ("t_act2_s", "t_dom2_s", "sill_lowering1_m", "sill_lowering2_m", "delivered1_m3", "delivered2_m3")}
    out["equivalence"].append(dict(scenario=scn_id, legacy_evals=a["nfev"], corrected_evals=b["nfev"], legacy_wall=a["wall"], corrected_wall=b["wall"], differences=d, max_Q_diff_over_I=float(max(np.abs(np.interp(a["t"], b["t"], b["Q1"]) - a["Q1"]).max(), np.abs(np.interp(a["t"], b["t"], b["Q2"]) - a["Q2"]).max()) / scn["inflow_m3_s"])))
    print(scn_id, out["equivalence"][-1], flush=True)
json.dump(out, open(os.path.join(HERE, "records_blocker.json"), "w"), indent=1, default=float)
