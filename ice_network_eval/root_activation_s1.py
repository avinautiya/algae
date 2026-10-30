"""First-overtopping stress test on the frozen S1 geometry: reduced model (9 stations, dt_geom 900 s) for the first 3 h, raw law (h_floor = 0), floors 1e-4/1e-3/1e-2 and the unchanged bisection reference law.
Each run has an evaluation cap and a wall limit; a hit is recorded as a stall with time and state. One variant per process: python root_activation_s1.py VARIANT"""
import json, os, signal, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import run_eval as R
import ice_network.evolution.network as NW
import ice_network.hydraulics.channel_ref as CR
CAP = 60000; WALL = 240
variant = sys.argv[1]; hf = {"raw": 0.0, "floor1e-4": 1e-4, "floor1e-3": 1e-3, "floor1e-2": 1e-2, "raw_bisection": 0.0}[variant]
state = dict(calls=0, t=0.0)
class Stop(Exception): pass
orig = CR.edge_flow_ref if variant == "raw_bisection" else NW.edge_flow
def counted(*a):
    state["calls"] += 1
    if state["calls"] > CAP: raise Stop("evaluation cap")
    return orig(*a)
NW.edge_flow = counted
R.T_END = 3 * 3600.0
def alarm(*a): raise Stop("wall limit")
signal.signal(signal.SIGALRM, alarm); signal.alarm(WALL)
scn = R.CFG["scenarios"][0]; out = dict(variant=variant, h_floor=hf, eval_cap=CAP, wall_limit_s=WALL); t0 = time.perf_counter()
def prog(t, n): state["t"] = t
try:
    m = R.run_reduced(scn, 9, 900, h_floor=hf, progress=prog)
    out.update(status="completed", wall_s=m["wall"], calls=m["nfev"] if False else state["calls"], t_act1_s=R.crossing(m["t"], m["Q1"], 0.01 * scn["inflow_m3_s"]), t_Q1_gt_1em6=R.crossing(m["t"], m["Q1"], 1e-6),
               delivered1_m3=float(m["delivered"][0]), stage_end_m=float(m["stage"][-1]))
except Stop as e:
    out.update(status="stalled_or_capped", reason=str(e), wall_s=time.perf_counter() - t0, calls=state["calls"], last_macrostep_end_s=state["t"])
except Exception as e:
    out.update(status="failed", error=f"{type(e).__name__}: {e}", wall_s=time.perf_counter() - t0, calls=state["calls"], last_macrostep_end_s=state["t"])
os.makedirs("records/root_activation", exist_ok=True); json.dump(out, open(f"records/root_activation/{variant}.json", "w"), indent=1); print(json.dumps(out))
