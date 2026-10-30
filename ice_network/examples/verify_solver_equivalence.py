"""Does the faster discharge calculation change the converged physical solution? Run the same coupled melting case with the
fast solver and with the unchanged bisection reference, same tolerances and time step, and compare end states and series."""
import json, os, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src")); sys.path.insert(0, os.path.join(HERE, "..", "tests"))
from helpers import chan, const, lake
import ice_network.evolution.network as NW
from ice_network.hydraulics.channel import edge_flow as fast
from ice_network.hydraulics.channel_ref import edge_flow_ref as ref


def run(fn, label):
    NW.edge_flow = fn
    e = chan(0, 1, 100.5, 1.0, 0.5, L=100.0, n=9)
    net = NW.Network([lake(99.0, 2000.0)], [dict(h_fixed=95.0)], [e])
    t0 = time.time()
    r = net.run([2000.0 * 1.55], [1.0], 3 * 3600.0, 900.0, const(0.05), const(1.0), rtol=1e-10, atol=1e-11)
    b = r["budgets"]
    return dict(label=label, seconds=round(time.time() - t0, 1), sill=float(net.edges[0]["z"][0]), width0=float(net.edges[0]["w"][0]),
                delivered=float(b["delivered_m3"][0]), heat=float(b["heat_to_ice_J"]), stage=[float(h["stage"][0]) for h in r["hist"]],
                Q=[float(h["Q"][0]) for h in r["hist"]], V=[float(h["V"][0]) for h in r["hist"]], budgets={k: (float(v) if np.ndim(v) == 0 else [float(x) for x in v]) for k, v in b.items()})


a = run(fast, "fast"); b = run(ref, "reference_bisection")
rel = lambda x, y: float(np.max(np.abs(np.array(x) - np.array(y)) / np.maximum(np.abs(np.array(y)), 1e-300)))
out = dict(fast=a, reference=b, speedup=b["seconds"] / a["seconds"],
           rel_diff=dict(sill_lowering=abs(a["sill"] - b["sill"]) / (100.5 - b["sill"]), width0=abs(a["width0"] - b["width0"]) / abs(b["width0"] - 0.5),
                         delivered=abs(a["delivered"] - b["delivered"]) / b["delivered"], heat=abs(a["heat"] - b["heat"]) / b["heat"],
                         stage_max=rel(a["stage"], b["stage"]), Q_max=rel(a["Q"], b["Q"]), V_max=rel(a["V"], b["V"])))
json.dump(out, open(os.path.join(HERE, "..", "records", "solver_equivalence.json"), "w"), indent=1)
print(json.dumps({k: v for k, v in out.items() if k in ("speedup", "rel_diff")}, indent=1))
print("seconds fast/ref", a["seconds"], b["seconds"])
