"""Reproduce and diagnose the challenge stall red_C009_r9 (stalled at simulated t = 900 s, 240 s no progress). Bounded: evaluation cap and wall limit; prints the state at the cap."""
import json, os, signal, sys, time
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, os.path.join(HERE, "..", "ice_network_eval"))
import run_eval as R
import ice_network.evolution.network as NW
from ice_network.hydraulics.channel import edge_flow_diag
case = next(c for c in json.load(open(os.path.join(HERE, "..", "ice_network_eval", "cases", "challenge.json"))) if c["label"] == sys.argv[1])
R.set_config(case["cfg_patch"]); scn = case["scn_override"]; ov = case["overrides"]
CAP = int(sys.argv[2]) if len(sys.argv) > 2 else 40000
state = dict(calls=0, last=None, hist=[], r1=None, r1hist=[]); orig = NW.edge_flow
class Stop(Exception): pass
def counted(h_i, h_j, z, w, n_man, dx, ke, msub, h, hf, ch):
    state["calls"] += 1
    if state["calls"] > CAP: raise Stop
    o = orig(h_i, h_j, z, w, n_man, dx, ke, msub, h, hf, ch)
    if state["calls"] % 500 == 0: state["hist"].append((state["calls"], float(h_i), float(h_j), float(o[0]), int(w.size), float(z[0])))
    if abs(z[0] - 100.0) < 1e-9:
        state['r1'] = (float(h_i), float(h_j), z.copy(), w.copy(), n_man, dx, ke, msub); state['r1hist'].append((state['calls'], float(h_i) - 100.0, float(o[0])))
    state["last"] = (float(h_i), float(h_j), z.copy(), w.copy(), n_man, dx, ke, msub); return o
NW.edge_flow = counted
def alarm(*a): raise Stop
signal.signal(signal.SIGALRM, alarm); signal.alarm(150)
tp = dict(t=0.0); trace = []
orig_rhs = NW.Network.rhs
def rhs_t(self, t, y, *a):
    trace.append((float(t), float(y[0]))); state['net'] = self; state['y'] = y.copy(); state['t'] = float(t); state['a'] = a; return orig_rhs(self, t, y, *a)
NW.Network.rhs = rhs_t
try:
    R.run_reduced(scn, case["stations"], case["dt_geom"], progress=lambda t, n: tp.__setitem__("t", t), **ov); print("completed (not reproduced)")
except Stop:
    print("stopped at", state["calls"], "evaluations; last macro-step end t =", tp["t"], "s")
print("rhs t samples (every 3000th):", [(round(a,3), round(b,2)) for a, b in trace[::3000]][:12]); print("last 4 rhs (t, V):", trace[-4:]); print("route-1 history (call, head over sill, Q) last 6 of", len(state["r1hist"])); [print("  ", x) for x in state["r1hist"][-6:]]
hi, hj, z, w, n, dx, ke, ms = state['r1']; h = np.zeros(z.size); hf = np.zeros(z.size - 1); ch = np.zeros(z.size, np.int64); dg = np.zeros(4)
sill = z[0]; print("last edge call: H_a - sill = %.3e m, H_b = %.3f, sill %.4f" % (hi - sill, hj, sill))
for head in np.concatenate([np.geomspace(1e-6, 3e-3, 14), [5e-3, 1e-2, 3e-2]]):
    q = edge_flow_diag(sill + head, hj, z, w, n, dx, ke, ms, h, hf, ch, dg)[0]; print("   head %.5f m  Q %.6e  brent_it %d  resid %.1e  control %d  choked %d" % (head, q, dg[0], dg[1], dg[3], ch.sum()))

# ---- conditioning diagnosis at the stalled state: finite-difference sensitivity of every right-hand-side component to the lake volume
state["calls"] = -10**9; net = state["net"]; y = state["y"]; t = state["t"]; a = state["a"]
f0 = orig_rhs(net, t, y, *a)
names = {0: "dV/dt(lake)"}; i = 0
for d in (1e-9, 1e-7, 1e-5, 1e-3):
    y2 = y.copy(); y2[0] += d; f1 = orig_rhs(net, t, y2, *a); df = (f1 - f0) / d
    big = np.argsort(-np.abs(df))[:4]
    print("dV = %.0e m3: largest |dF/dV| components %s  values %s  (state scale: V %.3f)" % (d, list(big), ["%.2e" % df[k] for k in big], y[0]))
print("F0 components (first 8):", ["%.2e" % x for x in f0[:8]], " n_state", y.size)

e0 = net.edges[0]; Tl = 0.0 + float(y[1] / max(y[0], 1e-9)); print("lake T (relative) %.3f; edge_state scan near the stall head (raw Q from edge_flow, floored Q and heat rate in segment 0):" % Tl)
sill0 = e0["z"][0]
for dh in np.linspace(6.70e-5, 6.80e-5, 11):
    st = net.edge_state(e0, sill0 + dh, net.outlets[0]["h_fixed"], Tl)
    h = np.zeros(e0["n"]); hf = np.zeros(e0["n"] - 1); chk = np.zeros(e0["n"], np.int64); dgn = np.zeros(4)
    qraw = edge_flow_diag(sill0 + dh, net.outlets[0]["h_fixed"], e0["z"], e0["w"], e0["n_man"], e0["dx"], e0["ke"], 16, h, hf, chk, dgn)[0]
    print("  head %.9e  Q_raw %.4e  Q_floored %.4e  h0 %.3e  q_ice0 %.4e W  phi_in %.3e  qmelt %.3e" % (dh, qraw, abs(st["Q"]), st["h"][0], st["q_ice"][0], st["phi_in"], st["q_melt"]))
