"""Synchronised animation of the S2 route change (R2 overtakes R1) from stored series: resolved reference (200 cells) and reduced network (r17 and r33), same initial and boundary data.
python animate_route_change.py  ->  figures/route_change_S2.gif   (idealised simulation, not a forecast)"""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import analyze_campaign as A
runs = A.load(lambda r: A.case_set(r["spec"]) == "frozen" and r["spec"]["scn"] == "S2" and r["status"] == "completed")
get = lambda pred: next(r for r in runs if pred(r["spec"]))
ref = A.series(get(lambda s: s["kind"] == "reference" and s["ncell"] == 400)); r17 = A.series(get(lambda s: s["kind"] == "reduced" and s["stations"] == 17)); r33 = A.series(get(lambda s: s["kind"] == "reduced" and s["stations"] == 33))
I = 0.05; OI = dict(black="#000000", blue="#0072B2", verm="#D55E00", green="#009E73", orange="#E69F00")
fig, ax = plt.subplots(2, 2, figsize=(10, 6), sharex=True); step = 6; idx = np.arange(0, len(ref["t"]), step)
th = ref["t"] / 3600.0
def series(a, k, t=None): return np.interp(ref["t"], a["t"], a[k])
for j, (k, ttl) in enumerate((("Q1", "R1 outflow / inflow (-)"), ("Q2", "R2 outflow / inflow (-)"))):
    ax[0, j].set_title(ttl); ax[0, j].set_ylim(0, 0.75)
    ax[1, j].set_title("R%d sill lowering (mm)" % (j + 1)); ax[1, j].set_ylim(0, 68); ax[1, j].set_xlabel("simulated time (h)")
lines = {}
for j, k in enumerate(("Q1", "Q2")):
    for lab, a, col in (("resolved reference, 400 cells", ref, OI["black"]), ("reduced r17", r17, OI["blue"]), ("reduced r33", r33, OI["verm"])):
        lines[(j, 0, lab)] = ax[0, j].plot([], [], color=col, lw=2.2 if a is ref else 1.2, label=lab)[0]
        lines[(j, 1, lab)] = ax[1, j].plot([], [], color=col, lw=2.2 if a is ref else 1.2)[0]
ax[0, 0].legend(fontsize=8, loc="lower right")
for a in ax.ravel(): a.set_xlim(0, 8)
txt = fig.text(0.5, 0.955, "", ha="center", fontsize=10); fig.text(0.5, 0.01, "Idealised simulation of a controlled two-outlet ice-water problem (S2); not a forecast for any real glacier.", ha="center", fontsize=8)
cur = [ax[i, j].axvline(0, color="0.5", lw=0.8) for i in range(2) for j in range(2)]
def upd(n):
    i = idx[n]; t = th[i]
    for j, k in enumerate(("Q1", "Q2")):
        for lab, a in (("resolved reference, 400 cells", ref), ("reduced r17", r17), ("reduced r33", r33)):
            q = series(a, k)[:i + 1] / I; s = series(a, "sill%d" % (j + 1)); s = (s[0] - s[:i + 1]) * 1e3
            lines[(j, 0, lab)].set_data(th[:i + 1], q); lines[(j, 1, lab)].set_data(th[:i + 1], s)
    for c in cur: c.set_xdata([t, t])
    dom = lambda a: "R2" if series(a, "Q2")[i] > series(a, "Q1")[i] else "R1"
    txt.set_text(f"t = {int(t):d} h {int(round((t % 1) * 60)):02d} min of 8 h   |   larger outflow: reference {dom(ref)}, r17 {dom(r17)}, r33 {dom(r33)}")
    return []
fig.tight_layout(rect=(0, 0.03, 1, 0.94)); an = FuncAnimation(fig, upd, frames=len(idx), blit=False)
os.makedirs(os.path.join(HERE, "figures"), exist_ok=True); an.save(os.path.join(HERE, "figures", "route_change_S2.gif"), writer=PillowWriter(fps=10), dpi=70); print("saved", len(idx), "frames")
