"""Figures drawn only from immutable campaign records (no solver is run here). Idealised controlled problems; not forecasts for any real glacier.
python figures_campaign.py hydro events incision cost"""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import analyze_campaign as A
OUT = os.path.join(HERE, "figures"); os.makedirs(OUT, exist_ok=True)
OI = dict(black="#000000", orange="#E69F00", sky="#56B4E9", green="#009E73", blue="#0072B2", verm="#D55E00", pink="#CC79A7")
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.25, "figure.dpi": 130})
NOTE = "Idealised simulations (controlled ice-water problems); not forecasts for a real glacier"
RED = {"r5_dt1800": OI["pink"], "r9_dt900": OI["sky"], "r17_dt450": OI["blue"], "r33_dt225": OI["verm"]}
I_OF = {s["id"]: s["inflow_m3_s"] for s in json.load(open(os.path.join(HERE, "competing_outlets_v1.json")))["scenarios"]}


def frozen():
    runs = A.load(lambda r: A.case_set(r["spec"]) == "frozen" and r["status"] == "completed")
    d = {}
    for r in runs:
        s = r["spec"]; key = (s["scn"], "ref%d" % s["ncell"] if s["kind"] == "reference" else s["label"].split("_", 2)[2]); d[key] = (r, A.series(r))
    return d


def fig_hydro(d):
    for scn in ("S1", "S2"):
        I = I_OF[scn]; fig, ax = plt.subplots(2, 2, figsize=(10.5, 6.4), sharex="col")
        for rr, col in (("ref100", OI["green"]), ("ref200", OI["orange"]), ("ref400", OI["black"])):
            r, s = d[(scn, rr)]; lw = 2.2 if rr == "ref400" else 1.3
            for j, k in enumerate(("Q1", "Q2")): ax[0, j].plot(s["t"] / 3600, s[k] / I, color=col, lw=lw, label={"ref100": "reference 100 cells", "ref200": "reference 200 cells", "ref400": "reference 400 cells (target)"}[rr])
            ax[1, 0].plot(s["t"] / 3600, (s["sill1"][0] - s["sill1"]) * 1e3, color=col, lw=lw); ax[1, 1].plot(s["t"] / 3600, (s["sill2"][0] - s["sill2"]) * 1e3, color=col, lw=lw)
        for k, col in RED.items():
            r, s = d[(scn, k)]
            for j, q in enumerate(("Q1", "Q2")): ax[0, j].plot(s["t"] / 3600, s[q] / I, color=col, lw=1.0, ls="--", label="reduced " + k)
            ax[1, 0].plot(s["t"] / 3600, (s["sill1"][0] - s["sill1"]) * 1e3, color=col, lw=1.0, ls="--"); ax[1, 1].plot(s["t"] / 3600, (s["sill2"][0] - s["sill2"]) * 1e3, color=col, lw=1.0, ls="--")
        ax[0, 0].set_title(f"{scn} route R1 outflow / inflow (-)"); ax[0, 1].set_title(f"{scn} route R2 outflow / inflow (-)")
        ax[1, 0].set_title("R1 sill lowering (mm)"); ax[1, 1].set_title("R2 sill lowering (mm)")
        for a in ax[1]: a.set_xlabel("simulated time (h)")
        ax[0, 0].legend(fontsize=7, loc="lower right"); fig.suptitle(f"{scn}: reduced network vs resolved reference, 8 h, inflow {I} m$^3$/s.  {NOTE}", fontsize=9)
        fig.tight_layout(); fig.savefig(os.path.join(OUT, f"hydro_{scn}.png")); plt.close(fig)


def fig_events(d):
    fig, ax = plt.subplots(1, 3, figsize=(11, 3.8))
    for si, scn in enumerate(("S1", "S2", "S3", "S4", "S5")):
        for j, (key, lab) in enumerate((("t_act2_s", "R2 activation t_act2"), ("t_dom2_s", "R2 overtakes R1 t_dom2"))):
            if j == 1 and scn not in ("S1", "S2"): continue
            vals = {rr: d[(scn, rr)][0]["metrics"][key] for rr in ("ref100", "ref200", "ref400")}
            a = ax[j]; x = si
            a.plot([x] * 3, [vals[r] / 3600 for r in ("ref100", "ref200", "ref400")], "-", color="0.6", lw=1)
            a.scatter([x] * 3, [vals[r] / 3600 for r in ("ref100", "ref200", "ref400")], c=[OI["green"], OI["orange"], OI["black"]], s=18, zorder=3)
            for k, col in RED.items():
                v = d[(scn, k)][0]["metrics"][key]
                if v is not None: a.scatter([x + 0.12 + 0.05 * list(RED).index(k)], [v / 3600], color=col, marker="s", s=16, zorder=3)
    for j, t in enumerate(("R2 activation time (h): first Q_R2 > 1 % inflow", "R2 overtakes R1 (h): first Q_R2 > Q_R1")): ax[j].set_title(t, fontsize=8); ax[j].set_xticks(range(5)); ax[j].set_xticklabels(["S1", "S2", "S3", "S4", "S5"])
    ax[0].set_ylim(0.8, 1.45); ax[1].set_xlim(-0.5, 1.7)
    ax[2].axis("off"); ax[2].text(0, 0.95, "Dots: reference 100 / 200 / 400 cells\n(green / orange / black).\nSquares: reduced r5, r9, r17, r33\n(pink, sky, blue, vermilion).\n\nS3-S5: R2 never overtakes R1 within 8 h\nin the 400-cell reference.\n\n" + NOTE, fontsize=8, va="top")
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "events.png")); plt.close(fig)


def fig_cost(d):
    S = json.load(open(os.path.join(HERE, "records", "campaign_summary_frozen.json")))
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
    for si, scn in enumerate(("S1", "S2", "S3", "S4", "S5")):
        v = S[f"{scn}|{{}}"]; mk = ["o", "s", "^", "D", "v"][si]
        for n, row in v["reference_cost_vs_error"].items(): ax[0].scatter(row["wall_s"], row["hydrograph_rmse_over_I"], marker=mk, color=OI["black"], s=22, label=None)
        for k, row in v["reduced"].items(): ax[0].scatter(row["solver_s"], row["hydrograph"]["hydrograph_rmse_over_I"], marker=mk, color=RED[k], s=22)
    ax[0].set_xscale("log"); ax[0].set_yscale("log"); ax[0].set_xlabel("solver wall time (s) - 2 concurrent workers, see serial timing"); ax[0].set_ylabel("hydrograph RMSE / inflow (-)")
    ax[0].set_title("Accuracy-cost, S1-S5 (marker = scenario)", fontsize=9)
    for k, col in RED.items(): ax[0].scatter([], [], color=col, label="reduced " + k)
    ax[0].scatter([], [], color=OI["black"], label="reference 100/200 cells"); ax[0].legend(fontsize=7)
    cfgs = list(RED); w = 0.16
    for si, scn in enumerate(("S1", "S2", "S3", "S4", "S5")):
        v = S[f"{scn}|{{}}"]
        for ci, k in enumerate(cfgs): ax[1].bar(si + (ci - 1.5) * w, abs(v["reduced"][k]["event_time_error_s"]["t_act2_s"]), w * 0.9, color=RED[k], label="reduced " + k if si == 0 else None)
    ax[1].set_yscale("log"); ax[1].set_xticks(range(5)); ax[1].set_xticklabels(["S1", "S2", "S3", "S4", "S5"]); ax[1].set_ylabel("|t_act2 error| vs 400-cell reference (s)"); ax[1].set_title("R2 activation-time error", fontsize=9); ax[1].legend(fontsize=7)
    fig.text(0.5, 0.01, NOTE, ha="center", fontsize=8); fig.tight_layout(rect=(0, 0.03, 1, 1)); fig.savefig(os.path.join(OUT, "accuracy_cost_frozen.png")); plt.close(fig)


if __name__ == "__main__":
    d = frozen()
    for w in sys.argv[1:]: {"hydro": fig_hydro, "events": fig_events, "cost": fig_cost}[w](d)
