"""Figures from stored records and the frozen predictions (idealised simulations; no real-glacier forecast).  python figures_capture.py"""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "ice_network_eval"))
import rtm, design, analyze_campaign as A, evaluate as EV
OUT = os.path.join(HERE, "figures"); os.makedirs(OUT, exist_ok=True)
OI = dict(black="#000000", orange="#E69F00", sky="#56B4E9", green="#009E73", blue="#0072B2", verm="#D55E00", pink="#CC79A7")
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.25, "figure.dpi": 130})
NOTE = "Idealised two-outlet simulations (not a forecast for any real glacier)"
runs = {r["spec"]["label"]: r for r in A.load(lambda r: A.case_set(r["spec"]) == "capture" and r["status"] == "completed")}
SC = json.load(open(os.path.join(HERE, "scores_test.json"))) if os.path.exists(os.path.join(HERE, "scores_test.json")) else {}


def fig1_mechanism():
    fig, ax = plt.subplots(1, 3, figsize=(13, 4))
    # (a) schematic
    a = ax[0]; a.axis("off"); a.set_xlim(0, 10); a.set_ylim(0, 6)
    a.add_patch(plt.Rectangle((0.3, 2.4), 3.2, 2.2, fc="#cfe8ff", ec="k")); a.text(1.9, 3.5, "lake\nA dH/dt = I - Q1 - Q2\nwarm water T_L", ha="center", fontsize=8)
    a.annotate("", xy=(8.8, 4.4), xytext=(3.5, 4.0), arrowprops=dict(arrowstyle="->", lw=2, color=OI["blue"])); a.text(6.2, 4.7, "outlet 1: sill s1, melt rate E1(Q1)", ha="center", fontsize=8)
    a.annotate("", xy=(8.8, 1.4), xytext=(3.5, 2.8), arrowprops=dict(arrowstyle="->", lw=2, color=OI["verm"])); a.text(6.2, 1.0, "outlet 2: sill s2 = s1 + gap, E2(Q2)", ha="center", fontsize=8)
    a.text(0.3, 0.2, "exact: d(gap)/dt = E1 - E2\nd ln(Q2/Q1)/dt = (lam2-lam1) dH/dt + (r2-r1)", fontsize=8)
    a.set_title("(a) mechanism", fontsize=9)
    # (b) phase lines for two geometries, thermal vs dissipation
    for ax_i, (name, base, col) in zip((ax[1],), (("geometry A (route 2 steeper, slope 0.02)", {}, OI["blue"]),)):
        for g, (nm, b, c) in enumerate((("A: slope2 0.02", {}, OI["blue"]), ("B: slope2 0.01", dict(S2=0.01), OI["verm"]))):
            q = rtm.full(dict(design.H12, **b)); tabs = (rtm.tabulate(q, 1, 17, "thermal"), rtm.tabulate(q, 2, 17, "thermal")); pd = rtm.gap_dynamics(tabs, q["I"], gaps=np.linspace(-0.05, 0.30, 281))
            ax_i.plot(pd["gaps"] * 1e3, pd["G"] * 3600 * 1e3, color=c, lw=2, label=f"{nm} (thermal)")
            for fp in pd["fixed_points"]: ax_i.plot([fp["gap"] * 1e3], [0], "o", color=c, mfc="white", ms=8)
            if pd["kappa"] is not None: ax_i.axvline(pd["kappa"] * 1e3, color=c, ls=":", lw=1)
        ax_i.axhline(0, color="k", lw=0.8); ax_i.set_xlabel("sill gap  s2 - s1 (mm)"); ax_i.set_ylabel("d(gap)/dt = E1 - E2  (mm/h)"); ax_i.legend(fontsize=7)
        ax_i.set_title("(b) quasi-steady phase line; open circle = capture threshold (unstable), dotted = equal flows", fontsize=8)
    # (c) symmetric pair, thermal vs dissipation exponents
    idx = json.load(open(os.path.join(HERE, "rtm_index.json")))
    names = ["D01", "D06", "D12", "D15", "T14"]; xs = np.arange(len(names)); w = 0.35
    ax[2].bar(xs - w / 2, [idx[n]["thermal_route1"]["index"] for n in names], w, color=OI["blue"], label="advected-heat (thermal)"); ax[2].bar(xs + w / 2, [idx[n]["dissipation_route1"]["index"] for n in names], w, color=OI["orange"], label="dissipation-only")
    ax[2].axhline(0, color="k", lw=0.8); ax[2].set_xticks(xs); ax[2].set_xticklabels(names); ax[2].set_ylabel("d ln(lambda E)/d ln Q  (< 0 stabilising at fixed lake level)"); ax[2].legend(fontsize=7); ax[2].set_title("(c) replicator stability index, symmetric pair", fontsize=8)
    fig.text(0.5, 0.005, NOTE, ha="center", fontsize=8); fig.tight_layout(rect=(0, 0.03, 1, 1)); fig.savefig(os.path.join(OUT, "fig1_mechanism.png")); plt.close(fig)


def fig2_boundary():
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    for a, (title, base, offs, ids) in zip(ax, (("geometry A (small lake, 3 C)", {}, np.linspace(0.10, 0.26, 17), ("T01", "T02", "T03")), ("geometry B, route-2 slope 0.01 (small lake, 3 C)", dict(S2=0.01), np.linspace(0.05, 0.19, 15), ("T04", "T05", "T06")))):
        q = rtm.full(dict(design.H12, **base, T=3.0, area=300.0)); tabs = (rtm.tabulate(q, 1, 17, "thermal"), rtm.tabulate(q, 2, 17, "thermal")); pd = rtm.gap_dynamics(tabs, q["I"]); dr = []
        for o in offs:
            t2 = dict(tabs[1]); t2["sill"] = q["sill1"] + o; r = rtm.rtm_run(dict(q, offset=o), "thermal", 17, tabs=(tabs[0], t2)); dr.append(rtm.events(r, q["I"])["gap_drift_mm"])
        a.plot(offs * 1e3, dr, color=OI["blue"], lw=2, label="RTM (frozen, lake lag included)")
        if pd["fixed_points"]: a.axvline(pd["fixed_points"][0]["gap"] * 1e3, color=OI["verm"], ls="--", label="quasi-steady-lake capture threshold")
        for cid in ids:
            for lab, mk, col in ((f"cap_{cid}_thermal_r17", "o", OI["black"]), (f"cap_{cid}_thermal_ref400", "s", OI["green"]), (f"cap_{cid}_thermal_ref200", "^", OI["sky"])):
                r = runs.get(lab)
                if r: s = EV.sim_summary(r); a.scatter([design.TEST[cid][1]["offset"] * 1e3], [s["gap_drift_mm"]], marker=mk, color=col, s=40, zorder=4, label={"o": "network r17", "s": "resolved ref 400", "^": "resolved ref 200"}[mk] if cid == ids[0] else None)
        a.axhline(0, color="k", lw=0.8); a.set_xlabel("initial sill gap s2 - s1 (mm)"); a.set_ylabel("gap drift over second half of 12 h (mm)"); a.set_title(title, fontsize=9); a.legend(fontsize=7)
    fig.text(0.5, 0.005, NOTE + ". Test points were run after the predictions were frozen.", ha="center", fontsize=8); fig.tight_layout(rect=(0, 0.03, 1, 1)); fig.savefig(os.path.join(OUT, "fig2_transition_boundary.png")); plt.close(fig)


def fig3_outcomes():
    fig, ax = plt.subplots(1, 3, figsize=(12, 3.8))
    P = json.load(open(os.path.join(HERE, "TEST_PREDICTIONS.json")))
    for cid, p in P.items():
        r = runs.get(f"cap_{cid}_thermal_r17")
        if not r: continue
        s = EV.sim_summary(r); th = p["rtm_thermal"]
        ax[0].scatter(s["phi2_deliv"], th["phi2_end"], color=OI["blue"], s=25); ax[1].scatter(s["d1"] * 1e3, th["sill1_mm"], color=OI["blue"], s=25); ax[1].scatter(s["d2"] * 1e3, th["sill2_mm"], color=OI["verm"], s=25)
        ax[2].scatter(s["gap_drift_mm"], th["gap_drift_mm"], color=OI["blue"], s=25)
        rr = runs.get(f"cap_{cid}_thermal_ref400") or runs.get(f"cap_{cid}_thermal_ref200")
        if rr: s2 = EV.sim_summary(rr); t33 = p["rtm_thermal33"]; ax[0].scatter(s2["phi2_deliv"], t33["phi2_end"], color=OI["green"], marker="s", s=30); ax[1].scatter(s2["d1"] * 1e3, t33["sill1_mm"], color=OI["green"], marker="s", s=30); ax[1].scatter(s2["d2"] * 1e3, t33["sill2_mm"], color=OI["green"], marker="s", s=30); ax[2].scatter(s2["gap_drift_mm"], t33["gap_drift_mm"], color=OI["green"], marker="s", s=30)
    for a, t, lim in zip(ax, ("outlet-2 delivered share", "inlet-sill lowering (mm), route 1 blue / route 2 red", "gap drift (mm)"), ((0, 1), (0, 450), (-60, 130))):
        a.plot(lim, lim, "k-", lw=0.8); a.set_xlabel("simulated"); a.set_ylabel("RTM prediction (frozen)"); a.set_title(t, fontsize=8)
    ax[0].scatter([], [], color=OI["blue"], label="RTM17 vs network r17"); ax[0].scatter([], [], color=OI["green"], marker="s", label="RTM33 vs resolved reference"); ax[0].legend(fontsize=7)
    fig.text(0.5, 0.005, NOTE, ha="center", fontsize=8); fig.tight_layout(rect=(0, 0.03, 1, 1)); fig.savefig(os.path.join(OUT, "fig3_predicted_vs_simulated.png")); plt.close(fig)


def fig4_trajectories():
    fig, ax = plt.subplots(2, 3, figsize=(13, 6.2)); P = json.load(open(os.path.join(HERE, "TEST_PREDICTIONS.json")))
    for j, cid in enumerate(("T07", "T05", "T13")):
        q = rtm.full(dict(design.H12, **design.TEST[cid][1])); r = rtm.rtm_run(q, "thermal", 17); I = q["I"]
        ax[0, j].plot(r["t"] / 3600, r["Q1"] / I, color=OI["blue"], lw=1.2, ls="--", label="RTM17 Q1"); ax[0, j].plot(r["t"] / 3600, r["Q2"] / I, color=OI["verm"], lw=1.2, ls="--", label="RTM17 Q2")
        ax[1, j].plot(r["t"] / 3600, (q["offset"] - (r["s2"] - r["s1"])) * 1e3, color=OI["black"], lw=1.2, ls="--", label="RTM17")
        for lab, ls, lw in ((f"cap_{cid}_thermal_r17", "-", 1.0), (f"cap_{cid}_thermal_ref400", "-", 2.2)):
            rr = runs.get(lab)
            if rr:
                s = A.series(rr); col = (OI["sky"], OI["orange"]) if "r17" in lab else (OI["blue"], OI["verm"]); nm = "net r17" if "r17" in lab else "resolved ref 400"
                ax[0, j].plot(s["t"] / 3600, s["Q1"] / I, color=col[0], lw=lw, label=nm + " Q1"); ax[0, j].plot(s["t"] / 3600, s["Q2"] / I, color=col[1], lw=lw, label=nm + " Q2")
                ax[1, j].plot(s["t"] / 3600, (s["sill2"] - s["sill1"]) * 1e3, color=col[0] if "r17" in lab else OI["black"], lw=lw, label=nm)
        ax[0, j].set_title(f"{cid}: {design.TEST[cid][0][:60]}", fontsize=7); ax[1, j].set_xlabel("simulated time (h)"); ax[0, j].set_ylabel("outflow / inflow (-)"); ax[1, j].set_ylabel("sill gap s2 - s1 (mm)")
    ax[0, 0].legend(fontsize=6); ax[1, 0].legend(fontsize=6)
    fig.text(0.5, 0.005, NOTE, ha="center", fontsize=8); fig.tight_layout(rect=(0, 0.03, 1, 1)); fig.savefig(os.path.join(OUT, "fig4_trajectories.png")); plt.close(fig)


if __name__ == "__main__":
    for f in sys.argv[1:] or ("1", "2", "3", "4"): {"1": fig1_mechanism, "2": fig2_boundary, "3": fig3_outcomes, "4": fig4_trajectories}[f]()
