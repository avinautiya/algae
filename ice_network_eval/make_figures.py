"""Figures from executed runs only. All cases are idealised controlled problems, not forecasts for a real glacier."""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, "..")
sys.path.insert(0, HERE); sys.path.insert(0, ROOT); sys.path.insert(0, os.path.join(ROOT, "ice_network", "src"))
import run_eval as E
OUT = os.path.join(HERE, "figures"); os.makedirs(OUT, exist_ok=True)
OI = dict(black="#000000", orange="#E69F00", sky="#56B4E9", green="#009E73", yellow="#F0E442", blue="#0072B2", verm="#D55E00", pink="#CC79A7")
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.25, "figure.dpi": 130})
NOTE = "Idealised controlled problem - not a forecast for any real glacier"


def fig_series():
    scn = E.CFG["scenarios"][0]; I = scn["inflow_m3_s"]
    ref = E.run_reference(scn, 200)
    reds = {rc["id"]: E.run_reduced(scn, rc["stations"], rc["dt_geom_s"]) for rc in E.CFG["reduced_model_configurations"] if rc["id"] in ("r9_dt900", "r17_dt450", "r33_dt225")}
    fig, ax = plt.subplots(2, 2, figsize=(10, 6.2), sharex=True)
    th = ref["t"] / 3600.0
    ax[0, 0].plot(th, ref["Q1"] / I, color=OI["black"], lw=2, label="resolved reference (200 cells)"); ax[0, 1].plot(th, ref["Q2"] / I, color=OI["black"], lw=2)
    cols = {"r9_dt900": OI["sky"], "r17_dt450": OI["blue"], "r33_dt225": OI["verm"]}
    for k, m in reds.items():
        ax[0, 0].plot(m["t"] / 3600.0, m["Q1"] / I, color=cols[k], lw=1.2, label="reduced " + k); ax[0, 1].plot(m["t"] / 3600.0, m["Q2"] / I, color=cols[k], lw=1.2)
        ax[1, 0].plot(m["t"] / 3600.0, (m["sill1"][0] - m["sill1"]) * 1e3, color=cols[k], lw=1.2); ax[1, 1].plot(m["t"] / 3600.0, (m["sill2"][0] - m["sill2"]) * 1e3, color=cols[k], lw=1.2)
    ax[1, 0].plot(th, (ref["sill1"][0] - ref["sill1"]) * 1e3, color=OI["black"], lw=2); ax[1, 1].plot(th, (ref["sill2"][0] - ref["sill2"]) * 1e3, color=OI["black"], lw=2)
    ax[0, 0].set_title("Route R1 (low sill, gentle): lake outflow / inflow"); ax[0, 1].set_title("Route R2 (sill +0.05 m, steep): lake outflow / inflow")
    ax[1, 0].set_title("R1 sill lowering (mm)"); ax[1, 1].set_title("R2 sill lowering (mm)")
    for a in ax[1]: a.set_xlabel("time (h)")
    ax[0, 0].legend(fontsize=7, loc="lower right"); fig.suptitle("Scenario S1 (inflow 0.05 m3/s at 1 C): reduced network vs resolved reference\n" + NOTE, fontsize=10)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig1_discharge_and_incision_S1.png")); plt.close(fig)


def fig_cost():
    S = json.load(open(os.path.join(HERE, "records", "eval_summary.json")))
    keys = ["sill_lowering1_m", "sill_lowering2_m", "delivered1_m3", "delivered2_m3"]
    fig, ax = plt.subplots(1, 2, figsize=(10, 4.2))
    cfgs = [c["id"] for c in E.CFG["reduced_model_configurations"]]
    med_err, med_wall, med_q = [], [], []
    for cid in cfgs:
        errs = []; walls = []; qs = []
        for sid, d in S.items():
            r = d["reduced"][cid]
            e = [abs(r["errors"][k]["rel"]) for k in keys if r["errors"][k].get("rel") is not None]
            errs.append(np.mean(e)); walls.append(r["wall_s"]); qs.append(np.mean([r["series"]["Q1"]["rmse_over_I"], r["series"]["Q2"]["rmse_over_I"]]))
            ax[0].scatter(r["wall_s"], np.mean(e), color=OI["blue"], alpha=0.25, s=14)
        med_err.append(np.median(errs)); med_wall.append(np.median(walls)); med_q.append(np.median(qs))
    ax[0].plot(med_wall, med_err, "o-", color=OI["blue"], lw=2, label="reduced network (median over scenarios)")
    for c, x, y in zip(cfgs, med_wall, med_err): ax[0].annotate(c, (x, y), fontsize=7, xytext=(3, 3), textcoords="offset points")
    # reference accuracy-cost: error of 100 and 200 cells against the Richardson value built from 200/400
    R = json.load(open(os.path.join(HERE, "records", "eval_results.json")))
    rw, re = {100: [], 200: []}, {100: [], 200: []}
    for sid, d in R.items():
        ext = S[sid]["reference_extrapolated"]
        for n in (100, 200):
            m = d["reference_metrics"][str(n)]
            e = [abs((m[k] - ext[k]) / ext[k]) for k in keys if ext[k]]
            re[n].append(np.mean(e)); rw[n].append(d["reference_wall_s"][str(n)])
    ax[0].plot([np.median(rw[n]) for n in (100, 200)], [np.median(re[n]) for n in (100, 200)], "s-", color=OI["black"], lw=2, label="resolved reference, 100 and 200 cells")
    ax[0].set_xscale("log"); ax[0].set_yscale("log"); ax[0].set_xlabel("wall-clock time per 8 h run, 1 thread (s)"); ax[0].set_ylabel("mean relative error of sill lowering and delivered volume")
    ax[0].set_title("Accuracy-cost (state metrics)"); ax[0].legend(fontsize=7)
    ax[1].plot(med_wall, med_q, "o-", color=OI["verm"], lw=2); ax[1].set_xscale("log"); ax[1].set_yscale("log")
    for c, x, y in zip(cfgs, med_wall, med_q): ax[1].annotate(c, (x, y), fontsize=7, xytext=(3, 3), textcoords="offset points")
    ax[1].set_xlabel("wall-clock time (s)"); ax[1].set_ylabel("RMSE of lake outflow / inflow vs 400-cell reference"); ax[1].set_title("Accuracy-cost (discharge series)")
    fig.suptitle("Reduced network accuracy-cost. The reference error is vs its own Richardson extrapolation (not independent)\n" + NOTE, fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig2_accuracy_cost.png")); plt.close(fig)


def fig_events():
    S = json.load(open(os.path.join(HERE, "records", "eval_summary.json")))
    fig, ax = plt.subplots(1, 2, figsize=(10, 4))
    cfgs = [c["id"] for c in E.CFG["reduced_model_configurations"]]
    for j, (ev, title) in enumerate((("t_act2_s", "Activation time of second route R2"), ("t_dom2_s", "Time R2 discharge first exceeds R1"))):
        for i, (sid, d) in enumerate(S.items()):
            ref = d["reference_extrapolated"][ev]; unc = d["reference_uncertainty"][ev]
            if ref is None: continue
            ax[j].errorbar(i, ref / 3600.0, yerr=(unc or 0) / 3600.0, fmt="s", color=OI["black"], capsize=3, label="reference (extrapolated, +-|400-200|)" if i == 0 else None)
            for k, cid in enumerate(cfgs):
                v = d["reduced"][cid]["errors"][ev].get("reduced")
                if v is not None: ax[j].scatter(i + 0.08 * (k + 1), v / 3600.0, color=[OI["orange"], OI["sky"], OI["blue"], OI["verm"]][k], s=18, label=cid if i == 0 else None)
        ax[j].set_xticks(range(len(S))); ax[j].set_xticklabels(list(S.keys())); ax[j].set_ylabel("time (h)"); ax[j].set_title(title)
    ax[0].legend(fontsize=7)
    fig.suptitle("Event times: reduced configurations vs resolved reference (scenarios without an event are omitted)\n" + NOTE, fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig3_event_times.png")); plt.close(fig)


def fig_incision():
    X = json.load(open(os.path.join(ROOT, "ice_network", "records", "incision_experiment.json")))
    fig, ax = plt.subplots(1, 3, figsize=(13, 4))
    for nu, c in ((0.25, OI["sky"]), (1.0, OI["blue"]), (4.0, OI["verm"])):
        rows = [r for r in X["E1"] if r["nu_factor"] == nu]
        ax[0].plot([r["exposure_a"] for r in rows], [r["ice_fraction_of_channel_dissipation"] for r in rows], "o", color=c, label="model, nu_factor %g" % nu)
    a = np.geomspace(0.02, 8, 100); ax[0].plot(a, 1 - (1 - np.exp(-a)) / a, "k--", label="uniform-source closed form")
    ax[0].set_xscale("log"); ax[0].set_xlabel("thermal exposure a = h_T P L / (rho c Q)"); ax[0].set_ylabel("fraction of channel dissipation reaching the ice"); ax[0].legend(fontsize=7)
    ax[0].set_title("Not all dissipation melts ice")
    tr = ["T0_no_melting", "T1_thermal", "T2_dissipation_all", "T3_dissipation_all_fixed_aspect"]; col = [OI["black"], OI["blue"], OI["orange"], OI["green"]]
    for k, (t, c) in enumerate(zip(tr, col)):
        ys = [[r for r in X["E2"] if r["treatment"] == t and r["T_inflow_C"] == T][0]["sill_lowering_m"] * 1e3 for T in (0.0, 0.3, 1.0)]
        ax[1].bar(np.arange(3) + 0.2 * k - 0.3, ys, 0.2, color=c, label=t.replace("_", " "))
    ax[1].set_xticks(range(3)); ax[1].set_xticklabels(["0 C", "0.3 C", "1.0 C"]); ax[1].set_yscale("symlog", linthresh=0.1); ax[1].set_ylim(0, 150); ax[1].set_ylabel("sill lowering after 8 h (mm)"); ax[1].set_xlabel("inflow temperature")
    ax[1].legend(fontsize=6); ax[1].set_title("Warm-water incision by treatment")
    cm = X["E3"]["closure_multiplier"]; ax[2].plot([r["nu_factor"] for r in cm], [r["sill_lowering_m"] * 1e3 for r in cm], "o-", color=OI["blue"], label="closure multiplier")
    ax[2].set_xscale("log"); ax[2].set_yscale("log"); ax[2].set_xlabel("heat-transfer multiplier nu_factor"); ax[2].set_ylabel("sill lowering (mm), T1, 1 C inflow")
    ax[2].set_xticks([0.25, 0.5, 1, 2, 4]); ax[2].set_xticklabels(["0.25", "0.5", "1", "2", "4"]); ax[2].minorticks_off()
    ax[2].set_yticks([20, 50, 100, 200]); ax[2].set_yticklabels(["20", "50", "100", "200"])
    ax[2].set_title("Closure sensitivity (numerical error ~1%)")
    fig.suptitle("Controlled incision experiment, 100 m reach, 8 h\n" + NOTE, fontsize=9); fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig4_incision_treatments.png")); plt.close(fig)


def fig_refine():
    C = json.load(open(os.path.join(ROOT, "ice_network", "records", "convergence_study.json")))
    fig, ax = plt.subplots(1, 3, figsize=(11, 3.6))
    for a, fam, h, lab in ((ax[0], "dt_geom", lambda r: r["dt"], "geometry step dt (s)"), (ax[1], "stations", lambda r: 100.0 / (r["n"] - 1), "station spacing (m)"), (ax[2], "msub", lambda r: 1.0 / r["msub"], "1 / substeps")):
        rows = C[fam]; hs = np.array([h(r) for r in rows])
        for key, c in (("sill_lowering_m", OI["blue"]), ("delivered_m3", OI["verm"])):
            v = np.array([r[key] for r in rows]); err = np.abs(v - v[-1]) / abs(v[-1]); a.loglog(hs[:-1], err[:-1], "o-", color=c, label=key)
        a.set_xlabel(lab); a.set_ylabel("relative difference from finest"); a.legend(fontsize=7)
    fig.suptitle("Numerical refinement of the coupled melting-channel case\n" + NOTE, fontsize=9); fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig5_refinement.png")); plt.close(fig)


if __name__ == "__main__":
    which = sys.argv[1:] or ["series", "cost", "events", "incision", "refine"]
    for w in which: {"series": fig_series, "cost": fig_cost, "events": fig_events, "incision": fig_incision, "refine": fig_refine}[w]()
    print("figures written to", OUT)
