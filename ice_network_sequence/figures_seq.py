"""Figures from stored records and the frozen predictions (idealised simulations; no real-glacier forecast).  python figures_seq.py [1 2 3 4 5]"""
import json, os, sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "ice_network_capture")); sys.path.insert(0, os.path.join(HERE, "..", "ice_network_eval"))
import forcing, sdesign as design, rtm_seq, rtm, theory_eta
OUT = os.path.join(HERE, "figures"); os.makedirs(OUT, exist_ok=True)
OI = dict(black="#000000", orange="#E69F00", sky="#56B4E9", green="#009E73", blue="#0072B2", verm="#D55E00", pink="#CC79A7")
plt.rcParams.update({"font.size": 9, "axes.grid": True, "grid.alpha": 0.25, "figure.dpi": 130})
NOTE = "Idealised two-outlet simulations (not a forecast for any real glacier)"
ED = json.load(open(os.path.join(HERE, "eval_dev.json"))); ET = json.load(open(os.path.join(HERE, "eval_test.json"))) if os.path.exists(os.path.join(HERE, "eval_test.json")) else {}
SC = json.load(open(os.path.join(HERE, "scores_test.json"))) if os.path.exists(os.path.join(HERE, "scores_test.json")) else {}
PR = json.load(open(os.path.join(HERE, "TEST_PREDICTIONS.json")))["predictions"]


def fig1_forcing():
    c = design.DEV["S01"][1]; seq, mean = design.physical(c); fig, ax = plt.subplots(1, 3, figsize=(13, 3.8)); cols = dict(AB=OI["blue"], BA=OI["verm"], AVG=OI["black"])
    p = dict(c["geom"], deficit=-0.10)
    for k, s in seq.items():
        Qf, Tf = forcing.make(s, 900.0); t = np.linspace(0, sum(x[0] for x in s), 2000); ax[0].plot(t / 3600, Qf(t), color=cols[k], label=k, lw=1.6); ax[1].plot(t / 3600, Tf(t), color=cols[k], lw=1.6)
        r = rtm_seq.run(p, s, lake_T0=mean["T_bar"]); ax[2].plot(r["t"] / 3600, r["TL"], color=cols[k], lw=1.6)
    ax[0].set_ylabel("inlet discharge Q (m$^3$/s)"); ax[1].set_ylabel("inlet temperature T (C)"); ax[2].set_ylabel("lake temperature T$_L$ (C, RTM)")
    for a in ax: a.set_xlabel("simulated time (h)")
    ax[0].legend(); ax[0].set_title("episodes A (warm, low Q) / B (cold, high Q), then common continuation", fontsize=8); ax[1].set_title("identical water and sensible-heat budgets (1e-13)", fontsize=8); ax[2].set_title("but different lake-temperature histories", fontsize=8)
    fig.text(0.5, 0.005, NOTE, ha="center", fontsize=8); fig.tight_layout(rect=(0, 0.03, 1, 1)); fig.savefig(os.path.join(OUT, "fig1_forcing_pairs.png")); plt.close(fig)


def fig2_mechanism():
    ids = list(ED); fig, ax = plt.subplots(1, 2, figsize=(12, 4)); w = 0.27; x = np.arange(len(ids))
    for j, (k, col, lab) in enumerate((("thermal_r17", OI["blue"], "thermal network"), ("fixed_r17", OI["green"], "fixed geometry (storage only)"), ("dissipation_r17", OI["orange"], "dissipation-only incision"))):
        ax[0].bar(x + (j - 1) * w, [ED[i][k]["AB"]["share_end"] - ED[i][k]["BA"]["share_end"] for i in ids], w, color=col, label=lab)
    ax[0].axhline(0.05, color="k", ls=":", lw=0.8); ax[0].axhline(-0.05, color="k", ls=":", lw=0.8); ax[0].set_xticks(x); ax[0].set_xticklabels(ids, rotation=45); ax[0].set_ylabel("sequence effect on end-of-run outlet-2 share (AB - BA)"); ax[0].legend(fontsize=7); ax[0].set_title("ablation: the partition effect needs thermal incision", fontsize=8)
    eta_net = [ED[i]["thermal_r17"]["AB"]["heat_to_ice_J"] / ED[i]["thermal_r17"]["BA"]["heat_to_ice_J"] - 1 for i in ids]; eta_pred = []
    for i in ids:
        c = design.DEV[i][1]; p = dict(c["geom"], deficit=-0.10); pp = rtm.full(p); V = pp["area"] * (pp["sill1"] - pp["deficit"] - pp["bottom"]); eta_pred.append(theory_eta.eta(p, c["pair"], c["tau_h"] * 3600, c["tail_h"] * 3600, V)["eta"])
    ax[1].scatter(eta_pred, eta_net, color=OI["blue"], label="development"); [ax[1].annotate(i, (a, b), fontsize=6) for i, a, b in zip(ids, eta_pred, eta_net)]
    if SC:
        for cid, r in SC.items():
            if r["net17"]: ax[1].scatter([r["predicted"]["eta"]], [r["net17"]["heat_rel"]], color=OI["verm"], marker="s", label="frozen test" if cid == "X01" else None)
    ax[1].plot([-0.05, 0.22], [-0.05, 0.22], "k-", lw=0.8); ax[1].set_xlabel("closed-form eta (thermal storage)"); ax[1].set_ylabel("network: relative heat difference to the channels (AB/BA - 1)"); ax[1].legend(fontsize=7); ax[1].set_title("lake thermal storage changes the heat delivered to the ice", fontsize=8)
    fig.text(0.5, 0.005, NOTE, ha="center", fontsize=8); fig.tight_layout(rect=(0, 0.03, 1, 1)); fig.savefig(os.path.join(OUT, "fig2_mechanism.png")); plt.close(fig)


def fig3_pred_vs_obs():
    fig, ax = plt.subplots(1, 3, figsize=(13, 4))
    for cid, r in SC.items():
        pe = r["predicted"]["share_end"]
        if r["net17"]: ax[0].scatter(r["net17"]["share_end"], pe, color=OI["blue"], s=28, zorder=3)
        ref = r["ref400"] or r["ref200"]
        if ref: ax[0].scatter(ref["share_end"], pe, color=OI["green"], marker="s", s=34, zorder=3)
        if r["net17"]: ax[1].scatter(r["net17"]["heat_rel"], r["predicted"]["heat_rel"], color=OI["blue"], s=28)
        if ref: ax[1].scatter(ref["heat_rel"], r["predicted"]["heat_rel"], color=OI["green"], marker="s", s=34)
        if r["net17"]: ax[2].scatter(r["net17"]["share_end"], r["predicted"]["product_rule"], color=OI["orange"], s=28)
    ax[0].plot([-0.06, 0.06], [-0.06, 0.06], "k-", lw=0.8); ax[0].axhspan(-0.05, 0.05, color="0.9"); ax[0].set_xlabel("observed share effect AB - BA"); ax[0].set_ylabel("RTM_seq prediction (frozen)"); ax[0].set_title("end-of-run outlet-2 share (grey: below the 0.05 threshold)", fontsize=8)
    ax[0].scatter([], [], color=OI["blue"], label="network r17"); ax[0].scatter([], [], color=OI["green"], marker="s", label="resolved reference"); ax[0].legend(fontsize=7)
    ax[1].plot([-0.05, 0.22], [-0.05, 0.22], "k-", lw=0.8); ax[1].set_xlabel("observed relative heat difference to channels"); ax[1].set_ylabel("RTM_seq prediction (frozen)"); ax[1].set_title("incision-energy effect", fontsize=8)
    ax[2].plot([-0.06, 0.06], [-0.06, 0.06], "k-", lw=0.8); ax[2].set_xlabel("observed share effect (network r17)"); ax[2].set_ylabel("closed-form product rule eta x dphi_inc"); ax[2].set_title("thermal-storage product rule (partial mechanism)", fontsize=8)
    fig.text(0.5, 0.005, NOTE + ". Predictions frozen before the test runs.", ha="center", fontsize=8); fig.tight_layout(rect=(0, 0.03, 1, 1)); fig.savefig(os.path.join(OUT, "fig3_predicted_vs_observed.png")); plt.close(fig)


def fig4_trajectories():
    ids = [i for i in ("X09", "X10", "X04") if i in ET and "thermal_r17" in ET[i]]; fig, ax = plt.subplots(1, len(ids), figsize=(4.4 * len(ids), 3.8), squeeze=False); cols = dict(AB=OI["blue"], BA=OI["verm"], AVG=OI["black"])
    for a, cid in zip(ax[0], ids):
        for sim, ls, lw in (("thermal_r17", "-", 1.0), ("thermal_ref200", "--", 1.8)):
            fl = ET[cid].get(sim)
            if not fl: continue
            for k in ("AB", "BA", "AVG"):
                if k in fl: a.plot(np.array(fl[k]["t_series"]) / 3600, fl[k]["share_series"], color=cols[k], ls=ls, lw=lw, label=f"{k} {'r17' if 'r17' in sim else 'ref200'}")
        a.axhline(0.5, color="0.5", lw=0.6); a.set_xlabel("simulated time (h)"); a.set_ylabel("outlet-2 share Q2/(Q1+Q2)"); a.set_title(f"{cid}: {ET[cid]['description'][:40]}", fontsize=7)
    ax[0][0].legend(fontsize=6)
    fig.text(0.5, 0.005, NOTE, ha="center", fontsize=8); fig.tight_layout(rect=(0, 0.03, 1, 1)); fig.savefig(os.path.join(OUT, "fig4_trajectories.png")); plt.close(fig)


def fig5_criterion():
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    for src, E, col, lab in ((ED, ED, OI["blue"], "development"), (ET, ET, OI["verm"], "test")):
        for cid, r in E.items():
            n = r.get("thermal_r17")
            if not n or "AB" not in n or "BA" not in n: continue
            c = r["config"]; p = dict(c["geom"], deficit=-0.10); pp = rtm.full(p); V = pp["area"] * (pp["sill1"] - pp["deficit"] - pp["bottom"]); Qm = r["mean"]["Q_bar"]; epsT = c["tau_h"] * 3600 * Qm / V
            ax[0].scatter(epsT, abs(n["AB"]["heat_to_ice_J"] / n["BA"]["heat_to_ice_J"] - 1), color=col, s=28, label=lab if cid in ("S01", "X01") else None); ax[1].scatter(epsT, abs(n["AB"]["share_end"] - n["BA"]["share_end"]), color=col, s=28)
    ax[0].axhline(0.10, color="k", ls=":", lw=0.8); ax[1].axhline(0.05, color="k", ls=":", lw=0.8); ax[0].set_xscale("log"); ax[1].set_xscale("log")
    for a, t in zip(ax, ("heat to channels: |AB/BA - 1| (threshold 0.10)", "routing: |share effect| (threshold 0.05)")): a.set_xlabel("epsilon_T = tau Q_bar / V  (episode / flushing time)"); a.set_title(t, fontsize=8)
    ax[0].legend(fontsize=7); fig.text(0.5, 0.005, NOTE, ha="center", fontsize=8); fig.tight_layout(rect=(0, 0.03, 1, 1)); fig.savefig(os.path.join(OUT, "fig5_effect_vs_epsilon.png")); plt.close(fig)


if __name__ == "__main__":
    for f in sys.argv[1:] or ("1", "2", "3", "4", "5"): {"1": fig1_forcing, "2": fig2_mechanism, "3": fig3_pred_vs_obs, "4": fig4_trajectories, "5": fig5_criterion}[f]()
