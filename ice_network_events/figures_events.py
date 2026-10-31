"""Figures 1-4 of the event study (matplotlib; fixed method colours, direct labels, one axis per panel)."""
import json, os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "ice_network_eval"))
import analyze_campaign as A, design, series as S
COL = dict(TR="#1a1a1a", AVG="#d95f02", PW6="#7570b3", PW3="#1b9e77", EA="#e7298a"); plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.alpha": 0.2})
E = {w: json.load(open(os.path.join(HERE, f"eval_{w}.json"))) for w in ("dev", "test")}; SC = json.load(open(os.path.join(HERE, "scores_test.json"))); SD = json.load(open(os.path.join(HERE, "scores_dev.json")))
c, dev, test = design.selection(); runs = {r["spec"]["label"]: r for r in A.load(lambda r: A.case_set(r["spec"]) == "events" and r["status"] == "completed")}

# Fig 1: mean misplaces the crossing (E04, E07, D03)
fig, ax = plt.subplots(2, 3, figsize=(11, 5.2), sharex=True)
for j, (cid, tab) in enumerate((("D03", dev), ("E04", test), ("E07", test))):
    cfg = tab[cid]; ser = S.make_series(cfg["forcing"], 36.0); t = np.array(ser["t"]) / 3600; Q = np.array(ser["Q"]); m = S.means(ser)
    a = ax[0, j]; a.plot(t, Q * 1e3, color=COL["TR"]); a.axhline(cfg["crit"]["Istar"] * 1e3, color="#555", ls="--"); a.axhline(m["Q_bar"] * 1e3, color=COL["AVG"], ls=":")
    a.text(0.5, cfg["crit"]["Istar"] * 1e3, " I*", va="bottom", color="#555"); a.text(35.5, m["Q_bar"] * 1e3, "mean ", va="bottom", ha="right", color=COL["AVG"]); a.set_title(f"{cid} ({cfg['category']}, R={cfg['crit']['R']:.1f})"); 
    if j == 0: a.set_ylabel("inflow (L/s)")
    b = ax[1, j]
    for mth in ("TR", "AVG", "PW3", "EA"):
        r = runs.get(f"evt_{cid}_{mth}_r17")
        if r: s = A.series(r); b.plot(s["t"] / 3600, s["Q2"] * 1e3, color=COL[mth], lw=1.4 if mth == "TR" else 1.0, label=mth)
    b.set_xlabel("time (h)")
    if j == 0: b.set_ylabel("outlet-2 discharge (L/s)")
    if j == 2: b.legend(frameon=False, ncol=2, loc="upper right")
fig.suptitle("Fig. 1  Mean inflow vs I*: the averaged run misplaces outlet-2 activation (network, r17)"); fig.tight_layout(); fig.savefig(os.path.join(HERE, "figures", "fig1_mean_vs_threshold.png"), dpi=150); plt.close(fig)

# Fig 2: criterion R vs observed averaging error
fig, a = plt.subplots(figsize=(6.2, 4.2))
for nm, sc in (("dev", SD), ("test", SC)):
    for cid, v in sc["rows"].items():
        mk = "o" if nm == "test" else "s"; col = "#c0392b" if v["positive"] else "#2c7fb8"
        a.scatter(max(v["R"], 3e-3), max(v["share_err_pp"]["AVG"], 0.05), c=col, marker=mk, s=40, edgecolor="k", lw=0.5); a.annotate(cid, (max(v["R"], 3e-3), max(v["share_err_pp"]["AVG"], 0.05)), fontsize=7, xytext=(3, 3), textcoords="offset points")
a.set_xscale("log"); a.set_yscale("log"); a.axvline(1, color="#555", ls="--"); a.axhline(5, color="#555", ls=":"); a.set_xlabel("screen R (R < 3e-3 plotted at 3e-3)"); a.set_ylabel("|share error| of averaged run (pp)")
a.set_title("Fig. 2  Screen R vs averaging error\nred = averaging misrepresents (frozen label); squares dev, circles test", fontsize=9); fig.tight_layout(); fig.savefig(os.path.join(HERE, "figures", "fig2_criterion.png"), dpi=150); plt.close(fig)

# Fig 3: outlet-2 volume relative to network TR, by method (configs with TR volume > 1 % of total)
ids = [k for k, v in SC["rows"].items() if v["rel_vol_err"]]; ms = ("AVG", "PW6", "PW3", "EA")
fig, a = plt.subplots(figsize=(7.5, 3.8)); w = 0.2
for i, m in enumerate(ms): a.bar(np.arange(len(ids)) + (i - 1.5) * w, [SC["rows"][k]["rel_vol_err"][m] * 100 for k in ids], w, color=COL[m], label=m)
a.set_xticks(range(len(ids))); a.set_xticklabels(ids); a.set_ylabel("outlet-2 volume error vs network TR (%)"); a.legend(frameon=False, ncol=4); a.set_yscale("symlog", linthresh=5)
a.set_title("Fig. 3  Delivered outlet-2 volume error (test, TR volume > 1 % of total)", fontsize=9); fig.tight_layout(); fig.savefig(os.path.join(HERE, "figures", "fig3_volume_error.png"), dpi=150); plt.close(fig)

# Fig 4: robustness - share_deliv (pp) per method: r17 vs ref200 and nu variants
rob = SC["robustness"]; fig, axs = plt.subplots(1, 2, figsize=(10, 3.8))
cfgs = [k for k in rob if "ref" in rob[k]]; 
for j, key in enumerate(("ref", "nu")):
    a = axs[j]; xs = []
    if key == "ref": items = [(k, "ref") for k in cfgs]
    else: items = [(k, t) for k in rob for t in ("net_nu0.5", "net_nu2") if t in rob[k]]
    for i, (k, t) in enumerate(items):
        d = rob[k][t]
        for mth in ("TR", "AVG", "EA"):
            if mth in d["share_deliv_pp"] and mth in d["share_deliv_r17_pp"]: a.scatter(d["share_deliv_r17_pp"][mth], d["share_deliv_pp"][mth], c=COL[mth], s=22, edgecolor="k", lw=0.4)
    lim = a.get_xlim(); a.plot(lim, lim, color="#999", lw=0.8); a.set_xlabel("outlet-2 delivered share, network r17 (%)"); a.set_ylabel("resolved ref200 (%)" if key == "ref" else "closure x0.5 / x2 (%)")
    a.set_title("(a) numerical refinement" if key == "ref" else "(b) heat-transfer multiplier", fontsize=9)
    for mth in ("TR", "AVG", "EA"): a.scatter([], [], c=COL[mth], label=mth, edgecolor="k", lw=0.4)
    a.legend(frameon=False)
fig.suptitle("Fig. 4  Robustness of the outlet-2 share to refinement and heat-transfer uncertainty", fontsize=10); fig.tight_layout(); fig.savefig(os.path.join(HERE, "figures", "fig4_robustness.png"), dpi=150); plt.close(fig)
print("ok")
