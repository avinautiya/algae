"""Figures 1-6 (7 only if supported) as vector PDF + PNG. `python -m ice_cliff_survival.reactivation.make_figures_m3 EVAL_DIR SPEC_DIR OUT_DIR`. Panel labels: IDEALISED (geometry/verification), DEVELOPMENT (conv/dev runs), EVALUATION (frozen E01-E10)."""
import glob, json, os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from . import metrics as M, cases
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, "results")
C = dict(bare="#555555", corridor="#1f5f8b", strip="#c8641e", porous="#3a8f5b", ref="#000000", red="#a33", cfd="#111111", cn="#1f5f8b", fr="#c8641e", fx="#3a8f5b")
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.alpha": 0.25, "figure.dpi": 130, "pdf.fonttype": 42})
def save(fig, out, name): fig.tight_layout(); fig.savefig(os.path.join(out, name + ".pdf")); fig.savefig(os.path.join(out, name + ".png"), dpi=170); plt.close(fig)
def runs(ev, suffix="medium"):
    o = {}
    for f in sorted(glob.glob(os.path.join(ev, "E*_%s.json" % suffix))):
        r = json.load(open(f))
        if r.get("status") == "complete": r["name"] = os.path.basename(f)[:-5].rsplit("_", 1)[0]; o[r["name"]] = r
    return o
def fig1(out):
    fig, ax = plt.subplots(3, 1, figsize=(9.5, 6.5), sharex=True)
    for a, lay in zip(ax, ("bare", "corridor", "strip")):
        s = cases.build(lay, dx=0.05); X = (np.arange(s.nx + 1)) * s.dx; Y = (np.arange(s.ny + 1)) * s.dx; Z = s.zs()[1:-1, 1:-1] - s.zb[1:-1, 1:-1] * 0
        m = a.pcolormesh(X, Y, (s.zb + np.where(s.solid, s.zd, 0))[1:-1, 1:-1].T, cmap="Blues_r", vmin=-1.0, vmax=0.3, shading="flat"); a.add_patch(Rectangle((s.i0 * s.dx - s.dx, 0), s.Lr, s.ny * s.dx, fill=False, ec="k", lw=0.8))
        a.set(ylabel="y [m]", title=f"IDEALISED geometry: {lay} (bed+deposit elevation z, m; datum: bypass floor = 0)"); a.set_aspect("equal"); fig.colorbar(m, ax=a, label="z [m]", fraction=0.025)
    ax[-1].set(xlabel="x [m]"); ax[0].text(0.05, 0.55, "inlet Q=0.07 m³/s, T=0.5 K →", fontsize=8, transform=ax[0].transAxes); ax[0].text(0.78, 0.55, "overfall (zero-gradient) →", fontsize=8, transform=ax[0].transAxes)
    ax[0].text(1.2 + 0.1, 0.12, "toe lane", fontsize=8); ax[0].text(1.2 + 0.1, 0.62, "bypass lane (divider wall between)", fontsize=8); save(fig, out, "fig1_geometry_bc")
def fig2(out):
    v = json.load(open(os.path.join(RES, "verify_swe_v5.json"))); vb = json.load(open(os.path.join(RES, "verify_swe_v5b_symmetric_plug.json"))); cv = json.load(open(os.path.join(RES, "convergence_dev.json")))
    fig, ax = plt.subplots(2, 3, figsize=(11, 6.2)); a = ax[0, 0]
    for lay, col in (("bare", C["bare"]), ("corridor", C["corridor"]), ("strip", C["strip"])):
        P = [p for p in v["partition"] if p["layout"] == lay]; a.plot([p["dx"] for p in P], [abs(p["err"]) for p in P], "o-", color=col, label=lay)
    a.axhline(0.05, color="k", ls="--", lw=0.8); a.set(xscale="log", xlabel="dx [m]", ylabel="|f_toe − reference| (abs.)", title="partition vs independent column-weir reference"); a.legend(frameon=False); a.invert_xaxis()
    a = ax[0, 1]; sy = [r for r in vb["symmetric"] if r["n"] == 0.001]; a.plot([r["dx"] for r in sy], [r["Cd_toe"] for r in sy], "o-", color=C["cn"]); a.axhline(1.0, color="k", lw=0.8); a.set(xscale="log", xlabel="dx [m]", ylabel="discharge coeff. vs ideal weir", title="symmetric lanes (frictionless): → 1"); a.invert_xaxis()
    a = ax[0, 2]; a.plot([p["dx"] for p in vb["plug"]], [abs(p["lanes"][0]["rel_err"]) for p in vb["plug"]], "o-", color=C["cn"], label="lane heat"); a.plot([p["dx"] for p in vb["plug"]], [p["lanes"][0]["lowering_profile_max_rel_err"] for p in vb["plug"]], "s-", color=C["fr"], label="lowering profile"); a.set(xscale="log", yscale="log", xlabel="dx [m]", ylabel="rel. error vs analytic exp(−NTU)", title="plug-flow heat-depletion limit"); a.legend(frameon=False); a.invert_xaxis()
    for a, key, lab in ((ax[1, 0], "dx", "dx [m]"), (ax[1, 1], "dtm", "bed update interval [h]"), (ax[1, 2], "cfl", "CFL")):
        for tag, col in (("strip_n", C["strip"]), ("strip_r8", C["corridor"])):
            rows = cv[tag][key]; a.plot([x["value"] for x in rows], [x["loss"] for x in rows], "o-", color=col, label=("no removal" if tag.endswith("n") else "removal 8 d"))
        a.set(xlabel=lab, ylabel="contact-loss day", title="DEVELOPMENT: event-time convergence (%s)" % key); a.set_xscale("log" if key != "dtm" else "linear")
    ax[1, 0].legend(frameon=False); save(fig, out, "fig2_verification_convergence")
def fig3(ev, out):
    R = runs(ev); fig, ax = plt.subplots(2, 2, figsize=(11, 6.5)); names = [("E04_strip_none", "none", "k"), ("E05_strip_r1", "early removal 1 d", C["corridor"]), ("E06_strip_r8", "removal 8 d", C["strip"]), ("E07_strip_r12", "late removal 12 d", C["porous"])]
    a = ax[0, 0]
    for n, lab, col in names:
        if n in R: x = M.arrays(R[n]); a.plot(x["day"], x["f"], color=col, label=lab)
    a.set(xlabel="day", ylabel="toe-route flow fraction f_toe", title="EVALUATION: solid strip, matched removal"); a.legend(frameon=False); a = ax[0, 1]
    for n, lab, col in names:
        if n in R: x = M.arrays(R[n]); a.plot(x["day"], x["crest_toe"], color=col, label=lab + " toe crest"); 
    if "E04_strip_none" in R: a.plot(M.arrays(R["E04_strip_none"])["day"], M.arrays(R["E04_strip_none"])["crest_byp"], "k--", label="bypass crest")
    a.set(xlabel="day", ylabel="controlling crest elevation [m]", title="controlling ice elevations"); a.legend(frameon=False, fontsize=7)
    for a, n in ((ax[1, 0], "E06_strip_r8"), (ax[1, 1], "E07_strip_r12")):
        if n not in R: continue
        recs = R[n]["records"]; xs = np.linspace(1.0, 3.0, len(recs[0]["toe_floor_x"]))
        for d, col in ((0, "#999999"), (4, "#6baed6"), (8, "#2171b5"), (12, "#08306b"), (21, "#e6550d")):
            i = int(round(d * 4)); i = min(i, len(recs) - 1); a.plot(xs, recs[i]["toe_floor_x"], color=col, label=f"toe floor day {d}"); a.plot(xs, recs[i]["toe_crest_x"], color=col, ls=":", lw=0.8)
        a.plot(xs, recs[-1]["byp_floor_x"], "k--", label="bypass floor day 21"); a.set(xlabel="x along reach [m]", ylabel="elevation [m]", title=f"EVALUATION interface evolution: {n} (dotted: crest incl. deposit)"); a.legend(frameon=False, fontsize=6)
    save(fig, out, "fig3_flow_interface_removal")
def fig4(ev, out):
    R = runs(ev); fig, ax = plt.subplots(1, 3, figsize=(11, 3.4)); col = [C["bare"], C["corridor"], C["strip"], C["porous"]]
    for (n, c) in zip(("E01_bare_none", "E02_corridor_none", "E04_strip_none", "E08_porous_none"), col):
        if n not in R: continue
        r = R[n]["records"]; d = np.array([x["day"] for x in r]); ax[0].semilogy(d, np.maximum(np.abs([x["water_resid"] for x in r]), 1e-18), color=c, label=n[4:]); ax[1].semilogy(d, np.maximum(np.abs([x["energy_resid_rel"] for x in r]), 1e-18), color=c)
        ax[2].plot(d, np.array([x["heat_in"] for x in r]) / 1e3, color=c, label=n[4:] + " in"); ax[2].plot(d, (np.array([x["heat_out"] for x in r]) + np.array([x["melt_heat_W"] for x in r]) + np.array([x["storage_W"] for x in r])) / 1e3, color=c, ls=":")
    ax[0].set(xlabel="day", ylabel="|water residual| / Q", title="water mass (incl. storage)"); ax[1].set(xlabel="day", ylabel="|energy residual| / heat in", title="sensible + latent + storage"); ax[2].set(xlabel="day", ylabel="heat rate [kW]", title="inflow (solid) vs out+latent+storage (dotted)"); ax[0].legend(frameon=False, fontsize=7); save(fig, out, "fig4_budgets")
def fig5(ev, out):
    R = runs(ev); fig, ax = plt.subplots(1, 3, figsize=(11, 3.5))
    for n, lab, c in (("E01_bare_none", "bare", C["bare"]), ("E02_corridor_none", "connected corridor", C["corridor"]), ("E04_strip_none", "solid strip", C["strip"]), ("E08_porous_none", "porous strip (drag)", C["porous"])):
        if n in R: x = M.arrays(R[n]); ax[0].plot(x["day"], x["f"], color=c, label=lab); ax[1].plot(x["day"], x["melt_toe"], color=c); ax[2].plot(x["day"], x["wet_toe"], color=c)
    ax[0].set(xlabel="day", ylabel="f_toe", title="EVALUATION: equal volume & footprint"); ax[1].set(xlabel="day", ylabel="cumulative local toe-reach melt [m³]", title="local ice loss (not glacier-wide)"); ax[2].set(xlabel="day", ylabel="wetted toe ice area [m²]", title="wetted area"); ax[0].legend(frameon=False, fontsize=7); save(fig, out, "fig5_equal_cover_connectivity")
def fig6(analysis, out):
    A = json.load(open(analysis)); comp = A["comparisons"]; fig, ax = plt.subplots(1, 2, figsize=(11, 4)); runs_ = sorted({c["run"] for c in comp}); meths = ["connectivity_reduced", "fractional_cover", "fixed_geometry", "bare_reduced(evolving)", "bare_fixed_geometry"]; cols = {"connectivity_reduced": C["cn"], "fractional_cover": C["fr"], "fixed_geometry": C["fx"], "bare_reduced(evolving)": C["cn"], "bare_fixed_geometry": C["fx"]}
    for i, n in enumerate(runs_):
        for j, c in enumerate([c for c in comp if c["run"] == n]):
            ax[0].bar(i + (j - 1) * 0.27, c["f_err"]["mean_abs"], 0.25, color=cols.get(c["method"], "#999"), label=c["method"] if i == 0 or c["method"] not in [k for k in cols if False] else None)
    ax[0].axhline(A["U_f"], color="k", ls="--", lw=0.8, label="CFD numerical uncertainty U_f"); ax[0].set_xticks(range(len(runs_))); ax[0].set_xticklabels([r[:3] for r in runs_]); ax[0].set(ylabel="mean |f_pred − f_CFD| over 21 d", title="EVALUATION: continuous partition error"); h, l = ax[0].get_legend_handles_labels(); u = dict(zip(l, h)); ax[0].legend(u.values(), u.keys(), frameon=False, fontsize=6)
    pts = [(c["run"], c["method"], c["loss_err"]) for c in comp if c["loss_err"] is not None]
    for k, (n, m, e) in enumerate(pts): ax[1].errorbar(k, e, yerr=A["U_loss_day"], fmt="o", color=cols.get(m, "#999"), capsize=3)
    ax[1].axhline(0, color="k", lw=0.8); ax[1].set_xticks(range(len(pts))); ax[1].set_xticklabels([f"{p[0][:3]}\n{p[1][:4]}" for p in pts], fontsize=6); ax[1].set(ylabel="loss-time error pred − CFD [d] (± CFD U)", title="contact-loss timing"); save(fig, out, "fig6_prediction_errors")
if __name__ == "__main__":
    ev, spec, out = sys.argv[1:4]; os.makedirs(out, exist_ok=True)
    for f, a in ((fig1, (out,)), (fig2, (out,)), (fig3, (ev, out)), (fig4, (ev, out)), (fig5, (ev, out)), (fig6, (os.path.join(RES, "analysis", "eval_analysis.json"), out))):
        try: f(*a); print(f.__name__, "ok")
        except Exception as e: import traceback; print(f.__name__, "FAILED", repr(e)); traceback.print_exc()
