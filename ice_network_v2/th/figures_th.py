"""Figures of the matched-state experiment, generated only from results/*.json, analysis.json, conductance_test.json."""
import json, os, sys
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); import baselines as BL
R = lambda n: json.load(open(os.path.join(HERE, "results", n + ".json"))); a = json.load(open(os.path.join(HERE, "analysis.json"))); ct = json.load(open(os.path.join(HERE, "conductance_test.json")))
H = 0.01; AL = 1e-6 / 7.0; F = H / AL; plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.alpha": 0.2}); CA, CB, CG, CO = "#1b6ca8", "#d95f02", "#555555", "#2e8b57"
st = R("step_n32_base"); sr = R("straight_n32_base"); A = st["A_dt10"]; B = st["B_dt10"]; s = A["Tb_sec"]["-1.0"] / B["Tb_sec"]["-1.0"]
# 1 matched profiles
fig, ax = plt.subplots(1, 3, figsize=(11, 3.6)); yA, TA, uA = map(np.array, A["prof"]["-1.0"]); yB, TB, uB = map(np.array, B["prof"]["-1.0"]); TbA = A["Tb_sec"]["-1.0"]
ax[0].plot(uA / np.mean(uA), yA - 1, color=CG, label="velocity (identical for A and B)"); ax[0].set_xlabel("u / mean u"); ax[0].set_ylabel("height above bed [H]"); ax[0].legend(frameon=False, fontsize=7)
ax[1].plot(TA / TbA, yA - 1, color=CA, label="A (cold wall from −3 H)"); ax[1].plot(s * TB / TbA, yB - 1, color=CB, label="B* (from −11 H, matched)"); 
for hs, dist, c in (("A", 2.0, CA), ("B", 10.0, CB)):
    Tb_g, qw_g, th, eta = BL.graetz_section(dist); ax[1].plot(th / Tb_g, eta, ":", color=c, lw=1)
ax[1].set_xlabel("T / T_b,A at x = −1 H (dotted: Graetz series)"); ax[1].legend(frameon=False, fontsize=7); ax[1].set_title("matched bulk T_b: profiles differ", fontsize=9)
ax[2].plot(A["Tb_sec"].keys() and [A["Tb_sec"][k] for k in ("-2.0", "-1.0", "-0.5")], [-2, -1, -0.5], "o-", color=CA, label="T_b(x) A"); ax[2].plot([s * B["Tb_sec"][k] for k in ("-2.0", "-1.0", "-0.5")], [-2, -1, -0.5], "s-", color=CB, label="s·T_b(x) B"); ax[2].set_xlabel("flow-weighted T_b"); ax[2].set_ylabel("x / H"); ax[2].legend(frameon=False, fontsize=7); ax[2].set_title("bulk matched at x = −1 H only", fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "figures", "th1_matched_profiles.png"), dpi=150); plt.close(fig)
# 2 heat flux
fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
for d, nm, ls in ((st, "step", "-"), (sr, "straight", "--")):
    A_, B_ = d["A_dt10"], d["B_dt10"]; sc = A_["Tb_sec"]["-1.0"] / B_["Tb_sec"]["-1.0"]; x = np.array(A_["x"]); m = x >= -1e-9; qa = np.array(A_["q"])[m] * F; qb = np.array(B_["q"])[m] * F * sc
    ax[0].plot(x[m], qa, ls, color=CA, label=f"{nm} A"); ax[0].plot(x[m], qb, ls, color=CB, label=f"{nm} B*"); ax[1].plot(x[m], qb / qa, ls, color=CO if nm == "step" else CG, label=nm)
ax[0].set_xlim(0, 30); ax[0].set_xlabel("x / H (bed, x ≥ 0)"); ax[0].set_ylabel("q ÷ (α/H) (T_b,A = 1 at section scale)"); ax[0].legend(frameon=False, fontsize=7, ncol=2); ax[0].set_title("n = 32 wall heat flux (step face at x = 0)", fontsize=9)
ax[1].set_xlim(0, 30); ax[1].set_xlabel("x / H"); ax[1].set_ylabel("q_B* / q_A"); ax[1].legend(frameon=False); ax[1].set_title("history effect on local flux", fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "figures", "th2_heatflux_step_straight.png"), dpi=150); plt.close(fig)
# 3 decision figure
fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
for j, X in enumerate((10, 30)):
    d = a["decision"][f"X{X}H"]; vals = [("step", d["delta_step"], a["step"]["U"][f"delta{X}"]["U_frozen"], a["step"]["U"][f"delta{X}"]["U_without_B_variants"]), ("straight", d["delta_straight"], a["straight"]["U"][f"delta{X}"]["U_frozen"], a["straight"]["U"][f"delta{X}"]["U_without_B_variants"])]
    for i, (nm, v, U1, U2) in enumerate(vals): ax[j].errorbar([i - 0.08], [v], yerr=[U1], fmt="o", color=CA, capsize=4, label="frozen U (incl. B-history spread)" if i == 0 else None); ax[j].errorbar([i + 0.08], [v], yerr=[U2], fmt="s", color=CB, capsize=4, label="numerical-only U" if i == 0 else None)
    ax[j].axhline(0, color=CG, lw=0.8); ax[j].text(1.45, 0.005, "L predicts 0", color=CG, fontsize=8, ha="right"); ax[j].set_xticks([0, 1]); ax[j].set_xticklabels(["step", "straight (baseline D)"]); ax[j].set_ylabel(f"δ_{X}H = (Q_A − Q_B*)/mean"); ax[j].set_title(f"bed heat 0–{X} H; step–straight excess {d['excess']*100:.1f} pp", fontsize=9)
ax[0].legend(frameon=False, fontsize=7); fig.tight_layout(); fig.savefig(os.path.join(HERE, "figures", "th3_history_effect_decision.png"), dpi=150); plt.close(fig)
# 4 convergence
fig, ax = plt.subplots(1, 2, figsize=(10, 3.6)); ns = [8, 16, 32]
for kind, c in (("step", CA), ("straight", CB)):
    for k, ls in (("delta10", "-"), ("delta30", "--")): ax[0].plot(ns, a[kind][f"mesh_seq_{k}"], "o" + ls, color=c, label=f"{kind} {k}")
ax[0].set_xscale("log", base=2); ax[0].set_xlabel("cells per H"); ax[0].set_ylabel("δ"); ax[0].legend(frameon=False, fontsize=7); ax[0].set_title("mesh sequence (step: not monotone-converged)", fontsize=9)
comp = ["u_mesh", "u_dt", "u_iter", "u_match", "u_loc_historyB", "u_loc_geometry"]; x = np.arange(len(comp))
ax[1].bar(x - 0.2, [a["step"]["U"]["delta30"][c] for c in comp], 0.4, color=CA, label="step δ30"); ax[1].bar(x + 0.2, [a["straight"]["U"]["delta30"][c] for c in comp], 0.4, color=CB, label="straight δ30"); ax[1].set_yscale("log"); ax[1].set_xticks(x); ax[1].set_xticklabels([c.replace("u_", "") for c in comp], rotation=30, fontsize=8); ax[1].legend(frameon=False, fontsize=7); ax[1].set_title("error-budget components (log)", fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "figures", "th4_convergence_budget.png"), dpi=150); plt.close(fig)
# 5 baselines and conductance shape
fig, ax = plt.subplots(1, 2, figsize=(10, 3.6)); L = a["local_closure"]; labs = []; vals = []
for kind in ("step", "straight"):
    for X in (10, 30): labs += [f"{kind}\n0–{X} H"]; vals += [(L[kind][f"relerr{X}_L_A"], L[kind][f"relerr{X}_L_B"])]
x = np.arange(len(labs)); ax[0].bar(x - 0.2, [v[0] for v in vals], 0.4, color=CA, label="L vs A"); ax[0].bar(x + 0.2, [v[1] for v in vals], 0.4, color=CB, label="L vs B*"); ax[0].set_xticks(x); ax[0].set_xticklabels(labs, fontsize=8); ax[0].axhline(0, color=CG, lw=0.8); ax[0].set_ylabel("(Q_L − Q_resolved)/Q_resolved"); ax[0].legend(frameon=False, fontsize=7); ax[0].set_title("local fully developed closure vs resolved bed heat", fontsize=9)
for k, c, ls in (("step_n32", CA, "-"), ("straight_n32", CB, "--")): mid = np.array(ct[k]["bin_edges_H"]); mid = 0.5 * (mid[1:] + mid[:-1]); ax[1].plot(mid, ct[k]["ratio_B_over_A"], "o" + ls, color=c, label=k)
ax[1].set_xlabel("x / H (3 H bins)"); ax[1].set_ylabel("conductance ratio g_B,i / g_A,i"); ax[1].legend(frameon=False); ax[1].set_title("common-shape assumption: ratio should be constant", fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(HERE, "figures", "th5_baselines_conductance.png"), dpi=150); plt.close(fig); print("ok")
