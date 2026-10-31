"""Figures of the V2 milestone-2 package, generated ONLY from stored result files (no simulation). python -m ice_network_v2.make_figures_v2 (from the repo root)"""
import json, os
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
H = os.path.dirname(os.path.abspath(__file__)); J = lambda p: json.load(open(os.path.join(H, p)))
C = dict(a="#1b6ca8", b="#d95f02", c="#2e8b57", d="#7b3294", g="#555555"); plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.alpha": 0.2})
def slope(ax, x0, y0, p, x1, label):
    xs = np.array([x0, x1]); ax.plot(xs, y0 * (xs / x0) ** p, color=C["g"], ls="--", lw=0.8); ax.text(x1, y0 * (x1 / x0) ** p, label, color=C["g"], fontsize=8, va="bottom")
# Fig 1: laminar verification errors (separate) and thermal-budget diagnosis
g = J("cfd/graetz_final.json")["space"]; ba = J("cfd/budget_audit.json"); fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
dy = np.array([r["dy_over_H"] for r in g])
for k, lab, c in (("err_u_over_Ub", "velocity profile (÷Ub)", C["a"]), ("err_T_shape_rms", "temperature profile shape (÷bulk excess)", C["b"]), ("err_wall_flux_rel", "wall heat flux (÷reference)", C["c"])): ax[0].loglog(dy, [r[k] for r in g], "o-", color=c, label=lab)
slope(ax[0], 0.1, 3e-3, 2, 0.025, "slope 2"); ax[0].set_xlabel("cell size dy / H"); ax[0].set_ylabel("error at x = 18 H"); ax[0].legend(frameon=False, fontsize=8); ax[0].set_title("Poiseuille/Graetz (OpenFOAM v1912), errors separately", fontsize=9)
n = [8, 16, 32]; ax[1].semilogy(n, [abs(ba[f"step_n{i}"]["inlet_cond_over_Pin"]) for i in n], "o-", color=C["a"], label="inlet axial conduction ÷ P_in (step)")
ax[1].semilogy(n, [abs(ba[f"step_n{i}"]["R_over_Pin"]) + 1e-13 for i in n], "s-", color=C["c"], label="residual after including it: step (orthogonal mesh)"); ax[1].semilogy(n, [abs(ba[f"ramp_n{i}"]["R_over_Pin"]) for i in n], "^-", color=C["b"], label="residual, ramp (non-orthogonal, max 63°)")
ax[1].set_xlabel("cells per H"); ax[1].set_ylabel("|budget residual| ÷ P_in"); ax[1].legend(frameon=False, fontsize=7.5); ax[1].set_title("Thermal budget: growth = inlet-corner conduction", fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(H, "figures", "fig1_laminar_and_budget.png"), dpi=150); plt.close(fig)
# Fig 2: Stefan in OpenFOAM ALE solver
s = J("foam/stefan1DFoam/stefan_of_results.json"); fig, ax = plt.subplots(1, 2, figsize=(10, 3.8))
for (Ste, tag), c, mk in ((((0.00631, "rhoI917_jump")), C["a"], "o"), ((0.00631, "equal_density"), C["b"], "s"), ((0.2, "rhoI917_jump"), C["c"], "^"), ((0.2, "equal_density"), C["d"], "v")):
    nx = [20, 40, 80, 160]; ax[0].loglog(nx, [s[f"Ste{Ste}_{tag}_nx{i}_dt2.5"]["s_rel_err"] for i in nx], mk + "-", color=c, label=f"Ste {Ste}, {tag.replace('_', ' ')}")
    dt = [2.5, 5.0, 10.0]; ax[1].semilogx(dt, [s[f"Ste{Ste}_{tag}_nx40_dt{d}"]["s_rel_err"] / s[f"Ste{Ste}_{tag}_nx40_dt2.5"]["s_rel_err"] for d in dt], mk + "-", color=c)
slope(ax[0], 20, 1e-3, -1, 160, "slope −1"); ax[0].set_xlabel("cells nx (dt = 2.5 s)"); ax[0].set_ylabel("front position rel. error at t = 3600 s"); ax[0].legend(frameon=False, fontsize=7); ax[0].set_title("stefan1DFoam vs exact Stefan solution (density jump included)", fontsize=9)
ax[1].set_xlabel("time step dt [s] (nx = 40)"); ax[1].set_ylabel("front error ÷ error at dt = 2.5 s"); ax[1].set_title("time-step refinement: front error unchanged", fontsize=9); ax[1].set_ylim(0.8, 1.2)
fig.tight_layout(); fig.savefig(os.path.join(H, "figures", "fig2_stefan_openfoam.png"), dpi=150); plt.close(fig)
# Fig 3: free-surface verification
hs = J("foam/freesurface/hydrostatic_results.json"); th = J("foam/freesurface/thermal_fs_results.json"); fig, ax = plt.subplots(1, 3, figsize=(12, 3.6))
dyh = [hs[k]["dx"] for k in hs]; ax[0].loglog(dyh, [hs[k]["2"]["max_U_water"] for k in hs], "o-", color=C["a"], label="water, t = 2 s"); ax[0].loglog(dyh, [hs[k]["2"]["max_U"] for k in hs], "s-", color=C["b"], label="air, t = 2 s")
ax[0].set_xlabel("cell size [m]"); ax[0].set_ylabel("max |U| at rest [m/s]"); ax[0].legend(frameon=False); ax[0].set_title("A. hydrostatic rest, stepped bed", fontsize=9)
d1 = [0.05 / n for n in (10, 20, 40)]; ax[1].loglog(d1, [abs(th[f"T1_n{n}_co0.2"]["amp_rel_err"]) for n in (10, 20, 40)], "o-", color=C["c"], label="cosine-mode amplitude"); ax[1].loglog(d1, [th[f"T1_n{n}_co0.2"]["profile_rms_err"] for n in (10, 20, 40)], "s-", color=C["d"], label="profile RMS (K)")
slope(ax[1], 0.005, 1.4e-3, 2, 0.00125, "slope 2"); ax[1].set_xlabel("dy [m]"); ax[1].set_title("T1. diffusion in water at rest, adiabatic walls", fontsize=9); ax[1].legend(frameon=False, fontsize=8)
d2 = [0.5 / n for n in (50, 100, 200)]; ax[2].loglog(d2, [th[f"T2_n{n}_co0.2"]["L2_err_water"] for n in (50, 100, 200)], "o-", color=C["a"], label="Co 0.2"); ax[2].loglog([0.005], [th["T2_n100_co0.1"]["L2_err_water"]], "s", color=C["b"], label="Co 0.1 (n=100)")
slope(ax[2], 0.01, 8e-2, 1, 0.0025, "slope 1"); ax[2].set_xlabel("dx [m]"); ax[2].set_ylabel("L2 error (water, K)"); ax[2].legend(frameon=False, fontsize=8); ax[2].set_title("T2. advection of a pulse, water/air cyclic", fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(H, "figures", "fig3_freesurface_verification.png"), dpi=150); plt.close(fig)
# Fig 4: rigid-lid step vs ramp wall heat flux
sp = J("cfd/step_ramp_post.json"); fig, ax = plt.subplots(1, 2, figsize=(10, 3.6))
for n, ls in ((8, ":"), (16, "--"), (32, "-")):
    ax[0].plot(sp[f"step_n{n}"]["x"], sp[f"step_n{n}"]["qn"], ls, color=C["a"], label=f"step n={n}"); ax[0].plot(sp[f"ramp_n{n}"]["x"], sp[f"ramp_n{n}"]["qn"], ls, color=C["b"], label=f"ramp n={n}")
ax[0].set_xlim(-12, 30); ax[0].set_ylim(0, 3); ax[0].set_xlabel("x / H"); ax[0].set_ylabel("wall heat flux ÷ (α ΔT / H)"); ax[0].legend(frameon=False, fontsize=7, ncol=2); ax[0].set_title("rigid-lid laminar, Re=50: wall heat flux", fontsize=9)
Tb = [(sp[f"step_n{n}"]["Tb_out"], sp[f"ramp_n{n}"]["Tb_out"]) for n in (8, 16, 32)]; ax[1].plot([8, 16, 32], [t[0] for t in Tb], "o-", color=C["a"], label="step"); ax[1].plot([8, 16, 32], [t[1] for t in Tb], "s-", color=C["b"], label="ramp")
ax[1].set_xlabel("cells per H"); ax[1].set_ylabel("outlet flow-weighted T_b"); ax[1].legend(frameon=False); ax[1].set_title("downstream enthalpy: step–ramp gap (6e-4) < mesh change (2e-3)", fontsize=9)
fig.tight_layout(); fig.savefig(os.path.join(H, "figures", "fig4_step_vs_ramp_rigidlid.png"), dpi=150); plt.close(fig)
# Fig 5: bore verification (if available)
bp = os.path.join(H, "foam/freesurface/bore_results.json")
if os.path.exists(bp):
    b = J("foam/freesurface/bore_results.json"); fig, ax = plt.subplots(1, 2, figsize=(10, 3.6)); ok = {k: v for k, v in b.items() if "failed" not in v and "not_run" not in v}
    for (k, v), c in zip(ok.items(), (C["a"], C["b"], C["c"], C["d"], C["g"])): x = np.linspace(0, 4, len(v["depth_profile"])); ax[0].plot(x, v["depth_profile"], color=c, label=f"dy = {v['dy']*1e3:.1f} mm ({k.split('_')[1]})")
    ax[0].axhline(0.05, color="k", lw=0.6); ax[0].axhline(0.1536, color="k", lw=0.6, ls="--"); ax[0].set_xlabel("x [m]"); ax[0].set_ylabel("depth [m] at t = 1.5 s"); ax[0].legend(frameon=False, fontsize=7); ax[0].set_title("bore depth profiles (dashed: design conjugate depth)", fontsize=9)
    ax[1].plot([v["dy"] * 1e3 for v in ok.values()], [v["h_rel_err"] for v in ok.values()], "o-", color=C["a"], label="depth behind vs conjugate(measured Fr)"); ax[1].plot([v["dy"] * 1e3 for v in ok.values()], [v["momentum_flux_rel_diff"] for v in ok.values()], "s-", color=C["b"], label="momentum-flux imbalance"); ax[1].plot([v["dy"] * 1e3 for v in ok.values()], [v["mass_flux_rel_diff"] for v in ok.values()], "^-", color=C["c"], label="mass-flux imbalance")
    ax[1].axhline(0, color="k", lw=0.6); ax[1].set_xlabel("cell size [mm]"); ax[1].set_ylabel("relative error"); ax[1].legend(frameon=False, fontsize=7); ax[1].set_title("bore jump conditions vs resolution (unresolved)", fontsize=9)
    fig.tight_layout(); fig.savefig(os.path.join(H, "figures", "fig5_bore.png"), dpi=150); plt.close(fig)
print("figures written")
