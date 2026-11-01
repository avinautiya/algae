"""Reproducible M2 figures from COMPLETED data only: python -m ice_cliff_survival.cfd.make_figures_m2  (reads cfd/results/m2/*.json, *.npz, shield_planar.json; writes cfd/figures_m2/*.png).
Labels: 'analytic' = closed-form/series reference; 'numerical' = flowStefanFoam / FV solution; 'assumed' = prescribed test inputs."""
import json, os, math
import numpy as np
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__)); RES = os.path.join(HERE, "results", "m2"); OUT = os.path.join(HERE, "figures_m2")
C = dict(a="#1f5f8b", b="#c8641e", c="#3a8f5b", g="#555555")
plt.rcParams.update({"font.size": 9, "axes.spines.top": False, "axes.spines.right": False, "axes.grid": True, "grid.alpha": 0.25, "figure.dpi": 130})
J = lambda n: json.load(open(os.path.join(RES, n + ".json")))
def fig_wedge():
    from ..geometry import debris_normal_thickness, debris_normal_thickness_v1_projected_horizontal as old
    from .. import shield_planar as S
    fig, ax = plt.subplots(1, 3, figsize=(11, 3.3)); zf, th, be = 0.012, 62.0, 33.0; z = np.linspace(0, zf, 100)
    ax[0].plot(z * 1e3, [debris_normal_thickness(v, zf, th, be) * 1e3 for v in z], color=C["a"], label="M2 corrected (ray-line)"); ax[0].plot(z * 1e3, [old(v, zf, th, be) * 1e3 for v in z], color=C["b"], ls="--", label="M1 (deprecated)")
    ax[0].set(xlabel="height z on face [mm]", ylabel="normal thickness d [mm]", title="debris normal thickness (θ=62°, β=33°, z_f=12 mm)"); ax[0].legend(frameon=False)
    h = lambda f: [0.5 / max(f(v, zf, th, be), 1e-4) for v in z[:-3]]
    ax[1].plot(z[:-3] * 1e3, h(debris_normal_thickness), color=C["a"], label="corrected"); ax[1].plot(z[:-3] * 1e3, h(old), color=C["b"], ls="--", label="M1"); ax[1].set(xlabel="z [mm]", ylabel="h = k/d [W m⁻² K⁻¹]", title="shield conductance (k_deb=0.5 assumed)"); ax[1].legend(frameon=False)
    rows = json.load(open(os.path.join(HERE, "results", "shield_planar.json")))["rows"]
    for d, col in ((0.002, C["c"]), (0.007, C["a"]), (0.02, C["b"]), (0.1, C["g"])):
        t = np.logspace(0, 5, 60); Tw = 1.0; k, rc = 0.5, 2.5e6; Q = S.Q_series(t, d, k, rc, Tw, nterms=2000); ax[2].loglog(t, 1 - Q / S.Q_massless(t, d, k, Tw), color=col, label=f"d={d*1e3:g} mm")
    ax[2].set(xlabel="t [s]", ylabel="relative heat deficit vs massless (analytic)", title="resolved debris storage (ρc=2.5e6)", ylim=(1e-4, 1.5)); ax[2].legend(frameon=False, fontsize=7); fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig1_wedge_shield.png")); plt.close(fig)
def fig_benchmarks():
    r1 = J("r1_graetz")["results"]; r2 = J("r2_stefan")["results"]; fig, ax = plt.subplots(1, 3, figsize=(11, 3.2))
    ny = [r["ny"] for r in r1]; ax[0].semilogy(ny, [r["err_Tb_rel_max"] for r in r1], "o-", color=C["a"], label="bulk T"); ax[0].semilogy(ny, [r["err_qwall_rel_max"] for r in r1], "s-", color=C["b"], label="wall flux"); ax[0].semilogy(ny, [r["err_u_rms_over_Ub"] for r in r1], "^-", color=C["c"], label="velocity profile")
    ax[0].set(xlabel="rows ny", ylabel="rel. error vs analytic Graetz", title="no melting: Graetz–Poiseuille (numerical vs analytic)"); ax[0].legend(frameon=False)
    for rho, col, m in ((1000.0, C["a"], "o"), (917.0, C["b"], "s")):
        sel = [r for r in r2 if r["rho_i"] == rho and r["dt"] == 40.0]; ax[1].loglog([r["ny"] for r in sel], [abs(r["rel_err_end"]) for r in sel], m + "-", color=col, label=f"ρ_i={rho:g}")
        sel = [r for r in r2 if r["rho_i"] == rho and r["ny"] == 40]; ax[2].loglog([r["dt"] for r in sel], [abs(r["rel_err_end"]) for r in sel], m + "-", color=col, label=f"ρ_i={rho:g}")
    ax[1].set(xlabel="rows ny (dt=40 s)", ylabel="|s−s_exact|/s_exact at 1 h", title="Stefan limit vs Neumann similarity (analytic)"); ax[2].set(xlabel="dt [s] (ny=40)", ylabel="|s−s_exact|/s_exact", title="Stefan limit: time step"); ax[1].legend(frameon=False); ax[2].legend(frameon=False)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig2_benchmarks.png")); plt.close(fig)
def fig_fields():
    d = np.load(os.path.join(RES, "fields_r4_917.npz")); T, Ux, Uy, b, Dx, H = d["T"], d["Ux"], d["Uy"], d["b"], float(d["Dx"]), float(d["H"]); nx, ny = T.shape
    X = np.zeros((nx + 1, ny + 1)); Y = np.zeros((nx + 1, ny + 1))
    for i in range(nx + 1):
        X[i] = i * Dx; Y[i] = b[i] + (H - b[i]) * np.linspace(0, 1, ny + 1)
    fig, ax = plt.subplots(3, 1, figsize=(9.5, 6.2), sharex=True)
    m = ax[0].pcolormesh(X * 100, Y * 1e3, T, cmap="viridis", shading="flat"); fig.colorbar(m, ax=ax[0], label="T [K above T_m]"); ax[0].set(ylabel="y [mm]", title="numerical: temperature at t = %.0f s (ρ_i = 917, bare ice)" % float(d["t"]))
    sp = np.hypot(Ux, Uy); m = ax[1].pcolormesh(X * 100, Y * 1e3, sp, cmap="magma", shading="flat"); fig.colorbar(m, ax=ax[1], label="|u| [m/s]"); ax[1].set(ylabel="y [mm]", title="numerical: speed")
    xs = np.arange(nx + 1) * Dx * 100; ax[2].plot(xs, b * 1e3, color=C["a"]); ax[2].set(xlabel="x [cm]", ylabel="floor position b [mm]", title="numerical: interface geometry (initially flat at 0; inert walls 0–8 cm and 28–34 cm)")
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig3_coupled_fields.png")); plt.close(fig)
def fig_budgets():
    r4 = J("r4_coupled")["results"]; r8 = J("r8_shield"); labs = ["discrete\nenergy", "physical\nenergy", "interface heat\nvs latent", "independent\nenergy (cum.)", "GCL\n(vol/V_sw)", "coupling\ndefect/step", "local\nmass (solver)", "local mass\n(independent)", "global mass\n(independent)"]
    fig, ax = plt.subplots(figsize=(10, 3.4)); w = 0.38
    for k, (r, col) in enumerate(zip(r4, (C["a"], C["b"]))):
        v = [abs(r["R_disc_over_scale"]), abs(r["R_phys_over_scale"]), abs(r["interface_heat_vs_latent_rel"]), r["indep_energy_cum_resid_rel"], r["gcl_max_over_Vsw"], r["max_coupling_defect_m"] / 1e-3, r["massLocal_solver"], r["indep_local_mass_resid_max_rel"], r["indep_global_mass_resid_rel_max"]]
        ax.bar(np.arange(len(v)) + (k - 0.5) * w, np.maximum(v, 1e-18), w, color=col, label=f"ρ_i = {r['rho_i']:g}")
    ax.set_yscale("log"); ax.set_xticks(range(len(labs))); ax.set_xticklabels(labs, fontsize=7); ax.set(ylabel="normalised residual (log; floor 1e-18)", title="coupled throughflow: budget residuals (solver ledger vs independent reconstruction)"); ax.legend(frameon=False)
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig4_budgets.png")); plt.close(fig)
def fig_convergence():
    r = J("r6_refine"); fig, ax = plt.subplots(1, 3, figsize=(11, 3.2)); s = r["space"]; ref = s[-1]
    h = [m["dx"] for m in s]; ax[0].loglog(h[:-1], [abs(m["Vsw"] - ref["Vsw"]) / ref["Vsw"] for m in s[:-1]], "o-", color=C["a"], label="melted volume"); ax[0].loglog(h[:-1], [abs(m["T_out"] - ref["T_out"]) / ref["T_out"] for m in s[:-1]], "s-", color=C["b"], label="outlet bulk T")
    ax[0].set(xlabel="Δx [m] (ny scaled)", ylabel="|f − f_finest|/f_finest", title="flow/thermal mesh (numerical)"); ax[0].legend(frameon=False)
    t = r["time"]; reft = t[-1]; dts = [m["dt"] for m in t]; ax[1].loglog(dts[:-1], [abs(m["Vsw"] - reft["Vsw"]) / reft["Vsw"] for m in t[:-1]], "o-", color=C["a"]); ax[1].set(xlabel="Δt [s]", ylabel="|V−V_Δt=2.5|/V", title="physical time step")
    c = r["coupling"]; ax[2].semilogx([m["tol"] if m["maxCoupling"] > 1 else 1.0 for m in c], [abs(m["Vsw"] - c[-1]["Vsw"]) / c[-1]["Vsw"] for m in c], "o-", color=C["a"]); ax[2].set(xlabel="coupling tolerance (1.0 = explicit predictor only)", ylabel="|V−V_tol=1e-9|/V", title="interface coupling"); ax[2].invert_xaxis()
    fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig5_convergence.png")); plt.close(fig)
def fig_cost():
    r = J("r6_refine"); fig, ax = plt.subplots(figsize=(4.8, 3.3)); s = r["space"]
    ax.loglog([m["cells"] for m in s], [m["wall_s"] for m in s], "o-", color=C["a"]); ax.set(xlabel="cells", ylabel="wall time [s] (dt=2.5 s, 1800 s)", title="measured cost (single core)"); fig.tight_layout(); fig.savefig(os.path.join(OUT, "fig6_cost.png")); plt.close(fig)
if __name__ == "__main__":
    os.makedirs(OUT, exist_ok=True)
    for f in (fig_wedge, fig_benchmarks, fig_fields, fig_budgets, fig_convergence, fig_cost):
        try: f(); print(f.__name__, "ok")
        except Exception as e: print(f.__name__, "FAILED", repr(e))
