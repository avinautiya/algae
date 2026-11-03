"""python3 make_figures.py  -> PNGs next to this file. Reads only committed result JSON."""
import json, os
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
H = os.path.dirname(os.path.abspath(__file__)); R = os.path.join(H, "..")
V = json.load(open(os.path.join(R, "viability", "viability_check.json"))); F = json.load(open(os.path.join(R, "verification", "refinement.json"))); A = json.load(open(os.path.join(R, "analytic_lane", "analytic_lane_results.json")))
# 1 viability
fig, ax = plt.subplots(1, 2, figsize=(11, 3.8)); cs = V["cases"]; ids = [c["id"] for c in cs]
ax[0].bar(range(len(cs)), [c["bound_21d_m"] * 1000 for c in cs], color="#1f77b4", label="rigorous 21-day upper bound on repair (mm)"); ax[0].axhline(cs[0]["required_m"] * 1000, color="#d62728", ls="--", label="required repair 61.1 mm"); ax[0].set_yscale("log"); ax[0].set_xticks(range(len(cs))); ax[0].set_xticklabels(ids, rotation=70, fontsize=6); ax[0].set_ylabel("gap repair over 21 d (mm, log)"); ax[0].legend(fontsize=7); ax[0].set_title("Bound vs requirement, 12 declared cases\n(not an extrapolated initial rate; ideal equal-K model)", fontsize=9)
n = V["ntu_dependence"]; ax[1].loglog([x["NTU_bypass_at_half_flow"] for x in n], [x["cap_over_film"] for x in n], "o-", label="bound/film limit = 1-(1-e^-N)/N"); ax[1].axhline(V["required_cap_over_film_for_21d"], color="#d62728", ls="--", label="ratio needed for 21 d (N=%.2f)" % V["required_NTU_for_21d"]); ax[1].set_xlabel("NTU of bypass route at half flow, N = hA/(ρc·Q/2)"); ax[1].set_ylabel("bound / film-limit rate"); ax[1].legend(fontsize=7); ax[1].set_title("Dependence on the assumption that controls the result\n(analysis only: declared parameters are not changed)", fontsize=9)
plt.tight_layout(); plt.savefig(os.path.join(H, "fig1_viability.png"), dpi=130); plt.close()
# 2 conservation / refinement
fig, ax = plt.subplots(1, 3, figsize=(12, 3.6))
for a, (name, key, lab) in zip(ax, [("space", "dx", "dx (m)"), ("outer_step", "dt", "outer step dt (s)"), ("inner_cfl", "cfl", "inner CFL")]):
    rows = F["controls"][name]; x = [r[key] for r in rows]
    a.semilogy(x, [max(r["max_over_steps"]["mass_kg"], 1e-18) for r in rows], "o-", label="max mass residual over steps (kg)"); a.semilogy(x, [max(r["max_over_steps"]["thermal_latent_J"], 1e-18) for r in rows], "s-", label="max sens.+latent residual over steps (J)")
    a.semilogy(x, [max(abs(r["geometry_minus_heat_m3"]), 1e-21) for r in rows], "^-", label="|geometric − heat-implied melt| (m³)"); a.set_xlabel(lab); a.set_title("control: " + name, fontsize=9); a.invert_xaxis() if key != "cfl" else None
ax[0].legend(fontsize=6); ax[0].set_ylabel("absolute defect")
fig.text(.5, -.03, "Per-step maxima, not only final cancellation. Short supported bare-ice case (20 s). Floor ~1e-18 plotted for exact zeros. Not event-time accuracy.", ha="center", fontsize=7)
plt.tight_layout(); plt.savefig(os.path.join(H, "fig2_conservation_refinement.png"), dpi=130, bbox_inches="tight"); plt.close()
# 2b convergence of outputs (relative difference to the finest level of each control)
fig, ax = plt.subplots(1, 3, figsize=(12, 3.4))
for a, (name, key) in zip(ax, [("space", "dx"), ("outer_step", "dt"), ("inner_cfl", "cfl")]):
    rows = F["controls"][name][:-1]; x = [r[key] for r in rows]
    a.plot(x, [100 * r["rel_diff_Q_toe_m3_s_to_finest"] for r in rows], "o-", label="toe discharge"); a.plot(x, [100 * r["rel_diff_melt_ice_ledger_m3_to_finest"] for r in rows], "s--", color="#d62728", label="ice melt volume")
    a.set_yscale("symlog", linthresh=1e-4); a.axhline(0, color="k", lw=.5); a.set_xlabel(key + " (finest level is the reference, not shown)"); a.set_title(name, fontsize=9)
ax[0].set_ylabel("relative difference to finest level (%)"); ax[0].legend(fontsize=7); ax[0].set_title("space: toe discharge differs by tens of % (not converged);\nmelt volume <1 %", fontsize=8)
fig.text(.5, -.03, "Short 20 s bare-ice case; engineering differences vs finest level, not confidence intervals and not event-time convergence.", ha="center", fontsize=7)
plt.tight_layout(); plt.savefig(os.path.join(H, "fig3_output_refinement.png"), dpi=130, bbox_inches="tight"); plt.close()
# 4 endpoints
ev = [c for c in A["cases"] if c["control"] == "evolving"]; fig, ax = plt.subplots(1, 2, figsize=(11, 3.8))
col = {"initially_disconnected": "#999999", "lost_without_reconnection": "#ff7f0e", "sustained_postpulse_recovery": "#2ca02c"}
for c in ev: ax[0].scatter(c["pulse_multiple"] * {6: 0.82, 24: 1.0, 72: 1.22}[c["pulse_h"]], c["toe_area_m2"], color=col[c["classification"]], s=60, marker={6: "v", 24: "o", 72: "^"}[c["pulse_h"]])
ax[0].set_xscale("log"); ax[0].set_xlabel("pulse multiple of Q0 (markers: ▼6 h, ● 24 h, ▲72 h)"); ax[0].set_ylabel("toe area (m²)"); ax[0].set_title("Endpoint class (grey initially disconnected, orange rewet then lost,\nno green = no sustained recovery in any of 36 cases)", fontsize=8)
ax[1].scatter([c["gap_change_m"] * 1000 for c in ev], [c["peak_toe_fraction_in_pulse"] for c in ev], c="#1f77b4"); ax[1].set_xlabel("gap change at end of 48 h post-pulse window (mm; + = gap grew)"); ax[1].set_ylabel("peak toe flow fraction during pulse"); ax[1].set_title("Temporary rewetting with no repair (all gap changes ≥ 0)", fontsize=8)
fig.text(.5, -.03, "Analytic ODE lane of the declared ideal model; inspected development parameters; not blind, not CFD, not field.", ha="center", fontsize=7); plt.tight_layout(); plt.savefig(os.path.join(H, "fig4_endpoints.png"), dpi=130, bbox_inches="tight"); plt.close(); print("ok")
