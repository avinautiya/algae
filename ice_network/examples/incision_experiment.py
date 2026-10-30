"""Controlled warm-water incision experiment: matched geometry, forcing and constants; four treatments.
T0 no channel melting | T1 thermally coupled (water-temperature transport, wetted-perimeter partition) | T2 explicit dissipation-to-melt
(all friction+entrance dissipation melts ice where released; walls/bed share by wetted perimeter) | T3 SaDS-type approximation (as T2, but the
section grows at a fixed width-to-depth ratio, hydraulics still Manning/GVF -- NOT a reimplementation of SaDS).
E1: fraction of channel dissipation that reaches the ice versus thermal exposure (T1 only), against the closed form 1-(1-e^-a)/a.
E2: incision under warm inflow; E3: numerical refinement and heat-transfer multiplier sensitivity of T1."""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src")); sys.path.insert(0, os.path.join(HERE, "..", "tests"))
from helpers import chan, const, lake
from ice_network.constants import CP_W, LF, RHO_I, RHO_W
from ice_network.evolution.network import Network
from ice_network.thermal.channel_heat import heat_transfer_coeff

OUT = os.path.join(HERE, "..", "records", "incision_experiment.json")
res = {"E1": [], "E2": [], "E3": {}}


def single(L, slope, n, w, T0, A, I, dt, t_end, stage_above=0.2, **kw):
    e = chan(0, 1, 100.5, slope * L, w, L=L, n=n)
    net = Network([lake(99.0, A)], [dict(h_fixed=95.0 - slope * L)], [e], **kw)
    v0 = A * (100.5 + stage_above - 99.0)
    r = net.run([v0], [T0], t_end, dt, const(I), const(T0), rtol=1e-10, atol=1e-11, evolve=kw.get("melting", True))
    return net, r


# ------------------------------------------------------------------ E1 partition versus thermal exposure (geometry frozen)
for L in (25.0, 50.0, 100.0, 200.0, 400.0):
    for nu in (0.25, 1.0, 4.0):
        net, r = single(L, 0.005, 17, 0.4, 0.0, 1e6, 0.0, 1800.0, 2 * 3600.0, nu_factor=nu)
        b = r["budgets"]; s = r["hist"][len(r["hist"]) // 2]
        q = abs(s["Q"][0]); e = net.edges[0]
        st = net.edge_state(e, s["stage"][0], s["stage"][1], 0.0)
        h = st["h"].mean(); wmean = e["w"].mean()
        a = heat_transfer_coeff(q, wmean, h, nu) * (wmean + 2 * h) * L / (RHO_W * CP_W * q)
        eta = b["heat_to_ice_J"] / b["dissipation_channel_J"]
        res["E1"].append(dict(L_m=L, nu_factor=nu, Q_m3_s=q, exposure_a=a, ice_fraction_of_channel_dissipation=eta, closed_form_1_minus_F_over_a=1 - (1 - np.exp(-a)) / a,
                              dissipation_channel_J=b["dissipation_channel_J"], dissipation_exit_J=b["dissipation_exit_J"], heat_to_ice_J=b["heat_to_ice_J"],
                              sensible_out_J=b["sensible_out_J"], energy_residual_rel=b["energy_residual_rel"]))
# ------------------------------------------------------------------ E2 incision, four treatments
def treatment(name, T0, **kw):
    net, r = single(100.0, 0.005, 17, 0.4, T0, 2000.0, 0.05, 450.0, 8 * 3600.0, stage_above=0.05, **kw)
    b = r["budgets"]; e = net.edges[0]
    z0 = 100.5 - 0.005 * 100.0 * np.linspace(0, 1, 17)
    last = r["hist"][-1]
    return dict(treatment=name, T_inflow_C=T0, sill_lowering_m=100.5 - e["z"][0], max_bed_lowering_m=float(np.max(z0 - e["z"])), mean_width_gain_m=float(np.mean(e["w"] - 0.4)),
                delivered_m3=float(b["delivered_m3"][0]), outlet_Q_end_m3_s=float(abs(last["Q"][0])), lake_stage_end_m=float(last["stage"][0]),
                heat_in_sensible_J=b["sensible_in_J"], heat_stored_sensible_change_J=b["sensible_storage_change_J"], heat_to_ice_latent_J=b["heat_to_ice_J"],
                heat_out_sensible_J=b["sensible_out_J"], dissipation_total_J=b["dissipation_J"], dissipation_channel_J=b["dissipation_channel_J"], dissipation_exit_J=b["dissipation_exit_J"],
                ice_mass_lost_kg=b["ice_mass_from_heat_kg"], liquid_mass_gained_kg=b["ice_mass_from_meltwater_kg"],
                water_res=b["water_residual_rel"], mass_res=b["mass_residual_rel"], energy_res=b["energy_residual_rel"], geometry_vs_latent_rel_max=b["geometry_vs_latent_rel_max"])

for T0 in (0.0, 0.3, 1.0):
    res["E2"].append(treatment("T0_no_melting", T0, melting=False))
    res["E2"].append(treatment("T1_thermal", T0))
    res["E2"].append(treatment("T2_dissipation_all", T0, melt_mode="dissipation_all"))
    res["E2"].append(treatment("T3_dissipation_all_fixed_aspect", T0, melt_mode="dissipation_all", geometry_mode="fixed_aspect"))
# ------------------------------------------------------------------ E3 refinement and closure multiplier (T1, T0 = 1 C)
def t1(n, dt, nu=1.0, msub=16):
    e = chan(0, 1, 100.5, 0.5, 0.4, L=100.0, n=n, msub=msub)
    net = Network([lake(99.0, 2000.0)], [dict(h_fixed=94.5)], [e], nu_factor=nu)
    r = net.run([2000.0 * 1.55], [1.0], 8 * 3600.0, dt, const(0.05), const(1.0), rtol=1e-10, atol=1e-11)
    b = r["budgets"]; return dict(n=n, dt=dt, nu_factor=nu, sill_lowering_m=100.5 - net.edges[0]["z"][0], delivered_m3=float(b["delivered_m3"][0]), heat_to_ice_J=b["heat_to_ice_J"])
res["E3"]["stations"] = [t1(n, 450.0) for n in (5, 9, 17, 33, 65)]
res["E3"]["dt"] = [t1(17, dt) for dt in (1800.0, 900.0, 450.0, 225.0, 112.5)]
res["E3"]["closure_multiplier"] = [t1(17, 450.0, nu) for nu in (0.25, 0.5, 1.0, 2.0, 4.0)]
json.dump(res, open(OUT, "w"), indent=1, default=float)
print("E1 (L, nu, a, eta, closed form)")
for r in res["E1"]: print("  %5.0f %5.2f a=%.3f eta=%.3f cf=%.3f" % (r["L_m"], r["nu_factor"], r["exposure_a"], r["ice_fraction_of_channel_dissipation"], r["closed_form_1_minus_F_over_a"]))
print("E2")
for r in res["E2"]: print("  T=%.1f %-34s sill %.5f maxbed %.5f deliv %.2f Qend %.5f latent %.3e out %.3e res %.1e" % (r["T_inflow_C"], r["treatment"], r["sill_lowering_m"], r["max_bed_lowering_m"], r["delivered_m3"], r["outlet_Q_end_m3_s"], r["heat_to_ice_latent_J"], r["heat_out_sensible_J"], r["energy_res"]))
for k, v in res["E3"].items():
    print(k); [print("  ", {a: (round(b, 7) if isinstance(b, float) else b) for a, b in x.items()}) for x in v]
