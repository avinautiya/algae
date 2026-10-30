"""Conservation (water, mass in kg, energy) on nonzero-throughput cases, kept separate from solution accuracy.
Cases isolate the dissipative-heating switch: heat_dissipation=False removes Phi from the thermal problem AND from the energy budget
without touching the friction law (head losses hf, hence discharge and storage, are unchanged)."""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src")); sys.path.insert(0, os.path.join(HERE, "..", "tests"))
from helpers import chan, const, lake
from ice_network.constants import RHO_I, RHO_W
from ice_network.evolution.network import Network


def case(name, T0, melting, evolve, heat_dissipation):
    e = chan(0, 1, 100.5, 1.0, 0.5, L=100.0, n=9)
    net = Network([lake(99.0, 2000.0)], [dict(h_fixed=95.0)], [e], melting=melting, heat_dissipation=heat_dissipation)
    r = net.run([2000.0 * 1.55], [T0], 6 * 3600.0, 900.0, const(0.05), const(T0), rtol=1e-10, atol=1e-11, evolve=evolve)
    b = r["budgets"]; Q = [h["Q"][0] for h in r["hist"]]
    vol_ice = b["ice_volume_geometry"]; vol_liq = b["melt_water_m3"]
    out = dict(case=name, T_inflow_C=T0, melting=melting, heat_dissipation=heat_dissipation,
               throughput=dict(input_m3=b["input_m3"], delivered_m3=float(b["delivered_m3"][0]), delivered_over_input=float(b["delivered_m3"][0] / b["input_m3"]), Q_min=float(min(Q)), Q_max=float(max(Q))),
               water=dict(residual_m3=b["water_residual_m3"], scale_m3=b["water_scale_m3"], relative=b["water_residual_rel"]),
               mass_kg=dict(in_kg=b["mass_in_kg"], out_kg=b["mass_out_kg"], stored_change_kg=b["mass_stored_change_kg"], residual_kg=b["mass_residual_kg"], scale_kg=b["mass_scale_kg"], relative=b["mass_residual_rel"]),
               ice_to_water=dict(ice_mass_lost_from_heat_kg=b["ice_mass_from_heat_kg"], ice_mass_lost_from_geometry_kg=b["ice_mass_from_geometry_kg"],
                                 ice_volume_lost_geometry_m3=vol_ice, ice_density=RHO_I, liquid_water_mass_gained_kg=b["ice_mass_from_meltwater_kg"],
                                 liquid_water_volume_gained_m3=vol_liq, water_density=RHO_W,
                                 liquid_volume_over_ice_volume=(vol_liq / vol_ice) if vol_ice else None, expected_ratio_rho_i_over_rho_w=RHO_I / RHO_W),
               energy_J=dict(sensible_in=b["sensible_in_J"], sensible_out_outlets=b["sensible_out_J"], sensible_storage_change=b["sensible_storage_change_J"],
                             latent_to_ice=b["heat_to_ice_J"], dissipation_total=b["dissipation_J"], dissipation_in_channel=b["dissipation_channel_J"], dissipation_at_exit_node=b["dissipation_exit_J"],
                             residual=b["energy_residual_J"], scale=b["energy_scale_J"], relative=b["energy_residual_rel"]),
               geometry_over_energy_ice_volume=b["geometry_energy_volume_ratio"] if b["ice_volume_energy"] else None,
               sill_lowering_m=100.5 - net.edges[0]["z"][0])
    return out


cases = [case("A melting, dissipation counted", 1.0, True, True, True),
         case("B melting, dissipation excluded", 1.0, True, True, False),
         case("C dissipation-only melting (inflow at 0 C)", 0.0, True, True, True),
         case("D no-melting control, dissipation counted", 1.0, False, False, True),
         case("E no-melting control, dissipation excluded", 1.0, False, False, False)]
json.dump(cases, open(os.path.join(HERE, "..", "records", "budget_tables.json"), "w"), indent=1, default=float)
for c in cases:
    print(c["case"]); print("  throughput", {k: round(v, 5) for k, v in c["throughput"].items()})
    print("  rel residuals water %.1e mass %.1e energy %.1e" % (c["water"]["relative"], c["mass_kg"]["relative"], c["energy_J"]["relative"]))
    print("  ice->water", {k: (round(v, 6) if isinstance(v, float) else v) for k, v in c["ice_to_water"].items()})
    print("  energy", {k: float("%.4g" % v) for k, v in c["energy_J"].items()}, "sill", round(c["sill_lowering_m"], 6), "geo/energy", c["geometry_over_energy_ice_volume"])
