"""Separately controlled refinement of the NATIVE physical-time bare-ice solver (phase/swe2d_v5.py), with per-step (intermediate) ledger defects.
Three independent controls, each varied with the other two fixed:
  space      : dx 0.25, 0.125, 0.0625 m   (dt=0.25 s outer split step, cfl=0.4)
  outer step : dt 1, 0.5, 0.25, 0.125 s    (this API ties hydraulic window, source epoch and bed update to one step, so this is the combined split/morphology-update control)
  inner step : cfl 0.4, 0.2, 0.1 (hydraulic sub-stepping inside each window; dx=0.25, dt=0.25 s)
Case = the delivered idealised bare connected ice case (Q 0.01 m3/s, Tw 2 K, heff 5000 W/m2/K, 20 s). SHORT SUPPORTED CASE ONLY: this is not event-time (21 d) accuracy and not a contact-convergence result.
Reports max over ALL steps of mass/sensible/latent residuals (not only the final cancellation), fixed-ROI (all interior ice cells) geometric melt vs heat-implied melt, per-route source/flow, positivity water.
python3 refinement.py OUT.json"""
import json, sys, time, os
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(os.path.join(HERE, "..", "..", ".."))
sys.path.insert(0, os.path.join(ROOT, "research", "cliff_repair_2026_11_01", "phase"))
from swe2d_v5 import Sim
RHO_I, LF = 917.0, 3.34e5
def run(dx, dt, cfl, T=20.0):
    s = Sim(dx=dx, Lpool=.5, Lr=.5, Ltail=.25, Wl=.25, Q=.01, Tw=2., heff=5000., ramp_ratio=0, cfl=cfl)
    z0 = s.zb.copy(); act = np.zeros_like(s.ice); act[s.interior] = s.ice[s.interior] & ~s.wall[s.interior]
    n = round(T / dt); worst = dict(mass_kg=0.0, sensible_J=0.0, thermal_latent_J=0.0, phase_mass_kg=0.0, volume_source_m3=0.0, latent_minus_window_J=0.0)
    t0 = time.perf_counter(); pos = 0.0
    for _ in range(n):
        r = s.step_physical(dt); L = r["ledgers"]; src = r["source"]
        worst["mass_kg"] = max(worst["mass_kg"], abs(L["combined_ice_water_mass_residual_kg"])); worst["sensible_J"] = max(worst["sensible_J"], abs(L["executed_state_sensible_residual_J"]))
        worst["thermal_latent_J"] = max(worst["thermal_latent_J"], abs(L["physical_split_thermal_latent_residual_J"])); worst["phase_mass_kg"] = max(worst["phase_mass_kg"], abs(src["phase_mass_residual_kg"]))
        worst["volume_source_m3"] = max(worst["volume_source_m3"], abs(src["volume_source_residual_m3"])); worst["latent_minus_window_J"] = max(worst["latent_minus_window_J"], abs(src["latent_minus_window_heat_J"])); pos = L["positivity_added_water_m3"]
    wall = time.perf_counter() - t0; L = s.combined_ledgers(); geo = float(((z0 - s.zb) * act).sum() * s.dx ** 2); heat = L["morphological_latent_J"] / (RHO_I * LF)
    d = s.diagnostics(); toe, byp = float(s.last_phase_source["toe_ice_volume_m3"]), float(s.last_phase_source["bypass_ice_volume_m3"])
    return dict(dx=dx, dt=dt, cfl=cfl, steps=n, cells=int(s.nx * s.ny), melt_ice_ledger_m3=L["ice_volume_lost_m3"], melt_geometric_roi_m3=geo, melt_heat_implied_m3=heat, geometry_minus_heat_m3=geo - heat, Q_toe_m3_s=float(s.Qsec(s.toe_cols())),
                f_toe=d["f_toe"], water_volume_m3=s.liquid_volume(), positivity_added_water_m3=pos, final_mass_residual_kg=L["combined_ice_water_mass_residual_kg"], final_thermal_latent_J=L["physical_split_thermal_latent_residual_J"],
                max_over_steps=worst, last_step_toe_ice_m3=toe, last_step_bypass_ice_m3=byp, last_step_route_source=s.last_phase_source["route_source_diagnostics"], wall_s=wall)
def main(out):
    res = dict(scope="short supported idealised bare-ice case, native physical time; NOT event-time accuracy; engineering differences, not confidence intervals", controls={})
    res["controls"]["space"] = [run(dx, 0.25, 0.4) for dx in (0.25, 0.125, 0.0625)]
    res["controls"]["outer_step"] = [run(0.25, dt, 0.4) for dt in (1.0, 0.5, 0.25, 0.125)]
    res["controls"]["inner_cfl"] = [run(0.25, 0.25, c) for c in (0.4, 0.2, 0.1)]
    for name, rows in res["controls"].items():
        ref = rows[-1]
        for r in rows:
            for k in ("melt_ice_ledger_m3", "Q_toe_m3_s", "f_toe", "water_volume_m3"):
                r["abs_diff_" + k + "_to_finest"] = None if r[k] is None else r[k] - ref[k]; r["rel_diff_" + k + "_to_finest"] = None if (r[k] is None or ref[k] == 0) else r[k] / ref[k] - 1
    json.dump(res, open(out, "w"), indent=1, allow_nan=False, default=float)
    for name, rows in res["controls"].items():
        print(name)
        for r in rows: print("  dx %.4f dt %.3f cfl %.2f melt %.6e Qtoe %.6e rel(melt) %+.2e  max mass %.1e kg sens %.1e J tl %.1e J  pos %.1e  geo-heat %.1e wall %.2fs" % (r["dx"], r["dt"], r["cfl"], r["melt_ice_ledger_m3"], r["Q_toe_m3_s"], r["rel_diff_melt_ice_ledger_m3_to_finest"], r["max_over_steps"]["mass_kg"], r["max_over_steps"]["sensible_J"], r["max_over_steps"]["thermal_latent_J"], r["positivity_added_water_m3"], r["geometry_minus_heat_m3"], r["wall_s"]))
if __name__ == "__main__": main(sys.argv[1])
