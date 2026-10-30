"""Spatial, time-step and tolerance refinement of the coupled melting-channel case, and conservation diagnostics.

Case (idealised, not a forecast): one storage node (2000 m2 prismatic basin) feeds a 100 m open channel (0.5 m wide,
1 m bed drop, Manning 0.04) draining to a free outfall. External inflow 0.05 m3/s at 1.0 C; ice temperate. 6 h.
Reported per refinement: sill (upstream bed) lowering, delivered volume, heat to ice, geometry/energy ice-mass ratio,
and the three conservation residuals, normalised by the cumulative absolute fluxes.
"""
import json, os, sys, time
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "src")); sys.path.insert(0, os.path.join(HERE, "..", "tests"))
from helpers import chan, const, lake
from ice_network.evolution.network import Network

T_END = 6 * 3600.0


def case(n=9, dt=900.0, rtol=1e-9, msub=16):
    e = chan(0, 1, 100.5, 1.0, 0.5, L=100.0, n=n, msub=msub)
    net = Network([lake(99.0, 2000.0)], [dict(h_fixed=95.0)], [e])
    t0 = time.time()
    r = net.run([2000.0 * 1.55], [1.0], T_END, dt, const(0.05), const(1.0), rtol=rtol, atol=1e-10)
    b = r["budgets"]
    return dict(n=n, dt=dt, rtol=rtol, msub=msub, sill_lowering_m=100.5 - net.edges[0]["z"][0], width_inlet_gain_m=net.edges[0]["w"][0] - 0.5,
                delivered_m3=float(b["delivered_m3"][0]), heat_to_ice_J=float(b["heat_to_ice_J"]), ice_mass_kg=float(b["ice_mass_from_heat_kg"]),
                geometry_over_energy=float(b["geometry_energy_volume_ratio"]), water_res=float(b["water_residual_rel"]),
                mass_res=float(b["mass_residual_rel"]), energy_res=float(b["energy_residual_rel"]),
                final_stage=float(r["hist"][-1]["stage"][0]), final_storage_m3=float(r["hist"][-1]["V"][0]),
                outlet_Q_hourly=[float(h["Q"][0]) for t, h in zip(r["t"], r["hist"]) if abs(t / 3600.0 - round(t / 3600.0)) < 1e-9],
                storage_hourly_m3=[float(h["V"][0]) for t, h in zip(r["t"], r["hist"]) if abs(t / 3600.0 - round(t / 3600.0)) < 1e-9],
                seconds=round(time.time() - t0, 1))


def main(out):
    study = {"stations": [case(n=n) for n in (5, 9, 17, 33)],
             "dt_geom": [case(dt=d) for d in (3600.0, 1800.0, 900.0, 450.0, 225.0)],
             "msub": [case(msub=m) for m in (4, 8, 16, 32)],
             "rtol": [case(rtol=r) for r in (1e-6, 1e-8, 1e-10)]}
    # error against the finest member of each series, and observed order from consecutive differences (only for ordered refinements)
    analysis = {}
    for k, rows in study.items():
        ref = rows[-1] if k in ("dt_geom", "stations", "msub", "rtol") else rows[-1]
        if k == "dt_geom":
            ref = rows[-1]
        a = []
        for r in rows[:-1]:
            a.append(dict(n=r["n"], dt=r["dt"], msub=r["msub"], rtol=r["rtol"],
                          err_sill_lowering=abs(r["sill_lowering_m"] - ref["sill_lowering_m"]) / abs(ref["sill_lowering_m"]),
                          err_delivered=abs(r["delivered_m3"] - ref["delivered_m3"]) / ref["delivered_m3"],
                          err_final_storage=abs(r["final_storage_m3"] - ref["final_storage_m3"]) / ref["final_storage_m3"],
                          err_outlet_Q_6h=abs(r["outlet_Q_hourly"][-1] - ref["outlet_Q_hourly"][-1]) / abs(ref["outlet_Q_hourly"][-1]),
                          err_outlet_Q_max_hourly=float(np.max(np.abs(np.array(r["outlet_Q_hourly"]) - np.array(ref["outlet_Q_hourly"])) / np.abs(np.array(ref["outlet_Q_hourly"]))))))
        analysis[k] = a
    study["analysis_vs_finest"] = analysis
    json.dump(study, open(out, "w"), indent=1)
    for k, rows in study.items():
        if k == "analysis_vs_finest":
            continue
        print(k)
        for r in rows:
            print("  ", {a: (round(b, 8) if isinstance(b, float) else b) for a, b in r.items() if not a.endswith("hourly") and not a.endswith("hourly_m3")}, flush=True)


if __name__ == "__main__":
    main(sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "records", "convergence_study.json"))
