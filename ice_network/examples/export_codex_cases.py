"""Export numerical results of the production machinery for the Codex limiting cases that the public interface can represent.

COMPATIBLE (run here): linear_reservoir, finite_drainage, symmetric_outlets, unequal_outlets [test-only prescribed edge laws on the production
storage/ODE/ledger code]; pre_overtopping [production Manning/GVF edge, real physics]; heat_to_melt [production segment kernel with the manifest's
cp, which differs from the production constant 4217 J/kg/K; no production constant is changed].
INCOMPATIBLE: thermal_front (needs a resolved advection solver, which the reduced model does not have; exported by the resolved reference solver
in ice_network_reference instead).
"""
import csv, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.join(HERE, "..", "..")
sys.path.insert(0, os.path.join(HERE, "..", "src")); sys.path.insert(0, os.path.join(HERE, "..", "tests")); sys.path.insert(0, ROOT)
from adapters import make_reservoir
from helpers import chan, const, lake
from ice_network.evolution.network import Network
from ice_network.thermal.channel_heat import segment_heat
from ice_network_verification.references import analytical_cases

OUT = sys.argv[1] if len(sys.argv) > 1 else os.path.join(HERE, "..", "records", "codex_candidate_exports")
os.makedirs(OUT, exist_ok=True)


def write(name, fields, rows):
    with open(os.path.join(OUT, name + ".csv"), "w", newline="") as fh:
        w = csv.writer(fh, lineterminator="\n"); w.writerow(["time_s"] + fields)
        for r in rows: w.writerow([f"{x:.17g}" for x in r])


cases = {c.name: c for c in analytical_cases()}


def run(net, v0, c, inflow=0.0):
    dt = c.duration_s / (c.samples - 1)
    return net.run([v0], [0.0], c.duration_s, dt, const(inflow), const(0.0), evolve=False, rtol=1e-12, atol=1e-14)


# linear reservoir and finite drainage
for name, law in (("linear_reservoir", lambda p: ("linear", p["conductance_m2_s"])), ("finite_drainage", lambda p: ("sqrt", p["coefficient_m5_2_s"]))):
    c = cases[name]; p = c.parameters; net = make_reservoir(p["area_m2"], [law(p)])
    r = run(net, p["area_m2"] * p["initial_head_m"], c)
    write(name, ["level_m", "storage_m3", "discharge_m3_s", "outflow_m3"],
          [(t, h["stage"][0], h["V"][0], h["Q"][0], h["C"][0]) for t, h in zip(r["t"], r["hist"])])
# parallel outlets
for name in ("symmetric_outlets", "unequal_outlets"):
    c = cases[name]; p = c.parameters; net = make_reservoir(p["area_m2"], [("linear", p["conductance_1_m2_s"]), ("linear", p["conductance_2_m2_s"])])
    r = run(net, p["area_m2"] * p["initial_head_m"], c)
    write(name, ["level_m", "storage_m3", "outlet_1_m3_s", "outlet_2_m3_s", "outlet_1_volume_m3", "outlet_2_volume_m3"],
          [(t, h["stage"][0], h["V"][0], h["Q"][0], h["Q"][1], h["C"][0], h["C"][1]) for t, h in zip(r["t"], r["hist"])])
# pre-overtopping through a real Manning/GVF edge (sill at 0.5 m)
c = cases["pre_overtopping"]; p = c.parameters
e = chan(0, 1, p["sill_height_m"], 0.5, 0.5, L=10.0, n=5); e["z"] = p["sill_height_m"] - 0.5 * np.linspace(0, 1, 5); e["zs"] = e["z"] + 1.0
net = Network([lake(0.0, p["area_m2"])], [dict(h_fixed=-5.0)], [e], melting=False)
dt = c.duration_s / (c.samples - 1)
r = net.run([p["area_m2"] * p["initial_head_m"]], [0.0], c.duration_s, dt, const(p["inflow_m3_s"]), const(0.0), evolve=False, rtol=1e-12, atol=1e-14)
write("pre_overtopping", ["level_m", "storage_m3", "inflow_volume_m3", "discharge_m3_s"],
      [(t, h["stage"][0], h["V"][0], h["Win"], abs(h["Q"][0])) for t, h in zip(r["t"], r["hist"])])
# heat_to_melt: production segment kernel, manifest constants
c = cases["heat_to_melt"]; p = c.parameters
C = p["water_density_kg_m3"] * p["water_cp_j_kg_k"] * p["discharge_m3_s"]
th_out, q = segment_heat(C, p["wall_conductance_w_k"], p["inlet_temperature_c"], 0.0)
rows = []
for i in range(c.samples):
    t = c.duration_s * i / (c.samples - 1); E = q * t; m = E / p["latent_heat_j_kg"]
    rows.append((t, th_out, E, m, m / p["ice_density_kg_m3"], m / p["water_density_kg_m3"]))
write("heat_to_melt", ["outlet_temperature_c", "heat_to_ice_j", "melt_mass_kg", "ice_volume_lost_m3", "meltwater_volume_m3"], rows)
print("exported to", OUT)
