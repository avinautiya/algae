"""Tests of the baseline predictor: restart equivalence, rejection of incompatible restarts, permitted horizon extension, budget reconstruction.  python -m ice_network_v2.test_predict (from the repo root)"""
import copy, json, os, shutil, subprocess, sys, tempfile
import numpy as np
ROOT = os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."); CASE = json.load(open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "cases", "baseline_synthetic.json")))
CASE["t_end_s"] = 3600; CASE["n_cells"] = 100; CASE["checkpoint_every_s"] = 900
def run(case, out, *args):
    p = os.path.join(out, "case.json"); os.makedirs(out, exist_ok=True); json.dump(case, open(p, "w"))
    r = subprocess.run([sys.executable, "-m", "ice_network_v2.predict", p, out, *args], cwd=ROOT, capture_output=True, text=True); return r
tmp = tempfile.mkdtemp(); res = {}
full = run(CASE, f"{tmp}/full"); assert full.returncode == 0, full.stderr
r1 = run(CASE, f"{tmp}/rs", "--stop-at", "1800"); r2 = run(CASE, f"{tmp}/rs", "--restart", f"{tmp}/rs/checkpoint.npz"); assert r2.returncode == 0, r2.stderr
a = np.load(f"{tmp}/full/snapshots.npz"); b = np.load(f"{tmp}/rs/snapshots.npz"); res["restart_bitwise_identical"] = bool(all(np.array_equal(a[k], b[k]) for k in a.files) and open(f"{tmp}/full/timeseries.csv").read() == open(f"{tmp}/rs/timeseries.csv").read())
# incompatible restarts must be refused
bad = {}
for name, mod in {"n_cells": lambda c: c.update(n_cells=120), "nu_factor": lambda c: c.update(nu_factor=2.0), "width": lambda c: c.update(width_m=0.6), "forcing_past": lambda c: c["inflow"].update(Q=[[0, 0.06], [86400, 0.06]]), "n_manning": lambda c: c.update(n_manning=0.05), "outlet": lambda c: c["outlet"].update(type="tailwater", stage_m=99.0)}.items():
    c = copy.deepcopy(CASE); mod(c); r = run(c, f"{tmp}/rs_{name}", "--restart", f"{tmp}/rs/checkpoint.npz"); bad[name] = (r.returncode != 0 and "REFUSED restart" in (r.stderr + r.stdout))
res["incompatible_restarts_refused"] = bad
# permitted: extend horizon (t_end 3600 -> 5400) from the 1800 s checkpoint of a run whose own t_end was 3600
c = copy.deepcopy(CASE); c["t_end_s"] = 5400; r3 = run(c, f"{tmp}/ext", "--restart", f"{tmp}/rs/checkpoint.npz")
# the restart checkpoint above was overwritten at the end of the run (t = 3600); extension continues from there
ful = run(c, f"{tmp}/full54"); ok = r3.returncode == 0 and ful.returncode == 0
if ok:
    a = json.load(open(f"{tmp}/full54/summary.json"))["outputs"]; b = json.load(open(f"{tmp}/ext/summary.json"))["outputs"]; res["extension_matches_uninterrupted_longer_run"] = all(abs(a[k] - b[k]) <= 1e-9 * max(abs(a[k]), 1) for k in a)
else: res["extension_matches_uninterrupted_longer_run"] = False; res["extension_error"] = (r3.stderr + ful.stderr)[-300:]
# budget reconstruction
s = json.load(open(f"{tmp}/full/summary.json")); B = s["budget"]
res["water_budget_closes"] = abs(B["initial_water_storage_m3"] + B["supplied_inflow_m3"] + B["meltwater_added_m3"] - B["delivered_m3"] - B["final_water_storage_m3"]) < 1e-8 * B["supplied_inflow_m3"]
res["energy_budget_closes"] = abs(B["initial_sensible_storage_J"] + B["supplied_sensible_J"] + B["dissipation_heating_J"] - B["exported_sensible_J"] - B["heat_to_ice_J"] - B["final_sensible_storage_J"]) < 1e-8 * B["supplied_sensible_J"]
res["latent_equals_heat_to_ice"] = abs(B["latent_energy_of_melted_ice_J"] - B["heat_to_ice_J"]) < 1e-9 * B["heat_to_ice_J"]
print(json.dumps(res, indent=1)); json.dump(res, open(os.path.join(os.path.dirname(os.path.abspath(__file__)), "test_predict_results.json"), "w"), indent=1)
shutil.rmtree(tmp); sys.exit(0 if res["restart_bitwise_identical"] and all(bad.values()) and res["extension_matches_uninterrupted_longer_run"] and res["water_budget_closes"] and res["energy_budget_closes"] and res["latent_equals_heat_to_ice"] else 1)
