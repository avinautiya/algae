"""Sensitivity of the model albedo to the sky condition, at the fitted Baseline-A parameters.

The S6 'albedo' is ALDN/ALUP with a remote cosine receptor (Cook 2020 archive, field logsheets: ALUP/ALDN pairs,
'Fore-optic RCR/8degree', natural light, sky from clear to thin/thick cloud) -> hemispherical albedo under the actual
sky. Baseline A models direct-beam illumination at noon SZA. This script quantifies, at the fitted parameters, how far
the same surface's albedo moves under fully diffuse illumination. It does not refit.
"""
import json, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, "..", "phase2"))
import biosnicar_bridge as BB
from forward import FixedOpticsForward
import operators as O

cfg = json.load(open(os.path.join(HERE, "config_baseline_A.json")))
fw = FixedOpticsForward(os.environ["BIOSNICAR_PATH"], cfg)
fd = FixedOpticsForward(os.environ["BIOSNICAR_PATH"], cfg)
fd.runner = BB.BioSNICARRunner(fd.root, incoming=3, direct=0)
m = (fw.wl_nm >= 400) & (fw.wl_nm <= 1300)
rows = []
for line in open(os.path.join(HERE, "..", "records", "inversion_m1", "baseline_A", "per_plot_fits.jsonl")):
    r = json.loads(line)
    if r["status"] != "ok": continue
    t = r["theta"]
    a_dir = fw.albedo(*t, r["sza"]); a_dif = fd.albedo(*t, r["sza"])
    d = (a_dif - a_dir)[m]
    rows.append(dict(sample=r["sample"], sza=r["sza"], rms_diffuse_minus_direct=float(np.sqrt(np.mean(d**2))), mean_diff=float(d.mean())))
rms = np.array([x["rms_diffuse_minus_direct"] for x in rows]); md = np.array([x["mean_diff"] for x in rows])
out = dict(n=len(rows), median_rms=float(np.median(rms)), p90_rms=float(np.percentile(rms, 90)), median_mean_diff=float(np.median(md)),
           note="diffuse minus direct at fitted parameters, 400-1300 nm; compare with the fit residual RMSE (median 0.0097)")
json.dump(dict(summary=out, rows=rows), open(os.path.join(HERE, "..", "records", "inversion_m1", "baseline_A", "sky_sensitivity.json"), "w"), indent=1)
print(json.dumps(out, indent=1))
