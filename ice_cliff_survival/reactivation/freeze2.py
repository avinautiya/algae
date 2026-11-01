"""Production freeze v2 (M3 addendum): identical-footprint melt/heat comparison, corrected recovery classes, NEW removal-time bracket test, counterfactual without differential incision, finest-grid boundary check. v1 outputs are preserved untouched.
`python -m ice_cliff_survival.reactivation.freeze2 OUT_DIR` -> production_spec_v2.json (+sha256). campaign2.py refuses to run unless file hashes and spec hash match."""
import hashlib, json, os, platform, subprocess, sys
import numpy as np, scipy, numba
from . import reduced as R, metrics as M
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
FILES = ["swe2d.py", "cases.py", "metrics.py", "run_case.py", "queue.py", "reduced.py", "verify_swe.py", "diagnose.py", "campaign2.py"]
E = [("E01_bare_none", "bare", None), ("E02_corridor_none", "corridor", None), ("E03_corridor_r12", "corridor", 12.0), ("E04_strip_none", "strip", None), ("E05_strip_r1", "strip", 1.0), ("E06_strip_r8", "strip", 8.0), ("E07_strip_r12", "strip", 12.0), ("E08_porous_none", "porous_strip", None), ("E09_porous_r1", "porous_strip", 1.0), ("E10_porous_r12", "porous_strip", 12.0)]
BRACKET = (8.5, 9.0, 9.5, 10.0, 10.5)       # NEW removal days (never run before): reduced boundary t* = 9.50 d (f(t_r+1d)=1e-3), previous CFD data: 8 d recovers partially, 12 d not
def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()
def git(*a): return subprocess.run(["git", "-C", ROOT] + list(a), capture_output=True, text=True).stdout.strip()
def runs(tag="V2"):
    L = []
    for n, lay, rem in E: L.append(dict(name=tag + "C_" + n, role="v2_coarse_repeat_with_regions", layout=lay, dx=0.05, dtm_h=6.0, days=21.0, removal_day=rem, extra={}))
    L.append(dict(name=tag + "C_CF_corridor_none_uniform", role="counterfactual_coarse", layout="corridor", dx=0.05, dtm_h=6.0, days=21.0, removal_day=None, extra=dict(incision="uniform")))
    L.append(dict(name=tag + "C_CF_porous_none_uniform", role="counterfactual_coarse", layout="porous_strip", dx=0.05, dtm_h=6.0, days=21.0, removal_day=None, extra=dict(incision="uniform")))
    L.append(dict(name=tag + "M_strip_none", role="medium_regions_and_bed_snapshots", layout="strip", dx=0.025, dtm_h=6.0, days=15.0, removal_day=None, extra={}))
    for r in BRACKET: L.append(dict(name=f"{tag}M_strip_r{r:g}", role="NEW_bracket_test", layout="strip", dx=0.025, dtm_h=6.0, days=15.0, removal_day=r, extra={}))
    L.append(dict(name=tag + "M_CF_strip_none_uniform", role="counterfactual_medium", layout="strip", dx=0.025, dtm_h=6.0, days=15.0, removal_day=None, extra=dict(incision="uniform")))
    L.append(dict(name=tag + "M_CF_strip_r8_uniform", role="counterfactual_medium", layout="strip", dx=0.025, dtm_h=6.0, days=15.0, removal_day=8.0, extra=dict(incision="uniform")))
    return L
def predictions():
    ser = lambda a: ([x["day"] for x in a["series"]], [x["f_toe"] for x in a["series"]]); out = dict(reduced_boundary_day=9.5046, bracket={}, footprint={}, counterfactual={})
    ref = R.fractional(None, 21.0, cover=0.0)["series"][0]["f_toe"]
    for r in BRACKET:
        a = R.topology_lanes("strip", r, days=15.0); d, f = ser(a); c = M.classify(np.array(d), np.array(f), r, f_ref=ref); out["bracket"][f"{r:g}"] = dict(f_at_removal_plus_1d=float(np.interp(r + 1, d, f)), cls=c["cls"], R=c["R"], f_end=f[-1])
    for n, lay, rem in E:
        if lay == "bare": continue
        if lay == "porous_strip": out["footprint"][n] = None; continue
        a = R.topology_lanes("corridor" if lay == "corridor" else "strip", rem, days=21.0); out["footprint"][n] = dict(melt_m3=a["melt_m3"], area_m2=a["area_m2"], mean_heat_W={k: a["melt_m3"][k] * 917 * 3.34e5 / (21 * 86400) for k in a["melt_m3"]})
    for lay, rem, days in (("strip", None, 15.0), ("strip", 8.0, 15.0), ("corridor", None, 21.0)):
        a = R.topology_lanes("corridor" if lay == "corridor" else "strip", rem, days=days, uniform_incision=True); d, f = ser(a); out["counterfactual"][f"{lay}_{rem}"] = dict(f0=f[0], f_end=f[-1], min_f=min(f), contact_lost=bool(min(f) < 1e-3))
    return out
def main(out, tag="v2"):
    os.makedirs(out, exist_ok=True)
    spec = dict(status="FROZEN " + tag + " BEFORE " + tag + " RUNS (v1/v2 outputs preserved). v3 = the v2 specification re-run with solver v3 (tail-clip correction: v1/v2 lifted the overfall tail floor to pool_floor after the first bed update, leaving hydraulics unconverged, storage up to 3% of Q; effect on coarse partition <=2e-4 in f). Predictions and criteria unchanged from v2. Development basis: v1 evaluation results are now DEVELOPMENT information for this addendum; the bracket days, counterfactual and finest-grid check are new.",
        model_class="2-D depth-averaged free-surface hydrodynamics (hydrostatic shallow water) with coupled thermal incision; idealised obstacles (solid raised bed / Ergun drag); not natural boulder behaviour.",
        solver=dict(git_commit=git("rev-parse", "HEAD"), file_sha256={f: sha(os.path.join(HERE, f)) for f in FILES}, v1_spec_sha256=open(os.path.join(HERE, "frozen", "production_spec.sha256")).read().strip(), toolchain=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__, numba=numba.__version__),
                    code_change_v2="ADDITIVE: region heat/melt accounting, bed snapshots, optional incision='uniform' counterfactual. Default physics unchanged (verified: v2 coarse strip f_toe series bit-identical to v1 for the first 6 updates; full 21-d identity is checked in analysis and v1 results kept)."),
        runs=runs(tag.upper()), geometry_and_physics="identical to v1 spec (frozen/production_spec.json) except where 'extra' says otherwise",
        counterfactual=dict(name="differential incision disabled", definition="the same total melt volume of each bed update is spread uniformly over ALL ice cells (reach + entrance ramps) so routes/shielded and open cells lower equally; heat extraction, latent ledger and hydraulics are unchanged. Prediction: no stage-lowering relative to the strip crest => no contact loss and full recovery if removal."),
        identical_footprint_comparison=dict(regions="toe_fp = deposit footprint cells of the toe reach (0.5 m2 for strip and corridor), toe_open = other toe-reach cells (corridor: 0.5 m2 = reduced 'bare corridor lane'; strip: not represented in the reduced model), byp = bypass reach (1.0 m2); ramps excluded (ice, in neither reduced model)", compare="cumulative melt volume and mean heat delivery per region, CFD vs reduced topology_lanes() (identical to imported topology.run), at 21 d", note="reduced and CFD use the SAME insulation closure (series resistance), so melt agreement tests water delivery to the footprint, not the closure"),
        recovery_classes=dict(definition=M.classify.__doc__, thresholds=dict(full=M.RECOVERY_FULL, partial=M.RECOVERY_PARTIAL), f_ref="bare-route share from the bare control run at the same resolution (CFD E01: 0.735-0.738; reduced 0.749)", applied_to_v1="E03 maintained; E05 maintained(+decline); E06 lost-recovered-partial; E07/E10 lost-not-recovered; E09 maintained(+decline)"),
        bracket_test=dict(removal_days=list(BRACKET), horizon_days=15.0, grid="medium dx=0.025", pass_criteria="CFD boundary (smallest removal day with f(t_r+1d)<1e-3, linear interpolation over the 5 runs and v1 points 8 d and 12 d) within +-1.0 d of reduced t*=9.50 d; class sequence monotone in removal day"),
        finest_grid_boundary_check=dict(method="fixed-bed post-removal steady flows: bed snapshots of V2M_strip_none (medium) at removal days 8..12 step 0.25, deposit removed, resolved at dx = 0.05, 0.025, 0.0125 (bilinear transfer of the bed); records f_toe at removal+0 per grid; boundary = day where f crosses 1e-3", pass_criteria="boundary shift between dx=0.025 and 0.0125 < 0.5 d; if not, report the finest-grid boundary and flag unresolved. LIMITATION: the bed history is evolved at dx<=0.025 only (bed history spatial sensitivity <1% in melt and 0.25 d in loss time)"),
        predictions=predictions(),
        acceptance="water/energy residual < 1e-10 per update; failed runs reported, not dropped", failure_rules="as v1; wall limit 90 min/run", runtime=dict(workers=4, expected_wall_min="~75", expected_medium_run_min="~25-35 under 4 concurrent jobs"))
    p = os.path.join(out, f"production_spec_{tag}.json"); json.dump(spec, open(p, "w"), indent=1, default=float); h = sha(p); open(os.path.join(out, f"production_spec_{tag}.sha256"), "w").write(h + "\n"); print(h); print(json.dumps(spec["predictions"], indent=1, default=float)[:2500])
if __name__ == "__main__": main(sys.argv[1], sys.argv[2] if len(sys.argv) > 2 else "v2")
