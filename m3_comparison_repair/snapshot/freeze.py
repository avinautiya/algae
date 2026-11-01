"""Production freeze generator (M3). `python -m ice_cliff_survival.reactivation.freeze OUT_DIR` writes production_spec.json (+ SHA256SUMS). Run it, COMMIT the output, then launch the evaluation campaign from the committed spec (campaign.py refuses to start if the spec hash or git HEAD differs)."""
import hashlib, json, os, platform, subprocess, sys
import numpy as np, scipy, numba
HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
FILES = ["swe2d.py", "cases.py", "metrics.py", "run_case.py", "queue.py", "reduced.py", "verify_swe.py", "diagnose.py"]
EVAL = [  # name, layout, removal_day (None=none)
    ("E01_bare_none", "bare", None), ("E02_corridor_none", "corridor", None), ("E03_corridor_r12", "corridor", 12.0),
    ("E04_strip_none", "strip", None), ("E05_strip_r1", "strip", 1.0), ("E06_strip_r8", "strip", 8.0), ("E07_strip_r12", "strip", 12.0),
    ("E08_porous_none", "porous_strip", None), ("E09_porous_r1", "porous_strip", 1.0), ("E10_porous_r12", "porous_strip", 12.0)]
def sha(p): return hashlib.sha256(open(p, "rb").read()).hexdigest()
def git(*a): return subprocess.run(["git", "-C", ROOT] + list(a), capture_output=True, text=True).stdout.strip()
def predictions():
    from . import reduced as R, metrics as M
    out = {}
    for n, lay, rem in EVAL:
        if lay == "bare": ms = [R.fractional(None, 21.0, cover=0.0), R.fixed(None, 21.0, cover=0.0)]; ms[0]["method"] = "bare_reduced(evolving)"; ms[1]["method"] = "bare_fixed_geometry"
        elif lay == "porous_strip": ms = []        # reduced models have no seepage closure: no prediction made (stated, not invented)
        else: ms = R.all_methods("corridor" if lay == "corridor" else "strip", rem, 21.0)
        out[n] = [dict(method=m["method"], loss=M.events(np.array([x["day"] for x in m["series"]]), np.array([x["f_toe"] for x in m["series"]]), 1e-3, rem)["loss"], sustained_reconnection=M.events(np.array([x["day"] for x in m["series"]]), np.array([x["f_toe"] for x in m["series"]]), 1e-3, rem)["sustained_reconnection"], f0=m["series"][0]["f_toe"], f_end=m["series"][-1]["f_toe"], toe_melt_m3=m["toe_melt_m3"], series_days=[x["day"] for x in m["series"]][::4], series_f=[x["f_toe"] for x in m["series"]][::4]) for m in ms]
    return out
def main(out):
    os.makedirs(out, exist_ok=True); unc = json.load(open(os.path.join(HERE, "results", "numerical_uncertainty.json")))
    spec = dict(
        status="FROZEN BEFORE EVALUATION RUNS. Development configurations (dev*/conv* runs, removal at 1/4/8 d on strip at coarse resolution) were used to set geometry, removal schedule and metrics; evaluation configurations E01-E10 at the medium resolution are run only after this commit.",
        model_class="2-D depth-averaged free-surface hydrodynamics (shallow-water, hydrostatic) with coupled thermal incision (advected sensible heat; ice-floor lowering by the same exchange flux). No vertical flow, grain-scale underflow or 3-D turbulence; porous deposit = Ergun drag closure, not resolved seepage.",
        solver=dict(git_commit=git("rev-parse", "HEAD"), file_sha256={f: sha(os.path.join(HERE, f)) for f in FILES}, verification_results={f: sha(os.path.join(HERE, "results", f)) for f in sorted(os.listdir(os.path.join(HERE, "results"))) if f.endswith(".json") and f.startswith(("verify", "diagnose", "numerical", "partition"))},
                    toolchain=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__, numba=numba.__version__, platform=platform.platform()), imported_reduced_package_commit=open(os.path.join(ROOT, "imported", "ice_stream_reactivation_3e51620", "PINNED_COMMIT")).read().strip()),
        geometry=dict(dx=0.025, grid="uniform square cells; pool Lpool=3.0 m (entrance ramps run:rise 3 over the last 3 m, floor to -1.0 m at inlet), two lanes 0.5 m wide separated by a one-cell divider wall, reach 2.0 m, free-overfall tail (1.5 m drop, slope 0.5, 1.0 m)", z_toe_floor=-0.08, z_bypass_floor=0.0, depth_datum="bypass lane floor at t=0 is z=0",
                      deposit=dict(solid_volume_m3=0.1, footprint_m2=0.5, strip="1.0 m x 0.5 m (full toe-lane width), centred in the reach, trapezoidal ends (taper 0.1 m)", corridor="0.25 m wide at the toe-lane side wall, full reach length, trapezoidal ends and free edge (taper 0.1 m)", porous="strip footprint; Ergun drag phi=0.4, D=0.2 m, shielding k_debris=2 W/m/K", scaling="shape scaled so that all layouts have identical volume 0.1 m3 (peak heights 0.22 strip / 0.26 corridor at taper 0.1)")),
        physics=dict(Q_in=0.07, T_in=0.5, h_eff_W_m2K=200.0, k_debris=2.0, manning_n=0.02, rho_w=1000.0, cp=4217.0, rho_i=917.0, Lf=3.34e5, g=9.81, debris_temperature="T_m (export carries zero enthalpy relative to T_m)", initial="stage 0.15 m, zero velocity, uniform T=T_in, ice at T_m"),
        boundary=dict(inlet="prescribed discharge per width, T_in", outlet="zero-gradient on a supercritical overfall tail (outlet-length independence verified)", lateral="reflective walls", bed_exchange="h_eff*T on wet ice cells (shielded by debris resistance where covered)"),
        numerics=dict(scheme="first-order HLL with hydrostatic reconstruction, explicit Euler, CFL 0.4, quasi-steady hydraulics warm-started between bed updates (partition drift <2e-4 Q per 4 s window, storage <1e-3 Q)", bed_update_interval_h=6.0, horizon_days=21.0, hydraulic_cap_s=60.0),
        schedules=dict(removal_days={n: r for n, _, r in EVAL}, removal_basis="early=1.0 d (before contact loss at ~2.7 d in the strip development run); 8 d = geometric timescale (bypass incision rate x time ~ toe-bump height + critical head); 12 d = 1.5 x that. Exploratory status; not tuned after seeing evaluation output. Removal is prescribed, NOT a predicted flood result."),
        evaluation_runs=[dict(name=n, layout=l, removal_day=r, dx=0.025, dtm_h=6.0, days=21.0) for n, l, r in EVAL],
        reduced_predictions=predictions(),
        metrics=dict(contact="f_toe=Q_toe/Q_in at the reach mid-section >= theta; theta primary 1e-3, sensitivity (1e-4,1e-2,5e-2)", reconnection="f_toe>=theta continuously >=1 d after removal following a loss; brief rewetting = shorter episodes", continuous=["mean |f_pred-f_ref| over [0,21 d]", "delivered toe water fraction", "toe-reach cumulative melt", "controlling crest elevations"], local_not_global="toe melt is local lane-reach loss, not whole-cliff or glacier loss"),
        numerical_uncertainty=unc,
        acceptance=dict(budgets="max |water_resid| < 1e-10 and max |energy_resid_rel| < 1e-10 each update", timing="a method 'beats' another on event timing only if its error is smaller than the other's by more than the CFD numerical uncertainty on that event time", continuous="partition errors compared against U_f (numerical uncertainty), melt errors against relative uncertainty from convergence study"),
        failure_rules="run marked FAILED (kept in the table, not dropped) if exception, NaN, budget residual above acceptance, or wall time > 90 min; a failed run is repeated at most once after documenting the cause; the spec is then versioned (v2) and the affected runs identified.",
        runtime=dict(per_run_limit_min=90, workers=4, memory_note="~150-400 MB per run; queue admits a job only if MemAvailable > 1.5 GB"),
        analysis_plan="reduced-model forecasts (committed here) are compared with the medium-resolution CFD runs; no parameter of any reduced model is fitted to evaluation output; development-run information used only: removal times and the fact that contact loss occurs near 2.7 d in the strip layout.")
    p = os.path.join(out, "production_spec.json"); json.dump(spec, open(p, "w"), indent=1, default=float); h = sha(p); open(os.path.join(out, "production_spec.sha256"), "w").write(h + "\n"); print(h)
if __name__ == "__main__": main(sys.argv[1])
