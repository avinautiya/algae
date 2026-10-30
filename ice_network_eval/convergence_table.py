"""Publishes CONVERGENCE_TABLE.md from records/convergence_summary.json (generated from per-run records)."""
import json, os
HERE = os.path.dirname(os.path.abspath(__file__))
K = json.load(open(os.path.join(HERE, "records", "convergence_summary.json")))
def f(x, n=3):
    return "-" if x is None else (x if isinstance(x, str) else f"{x:.{n}g}")
L = []
P = L.append
P("""# Convergence study (frozen scenarios S1-S4; records `records/runs/*conv*`, summary `records/convergence_summary.json`)

Idealised controlled problems, not a forecast. All values are generated from immutable per-run records by `convergence_analysis.py` / `convergence_table.py`.

## Definitions
* **R1, R2**: the two routes out of the lake. R1: 100 m long, 0.4 m wide, bed slope 0.005, lower sill. R2: 60 m long, 0.4 m wide, bed slope 0.02, sill above R1 by 0.05 m (S1, S4, S5), 0.02 m (S2) or 0.10 m (S3).
* **Incision** = lowering of the route's sill (bed elevation at the lake end of the route) after 8 h, in mm. **Incision error** = reduced minus reference; **relative incision error** = that difference divided by the reference sill lowering of the same route at the same time (the 800-cell reference in this table; the 400-cell reference in `FROZEN_S1_S5_ANALYSIS.md`).
* **Hydrograph RMSE/I**: root-mean-square difference between reduced and reference outflow of both routes on the reference's 60 s grid, combined as sqrt(R1^2 + R2^2), divided by the (mean) inflow.
* **Delivered-volume error**: difference of the 8 h volume through each route divided by the reference volume of that route. **Route (outlet choice)**: the route with the larger delivered volume. **t_act2**: first time Q_R2 > 1 % of the inflow. **t_dom2 (switching time)**: first time Q_R2 > Q_R1; "-" when it does not occur in either run (no switch).
* Baseline configuration: 17 stations per route, geometry step 450 s, 16 hydraulic substeps per segment, ODE rtol 1e-7 / atol 1e-9. One factor is varied at a time: stations 5/9/17/33/65 (spacing 25/12.5/6.25/3.1/1.6 m for R1, 15/7.5/3.75/1.9/0.94 m for R2); geometry step 112.5-1800 s; substeps 2-32; rtol 1e-9 to 1e-5. Reference: 100/200/400/800 cells at CFL 0.4, and 200 cells at CFL 0.2 and 0.1 (time-step test); target = 800 cells.

## How the reference uncertainty was estimated
The reference is a first-order upwind scheme (observed order 0.5-1.06 in the 100/200/400 frozen runs). The uncertainty of the 800-cell target is estimated as the change between the 400- and 800-cell references, |Q_800 - Q_400|: for an order-1 scheme this equals the 800-cell discretisation error; for an order-0.5 scheme the error would be 2.4 times larger, so the values below are a lower-range estimate (both are given for the hydrograph). It is an estimate from two meshes, not a bound; the 200-cell CFL test shows the reference is mesh-limited, not time-step-limited (hydrograph error vs 800 cells 0.0129 / 0.0125 / 0.0123 at CFL 0.4 / 0.2 / 0.1 in S1).
""")
for scn, v in K.items():
    u = v["reference_uncertainty_400_to_800"]; t = v["target"]
    P(f"## {scn}\n")
    P(f"Target (800 cells): route {t['route']}; t_act2 {f(t['t_act2_s'],5)} s; t_dom2 {f(t['t_dom2_s'],5)} s; sill lowering R1/R2 {t['sill_lowering_m'][0]*1e3:.2f} / {t['sill_lowering_m'][1]*1e3:.2f} mm; delivered R1/R2 {t['delivered_m3'][0]:.0f} / {t['delivered_m3'][1]:.0f} m3.")
    P(f"**Reference uncertainty (|400 - 800|)**: hydrograph RMSE/I {f(u['hydrograph_rmse_over_I'],2)} (x2.4 if order 0.5: {f(2.4*u['hydrograph_rmse_over_I'],2)}); t_act2 {f(u['t_act2_s'],2)} s; t_dom2 {f(u['t_dom2_s'],3)} s; sill lowering R1/R2 {u['sill_lowering1_m']*1e3:.2f} / {u['sill_lowering2_m']*1e3:.2f} mm.\n")
    P("| stations | route ok | hydrograph RMSE/I | order vs prev. | t_act2 err (s) | t_dom2 err (s) | delivered err R1/R2 (%) | incision err R1/R2 (mm) | rel. incision err R1/R2 (%) | runtime* (s) | law evals | s/eval (ms) |"); P("|" + "---|" * 12)
    for x in v["sweeps"]["stations"]:
        P(f"| {int(x['value'])} | {'yes' if x['route_match'] else 'NO'} | {f(x['hydrograph_rmse_over_I'],3)} | {f(x.get('observed_order_vs_previous'),2)} | {f(x['t_act2_err_s'],4)} | {f(x['t_dom2_err_s'],5)} | {x['delivered_rel_err'][0]*100:+.1f} / {x['delivered_rel_err'][1]*100:+.1f} | {x['sill_err_m'][0]*1e3:+.2f} / {x['sill_err_m'][1]*1e3:+.2f} | {x['sill_rel_err'][0]*100:+.1f} / {x['sill_rel_err'][1]*100:+.1f} | {f(x['solver_s'],3)} | {x['law_evals']} | {1e3*x['solver_s']/x['law_evals']:.2f} |")
    P("\n| factor | value | route ok | hydrograph RMSE/I | t_act2 err (s) | incision err R1/R2 (mm) | rel. R2 incision err (%) | runtime* (s) | law evals |"); P("|" + "---|" * 9)
    for fac in ("dt_geom", "msub", "rtol"):
        for x in v["sweeps"][fac]:
            P(f"| {fac} | {x['value']:g} | {'yes' if x['route_match'] else 'NO'} | {f(x['hydrograph_rmse_over_I'],3)} | {f(x['t_act2_err_s'],4)} | {x['sill_err_m'][0]*1e3:+.2f} / {x['sill_err_m'][1]*1e3:+.2f} | {x['sill_rel_err'][1]*100:+.1f} | {f(x['solver_s'],3)} | {x['law_evals']} |")
    P("\n| reference mesh / CFL | hydrograph RMSE/I vs 800 | runtime* (s) | steps | t_act2 (s) | R2 sill lowering (mm) |"); P("|---|---|---|---|---|---|")
    for k, r in v["reference"].items():
        P(f"| {k.replace('_cfl', ' cells, CFL ')} | {f(r['hydrograph_rmse_over_I_vs800'],3)} | {f(r['solver_s'],3)} | {r['nsteps']} | {f(r['t_act2_s'],5)} | {r['sill_lowering2_m']*1e3:.2f} |")
    P("")
P("*Runtimes were measured with two campaign workers running concurrently on a four-core machine (and, for part of the study, one further job); they are indicative, not matched-resource timings. The serial timing set is reported separately.\n")
P("""## Which errors refine and which persist (numbers from the tables above)
* **Outlet choice**: wrong only at 5 and 9 stations in S1 (a near-tie: reference delivered 546 against 490 m3, 5.6 % margin); correct at every setting in S2-S4 and from 17 stations in S1. Not affected by geometry step, substeps or rtol.
* **Switching time (S1, S2 only)**: error shrinks with stations in S1 from -10 177 s (5) to -726 s (65), against a reference uncertainty of 680 s: resolved to the reference uncertainty only at 65 stations. S2: -2404 s to -96 s (reference uncertainty 39 s: not yet resolved at 65). Always early (reduced model overtakes sooner).
* **Activation time of R2**: -226 to -310 s at 5-17 stations, +14 to -4 s at 33-65 stations; reference uncertainty 1-2 s. Refines, and flips sign between 17 and 33 stations.
* **Discharge (hydrograph)**: 0.054 -> 0.0050 (S1) with stations; second order between 17 and 33 stations, then flat (0.0058 -> 0.0050) at roughly the size of the reference uncertainty (0.0044; 0.011 if the reference converges at order 0.5). Independent of geometry step, substeps and rtol (changes below 0.001 absolute between the extremes).
* **Delivered volume**: R2 over-delivered by +12.4 % (5 stations) -> +0.8 % (65) in S1, +9-12 % -> +0.7-0.9 % in S2-S4; R1 within -5.6 % to +1.2 %. Settles to +0.7-1.1 % (R1 and R2), above any estimated reference uncertainty of volume but small.
* **Incision**: the R2 over-prediction **decreases monotonically with station refinement**, so it is a resolution error and not a structural offset in the tested range: S1 21.7 -> 16.7 -> 10.2 -> 3.8 % (5, 17, 33, 65 stations), S2 20.0 -> 3.6 %, S3 18.0 -> 3.5 %, S4 26.0 -> 4.5 % (error in mm: 11.97 -> 2.12 for S1). It does **not** depend on the geometry step (16.5-16.8 % over 112-1800 s), substeps (16.2-16.7 %) or rtol (S1), so the first-order geometry update is not limiting and no higher-order coupling was implemented. Residual at 65 stations: +2.1 mm (R2) against a reference uncertainty of 0.28 mm, i.e. still 7x the estimated reference uncertainty; whether it reaches zero cannot be told from five stations values (observed order 33 -> 65: 1.4).
* **A possibly structural remainder is visible in R1 incision**: +0.83, +1.11, +0.86 mm (S1; 1.8-2.4 %) at 17, 33, 65 stations does not shrink, against a reference uncertainty of 0.11 mm. Candidate causes (not tested): the sill value in the reduced model is a station value averaged from adjacent segments while the reference uses a cell value; different wall-retreat idealisation (D' = D + dz) from the reference's per-cell melt; unequal thermal-entrance treatment. This is a ~2 % difference in R1 and is reported as unexplained.
* **Cost**: not monotone in stations. Time per discharge-law evaluation rises from 0.56 ms (5 stations) to 1.8 ms (17) to 8.2 ms (33) and 9.0 ms (65), while the number of evaluations varies 49 704 / 14 840 / 5 852 / 8 716 / 5 392; the 33-station run is the most expensive (71 s) and the 65-station run costs less (48 s). The profile (`records/profile/`) attributes 90-97 % of the time to the discharge-law solve, with mean Brent iterations rising from 4.9 (17 stations) to 7.6 (33) and max 18 -> 36 (bracket search at near-dry states). No optimisation was implemented.
""")
open(os.path.join(HERE, "CONVERGENCE_TABLE.md"), "w").write("\n".join(L)); print(len(L))
