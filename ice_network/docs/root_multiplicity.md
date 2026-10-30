# Multiple roots of the entrance equation near the dry limit

Records: `records/root_study/` (scan, ramp, timing, failed first attempt), `../ice_network_eval/records/root_activation/` (S1 first-overtopping stress test). Scripts: `examples/root_multiplicity_study_v2.py` (sectioned, one process per timing variant, wall limit 240 s, evaluation cap 150 000), `../ice_network_eval/root_activation_s1.py`. The first version of the study (one process, one output file written at the end, rtol 1e-10) produced no output in 31+ minutes and was lost in the container restart of 14:38 UTC (see `FAILURE_LOG.md` F6); its timing section also exceeded a 600 s wall limit when re-run on the production variant (`FAILED_v1_rtol1e-10_timing_5e-2_floor_1e-3.json`). The v2 timing sections use rtol 1e-8 (documented change; no result of the v1 study is used anywhere).

## Governing equation
For one edge with stage `H_a` above the sill, tailwater `H_b`, the discharge `Q` is a root of `R(Q) = H_req(Q) - H_a = 0`, where `H_req(Q)` is the upstream head required to pass `Q` through the entrance (loss `k_e`), the Manning reach (standard-step) and the exit control. `Q = min(Q_outlet, Q_inlet)`.

## What was mapped (63 combinations: 3 geometries x 7 heads over the sill from 1e-9 to 0.1 m x substeps 4/16/64; 600-point log grid in Q plus bisection of every sign change)
* **37 of 63 combinations have no root at all** (the head is too small for any `Q > 0` to satisfy the Manning/critical-brink law): `Q = 0` is returned. For the 100 m, slope 2e-4 channel this holds up to head 0.01-0.03 m (the ramp test finds the first non-zero discharge at **0.0308 m over the sill**, identically for the production and the reference solver). Consequence: the activation time of a flat route is set by the thin-film Manning law, not by the sill elevation alone.
* **26 combinations have roots, in pairs**: one upward crossing (R from negative to positive with increasing Q: the stable, physical branch) and one downward crossing (the fold partner). Of 72 roots, 36 are upward and 36 downward. Twenty-three combinations have exactly one pair, three have 3-6 pairs.
* **Downward (unstable) roots** have `h_min <= 2.8e-4 m`; 25 of 36 have minimum depth below 10 micrometres, where a rough-wall Manning law is meaningless. **Upward roots**: 11 of 36 have `h_min < 10 micrometres` (all at heads 1e-7 and 1e-5 m), the other 25 have `h_min` from 1.9e-5 to 5.0e-2 m and Froude number 1 at the control.
* **Several pairs occur only at heads of 1e-5 to 1e-4 m in the 1 m drop geometry**, with discharges 1e-9 to 7e-9 m3/s and depths of 1-3 micrometres; their number depends on the substep count (4/16/64 substeps: 0/12/2 roots at head 1e-5) and is therefore a discretisation artefact of the micro-depth regime, not a property of the physical problem.
* Admissibility under the model's own regime (subcritical gradually varied flow with a critical control, depths large enough for a fully rough law): only upward roots with depth above the molecular-film scale qualify. No threshold is universal; the scan reports `h_min`, Froude number and control stations for each root so the choice is auditable.

## Branch selection rule (production solver)
Largest upward root: the bracket search starts at `Q_inlet` and halves until `R < 0`, then Brent. It never returns the unstable partner. The unchanged bisection reference solver can return a different root where several exist; on the ramp test (400 heads, flat 100 m channel, 16 substeps) the production and reference discharges differ by at most 1.3e-15 m3/s and the largest relative step between consecutive heads is 0.256 for both (the law is continuous in the head, with a steep switch-on near 0.031 m, identical in both solvers). The production solver is 9.7x cheaper per call than the reference over the level/tailwater sweep (`records/solver_sweep_stats.json`) and 4.3x faster in the activation timing test below (39 s against 168 s for the same 600 s problem).

## Regularisation (numerical, not physical)
`smooth_activation` multiplies nothing and moves nothing for `h_c >= h_floor`; below it, the discharge is blended to zero with `S(x) = x^2 (3 - 2x)`, `x = h_c/h_floor`. It is a **numerical regularisation** of the fold, not a physical selection rule: it removes a discontinuity in the micro-flow range where Manning friction does not apply. The raw law (`h_floor = 0`) is retained as a code path and in all records.

First overtopping of the frozen S1 geometry (reduced model, 9 stations, `dt_geom` 900 s, first 3 h), evaluation cap 60 000 law evaluations, wall limit 240 s:

| variant | status | wall (s) | law evaluations | first-activation time `t_Q1 > 1 %` of inflow (s) | delivered R1 (m3) |
|---|---|---|---|---|---|
| raw law (`h_floor` 0) | **stalled at the cap** in the macro-step starting at t = 1800 s | 59 | 60 001 | n/a | n/a |
| floor 1e-4 m | **stalled at the cap**, same macro-step | 64 | 60 001 | n/a | n/a |
| floor 1e-3 m (production default) | completed | 15 | 27 168 | 2584.13 | 87.846435 |
| floor 1e-2 m | completed | 5.6 | 18 560 | 2584.13 | 87.846435 |
| raw law, unchanged bisection solver | see `records/root_activation/raw_bisection.json` | | | | |

A floor of 0.1 mm is **not** sufficient to remove the stall; 1 mm and 10 mm give identical event time and delivered volume to all printed digits, so within the tested range the choice of floor above 1 mm does not influence the S1 endpoints.

## Effect on activation timing and volume in the flat 100 m geometry (activation of a slowly filling basin, rtol 1e-8, inflow 5e-2 and 5e-3 m3/s; `records/root_study/timing_*`)
* Activation times (first `Q` above 1e-9, 1e-6, 1e-4 m3/s) are the same for raw fast, floor 1e-4, floor 1e-3 and the bisection reference: 506.8 s at 5e-2 m3/s and 5068 s at 5e-3 m3/s; the stage reaches the sill at 200 s and 2000 s, i.e. activation lags the sill crossing by 307 s and 3068 s because of the 0.031 m threshold head above.
* Delivered volume differs between floors 1e-4/1e-3/raw by < 1e-10 m3 and from the bisection reference by 1e-10 m3 (relative 2e-11). Floor 1e-2 changes the (tiny) delivered volume during the first 600 s from 0.0484 to 0.0238 m3 (0.08 % of the 30 m3 inflow) and does not shift the activation times within the printed digits except the 1e-4 m3/s threshold (513.2 against 506.8 s).
* No stall occurred in any v2 timing variant (rtol 1e-8); the stall is specific to the S1 geometry (steep downstream slope, small crest) at first overtopping.

## Limits of this study
Three geometries, one friction coefficient (0.04), one entrance-loss coefficient; rising-head ramps only for the flat channel (falling-head hysteresis is not exercised because the law is single-valued in `H_a`; this follows from the quasi-steady formulation, not from a test). The adaptive sampling requested for transitions was implemented as bisection of every sign change on a fixed log grid of 600 points per case, not as refinement in the parameter space. An evaluation in which the activation time is compared with a resolved reference (not only with the other branch of the same law) is the S1-S5 evaluation; its activation times agree within 8-300 s depending on resolution (see `ice_network_eval/EVALUATION_REPORT.md`).
