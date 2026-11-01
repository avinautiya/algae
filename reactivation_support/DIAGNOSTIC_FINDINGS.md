# Independent acceleration and verification findings

Pinned source: avinautiya/algae, commit `c70c55319232d1dae3d3a7ee856327101df995ea`, `ice_cliff_survival/reactivation/`. This package never modifies Claude's active source or checkpoints. Claude may already have corrected an item after this pin.

## 1. Main partition lead: the abrupt inlet and reference assumptions

Fixed-bed, isothermal tests reproduce the reported discrepancy. At dx=0.05, changing Manning n from 0.02 to zero changes bare toe fraction from 0.678273 to 0.678885, and strip fraction from 0.249281 to 0.244823. This is much smaller than the approximately 0.07/0.11 discrepancies. Recomputing the actual numerical face flux rather than cell-centered hu changes fractions by about 0.0009 or less in these tests. Neither alone explains the problem.

Removing only the 1.5 m outlet drop, retaining a descending tail, also barely changes the initial fixed-bed partition. Thus that outlet change alone does not resolve the mismatch.

A matched **verification-only** experiment with pool length 3 m, no Manning friction and continuous descending tail gives:

| Layout / approach | dx | Numerical-face toe fraction | Ideal zero-loss control reference | Difference |
|---|---:|---:|---:|---:|
| Bare / abrupt approach | 0.05 | 0.679178 | 0.747107 | -0.067929 |
| Bare / graded approach | 0.05 | 0.728243 | 0.747107 | -0.018864 |
| Bare / graded approach | 0.025 | 0.735161 | 0.747107 | -0.011946 |
| Rectangular strip / abrupt approach | 0.05 | 0.244521 | 0.139771 | +0.104750 |
| Rectangular strip / graded approach | 0.05 | 0.151289 | 0.139771 | +0.011518 |
| Rectangular strip / graded approach | 0.025 | 0.135935 | 0.139771 | -0.003836 |

These results localize much of the mismatch to the entrance representation and its relation to the ideal reference. They do NOT prove a unique numerical root cause, prove a physical entrance-loss coefficient, or validate the actual abrupt-step experiment. The strip's signed error crosses zero under refinement; do not infer an asymptotic order from two points.

The graded geometry is a different physical problem. Do not silently substitute it into production to pass a tolerance. Use it to verify a matching reference, then resolve or bound the abrupt-step model difference. In a moving-bed run, a fixed approach ramp can become a new artificial sill after incision; do not merely set `tail_drop=0` or freeze the new ramp and assume the long-time boundary problem is solved.

The primary numerical literature documents limitations of hydrostatic reconstruction when local bed jumps exceed adjacent water depths: Delestre et al., DOI 10.1016/j.crma.2012.08.004, https://arxiv.org/pdf/1206.4986 . It motivates a local jump/depth diagnostic; it does not by itself establish the cause of this particular discrepancy.

Raw data: `records/partition_diagnostic.json`, `continuous_tail_diagnostic.json`, `smooth_control_diagnostic.json`, `smooth_control_refinement.json`.

## 2. Faster independent prismatic-channel reference

For flat wide-channel Manning flow integrated upstream from a critical exit, with unit discharge q and hc=(q²/g)^(1/3), direct integration yields

    n² q² L = (3/13)(h_up^(13/3)-hc^(13/3))
                -(3q²/(4g))(h_up^(4/3)-hc^(4/3)).

This comes from dx_up/dh = (1-Fr²)/Sf, Sf=n²q²/h^(10/3). The right side is strictly increasing for h>hc. Solve one monotone scalar root, then use H=h_up+q²/(2gh_up²). No ODE integration from a nearly singular critical depth is needed. A stable small-log-depth series avoids cancellation near hc.

`hydraulic_reference.py` matches an independently integrated GVF ODE in the unit tests. On the pinned two-lane reference, old/new toe fractions are 0.7465028617 and 0.7465028943. Local serial medians were 0.92466 s and 0.002220 s: approximately 416× for this reference routine, NOT for the whole CFD campaign. The formula assumes prismatic horizontal lanes, a critical exit, the wide-channel Manning approximation and no entrance loss. It is not an exact strip-layout oracle.

## 3. Phase mass source: real omission, separate from the partition issue

Pinned `Sim.lower_bed()` changes zb while leaving h unchanged and does not add a meltwater source in the hydraulic continuity equation. Consequently the within-window flux/storage ledger does not establish conservation across phase change.

A prescribed 100 W/m² melt flux over both ice lanes for one hour removes 0.0023508055 m³ of ice. It should generate 0.0021556886 m³ of water at rho_i/rho_w=0.917. The legacy geometry update generates zero. The isolated conservative operator reproduces the expected addition to about 1.3e-17 m³ and preserves sensible heat when latent heat is externally supplied.

This does NOT explain a 0.11 flow-partition discrepancy. For the nominal 0.07 m³/s, 0.5 K, 200 W/m²/K and 2 m² heated-area case, a 200 W total heat ceiling gives meltwater/input discharge <=8.6e-6. If neglect is intentional, declare and bound it as an approximation instead of claiming exact coupled phase-mass closure.

`phase_update.py` is an exchange operator with an explicit liquid-energy contract. In a closed adiabatic box use delta_E_liquid=-rho_i*Lf*dz. If surrounding boundary integration already supplies the latent heat, use zero; do not subtract latent heat a second time. Cold meltwater adds mass without sensible enthalpy at Tm and dilutes T=hT/h.

**Multirate caution:** `integration_adapter.py` only demonstrates this local update. It does not repair the full two-clock quasi-steady algorithm's cumulative physical-time ledger. Do not add both an integrated phase-volume jump and a continuous meltwater source without proving they are not double counted. Fast hydraulic spin-up fluxes are not automatically physical-time integrals. The caller must close that multirate ledger, or document a controlled leading-order approximation. Kinetic/potential energy and mixing dissipation are outside this sensible/latent exchange operator.

## 4. Event-time scheduling is currently quantized

Pinned `cases.run()` uses round(horizon/dtm), and removal occurs at the first regular update after its requested time. A real 0.5-day development fixture requested removal on day 0.3 and applied it on day 0.5: a 4.8-hour delay. The result metadata retains the requested time, while the export log holds the actual time.

`event_schedule.py` puts removal and the exact requested endpoint into the morphology grid. Record both requested and actual event times. `campaign_analysis.py` uses the actual export timestamp when appropriate for the constant-forcing prescribed-removal experiment. Flood experiments must explicitly provide the ordinary-flow return time.

## 5. Geometry realization and interpretation

`geometry_audit.py` shows that at dx=0.1 the nominal 50%-covered corridor realizes 0.4 m² covered area, while the strip realizes 0.5 m². At dx=0.05 and 0.025 both realize 0.5 m². Do not attribute their entire coarse-to-fine change to PDE error or estimate an order without accounting for this geometry change.

The one-cell divider represents different finite widths and pool widths on different grids. Report this approximation and distinguish it from an intentionally fixed physical wall thickness.

A porous deposit's depth×area is bulk volume. If phi=0.4 applies, its solid grain volume is 0.6×bulk volume. Do not label an equal-thickness solid/porous comparison as equal solid mass without the conversion. The first element of the pinned `porous` tuple is not an adjustable permeability: hydraulics replaces it with a Boolean activation value; Ergun drag actually depends on phi and D.

## 6. Analysis and connectivity helpers

The pinned contact detector does not gate its first upward crossing on removal time and does not distinguish a temporary rewet from sustained ordinary-flow contact. The new analysis reports complete active intervals, observation brackets, actual event times, persistence windows, horizon censoring, net delivered water and missing convergence metadata. It does not award production eligibility or invent uncertainty intervals.

Pinned `cases.crest()` minimizes the maximum elevation along straight columns. A winding route can have a lower geometric barrier. `connectivity.py` computes the true four-neighbor minimax barrier and returns its path. A tested synthetic grid has straight-column crest 1 but a winding path with crest 0. This is a geometry diagnostic, not a flow solution. Its algorithm is established graph theory, not a new theorem. It does not automatically explain the current rectangular-layout partition error.

The pinned runner also ignores returned hydraulic convergence flags and sets cumulative storage to zero. Preserve convergence/stop reasons per update and actual endpoint sensible storage; within-window residuals are not a substitute for a full physical-time ledger.
