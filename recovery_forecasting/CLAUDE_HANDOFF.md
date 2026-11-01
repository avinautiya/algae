# Tonight's bounded integration handoff

This is a separate additive toolkit. Do not interrupt existing computations or
replace the active spatial solver. Import the package into its own folder and
run `python3 verify.py` first. The 45 tests should pass without OpenFOAM or new
data downloads. No large computation is launched by the test runner.

## Finish in this order

1. Complete the hydraulic partition/reference diagnosis already underway. Record
   what changed and which original geometry remains supported. Do not silently
   treat a graded-inlet geometry as verification of a different abrupt inlet.
2. Confirm exact clearing/forcing timestamps, physical-time water/thermal/latent
   ledgers, and event convergence. Keep failed and interrupted results visible.
3. Use `dynamic/` as an independent geometry-aware REDUCED baseline. Map its
   ideal control/storage geometry explicitly before comparing it with CFD.
   It adds sensible-heat storage and downstream meltwater accounting; it does
   not resolve the earlier entrance, porous or lateral-flow questions.
4. Freeze test inputs/predictions, operational contact threshold, persistence
   duration and post-flood assessment time. Include recovering, nonrecovering
   and never-disconnected cases. Do not tune the target after seeing outcomes.
5. Run the existing bounded spatial campaign only after its gate passes. Compare
   peak-only, storage-only, leakage-aware and geometry-aware predictors on the
   SAME target and common cases. First zero-depth contact is not sustained
   ordinary-flow recovery. A necessary energy exclusion is not a point forecast.
6. Use `audit/` for explicit abstentions and interval/right-censored timing.
   Its readiness output checks supplied metadata only; passing the checklist is
   not independent review of the underlying solver records.
7. Generate final figures from saved records: geometry and water routing,
   clearing/ordinary-return events, sustained toe delivery, local heat and melt,
   paired prediction errors, and numerical/closure sensitivity. Distinguish
   synthetic demonstration, development and held-out figures.

## Safe physical interpretation

For dynamic forecasts, demonstrate that fixed reservoir hypsometry and
quasi-steady outlet controls apply. Report the supplied meltwater/route-inflow
ratio bound; neglected feedback is an approximation. Uniform deposit resistance
cannot establish coarse-grain or porous-blockage novelty. If the outlet lowers
to the reservoir floor, the dynamic run refuses rather than adding a hidden sill.

The static predictor's corner values are sensitivity samples, not confidence
intervals or proven extrema. The spatial default uses a fractional-width
relaxation; use whole-cell path costs only for the declared discrete model.
All water/head/heat bounds require their own physical provenance.

## Deliverable boundary for tonight

Aim for a reproducible release containing the verified solver, frozen case
definitions, completed bounded comparison, saved predictions, figures and an
evidence-matched report. If the hydraulic gate still fails, publish a diagnostic
release with the forecast tools and unresolved discrepancy instead of calling
the data final. A submission-level novelty claim depends on the actual verified
result and its prior-art comparison, not the amount of code or computation.
