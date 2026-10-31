# Independent contribution for later review

This directory is added only on `codex/unsteady-seal-model`. No Claude worktree, branch, job, checkpoint, physical parameter or previous frozen result was changed. Do not interrupt an active job to integrate it.

## Concrete result to inspect

The candidate model retains fast seal-shape evolution and eliminates the spatial hydraulic solve into a derived state/renewal law. MATCHED_LINEAR_CLOSURE.md removes the outer translation mode; MATCHING_THEOREM.md proves the canonical hydraulic pond limit. This is more specific than adding a generic delay or a neural network.

Review the proof before scheduling computation. Check especially the two different seal heights, canonical scaling, the downstream normal impedance, far-field intercept evolution and the distinction between a positive memory kernel and a new instability. STABILITY_AND_SLOW_LIMIT.md intentionally recovers the known canonical threshold. The mathematical work does not prove historical originality.

## Implementation to reproduce

Run the existing test suite and manufactured CLI cases in README.md. The fixed-frame finite interval and outer-matched linear engine are different boundary formulations. Their head and ledger labels must remain distinct. A small total-head bookkeeping residual is not proof of a physical water/enthalpy budget. A nonlinear hydraulic failure stops the prototype; do not clip velocities, substitute a frozen rating or silently change branch.

The serialized diagnostics include finite-cutoff normal mismatch, translation-identity defects, finite hydraulic replay differences, criticality margins and the stage/kinetic correction. Preserve them in any comparison. The hydraulic sensitivity required spatial refinement; a 129-station result alone did not meet the independent derivative tolerance.

## Scientific claim boundary

This supplies a candidate resolution of a **canonical linear part** of the explicitly identified fast-inner modeling issue, plus a finite nonlinear implementation. It does not supply the complete nonlinear outer kinetic relation. Large branch changes, dry/critical continuation, environmental forcing ranges and independent CFD comparison remain unestablished.

Do not rename existing negative heat-prediction results as successes for this model. Do not claim the first dynamic rating, first ice instability, first moving seal, first phase-change mass source, or first coupled hydraulic/ice solver. The existing dynamic and coupled models listed in PRIOR_ART.md are necessary comparators. The meaningful eventual contribution would be a derived nonlinear selection law and consequential predictive correction that those strong comparators miss.

No production campaign, training run, parameter scan or frozen acceptance rule is launched or authorized by this handoff alone. The present records are verification, not a forecast-accuracy study.
