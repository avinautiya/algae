# Independent thermal integration review — pre-fix findings

Scope: read-only review of the parent integration at the exact source hashes recorded in thermal_review_pre_fix.json. Twenty seeded warm synthetic configurations were written before execution. No model parameters were changed in response to outcomes. These are implementation/applicability checks, not field accuracy or independent scientific validation.

## Observed status

- 14 configurations completed, including 9 normal-slab and 5 lumped cases. Every completed case had at least one wet-debris step.
- 2 configurations correctly refused incompatible routing/contact head requirements.
- 4 configurations unexpectedly refused numerical representation/class allocation states.
- Maximum completed normalized global energy residual: 4.42e-16 of gross incoming/atmospheric energy. Maximum geometric phase-volume residual: 5.51e-15m3. Maximum water bookkeeping residual: 1.20e-9m3 (this bookkeeping is not independent mass verification).
- Cold-reservoir wet-debris input is correctly refused. Zero-inflow atmospheric melt is refused without a drainage contract and accepted with the declared synthetic free-outlet contract.

## Concrete bugs to fix before production freeze

1. T11/T12: wet+dry class areas round one ULP above the exact total supplied by geometry, causing a false material-contraction refusal. The discrepancy is−1.3877787807814457e-17m2 in both reproductions. Do not silently allow genuine contraction. Preserve the exact requested total as a canonical class-state reference while retaining per-cell enthalpy residuals.
2. T16: the geometry adapter correctly has zero bare-water area, but the parent algebraic subtraction yields 1.545e- 14J of bare-water energy. Derive conductance areas from the same authoritative class partition used for geometric allocation; do not discard the resulting phantom energy at the end.
3. T19: the geometry has zero dry-covered area, but covered_total minus wet_covered yields 1.137e-13J of dry-covered heat. The temperate phase energy round-trip energy/(rhoL)*rhoL can change one ULP. Preserve the input latent energy exactly when no fresh-ice warming is charged, and use the independently computed dry-column interface heat instead of subtracting two nearly equal quantities.

thermal_failure_details_pre_fix.json contains exact values and stack locations. The parent has been notified and is applying fixes. The reviewer did not edit parent files.

## Coupling accounting inspected

The upstream floor removes water heat once; the bare/covered face then shares its remaining water stream. The bypass receives its own routed discharge. Covered-interface energy is debited through phase change once. Separate wet/dry normal profiles are retained rather than overwritten with a global mean; the mean scalar is a bookkeeping projection. Area growth imports declared incoming material enthalpy. No hidden wet/dry heat exchange occurs at a fixed classification. These observations do not verify the imposed hydraulic or environmental closures.

`prescribed_all_cover` consistently supplies the declared external surface flux on both wet and dry covered areas. It requires a prepared external heat input excluding the separately modelled water-convective term; calling an inundated air-only flux environmental truth would be an input-model mistake. `fixed_cover_interval_m` avoids the physical material-transport ambiguity of a prescribed cover fraction following a lowering floor, but remains an Eulerian mask, not a mobile debris model.

## Output-state issue

The exported wet_debris_thermal.partition is post-heat/pre-geometry-remap, while exported scalar reservoirs and geometry are post-remap. Those areas/energies need not coincide. The actual final normal profiles should be separately exported with timing labels; otherwise output cannot reconstruct the state used by the next step.

## Physical applicability issue

The current small-meltwater-feedback bound divides by total external inflow. If most water bypasses an almost blocked toe, this does not bound generated cold runoff relative to the toe stream. Local phase water may affect contact/mixing even when the global fraction is small. Report/gate route-specific generated flow or require an explicit drainage/separation assumption. This is especially relevant near a contact threshold. A global conservation pass is not a certificate that local hydraulic feedback is negligible.

## Decision

Hold the production freeze until numerical class/area refusals are fixed and the same frozen cases are rerun. Keep original failures and clearly distinguish expected applicability refusals from implementation bugs. No award, historical originality, environmental transfer or full-cliff accuracy conclusion follows from this review.
