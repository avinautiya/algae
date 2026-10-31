# Unsteady ice-seal kinetics and state-dependent lake discharge

This independent addition develops a substantially larger model than the earlier retention identity. It contains continuum derivations, an infinite-pond hydraulic matching result, a canonical linear state/history law with the outer translation removed, and executable finite-interval linear and nonlinear forecasters.

The research target is specific: **derive discharge during rapid evolution of a lake's ice seal, where an instantaneous outlet law can omit the evolving inner geometry**. Schoof et al. explicitly identify this fast inner-region consistency problem in [JFM 2023, section 3.5](https://doi.org/10.1017/jfm.2023.130). Identifying that gap is not evidence that every present glacier model lacks time-dependent hydraulics.

## What this addition actually establishes

The governing extension retains the fast bed-shape state B(X,T), with hydraulic velocity and discharge selected by a spatial boundary problem. Discharge is determined by current head **and** that state. No relaxation time is fitted to manufacture a delay.

For a canonical nonuniform seal, the derivations establish:

1. A complete base profile, with algebraic approach to the upstream pond and exponential approach to downstream normal flow.
2. A positive hydraulic normalization that remains finite as the upstream boundary moves into the pond. Its leading upstream truncation correction is proportional to L^(-4).
3. A nonzero sensitivity of discharge to interior bed shape even when the upstream head, outer slopes and crest geometry are unchanged. This sensitivity persists in the hydraulic matching limit.
4. A well-posed full-line linear state evolution. Removing the outer translation/datum mode yields a closed shape/discharge model with a derived, nonzero renewal kernel. The smooth inner crest and outer line-intersection seal height remain distinct observables.
5. Conditional slow-limit and stability results stated in STABILITY_AND_SLOW_LIMIT.md. These recover the existing canonical threshold; no new oscillatory instability is invented.
6. A non-exponential-memory lower bound in the ideal infinite-pond model: the renewal kernel is at least a positive constant times T^(-6) at sufficiently large inner times. This is a lower bound, not an exact tail law or an observed glacier timescale; finite ponds and the ordering of the inner/outer limits matter. See MEMORY_TAIL.md.

The mathematical tools include asymptotic analysis, positive integral operators, transport semigroups, implicit hydraulic selection and Volterra equations. Those tools are established. Historical priority of the **specific model and derived results** has not been certified. The three agents independently examined prior work, derivation and rejection tests; they are not a substitute for external expert review.

## Executable models

All inputs and outputs are dimensionless. These are mathematical prototypes, not site-calibrated glacier forecasts.

| Engine | Implemented state and output | Principal boundary |
|---|---|---|
| `matched` | Canonical linear outer-relative head, internal shape, discharge, outer seal shift, inner crest/head, derived renewal kernel. | Discretizes the proved real-line linear law using finite cutoffs. Cutoff/quadrature defects are reported. |
| `linear` | Fixed-frame bed profile and total-head state, discharge, linear crest shift, matrix-exponential propagation. | General finite positive-velocity subcritical interval with imposed upstream bed perturbation and terminal normal-flow impedance. |
| `nonlinear` | Finite bed-profile evolution and repeated nonlinear hydraulic BVP solves. | Fixed frame, wet/subcritical branch only; no outer nonlinear shock-branch selection. |

```sh
python3 -m unittest discover -s ice_unsteady_seal/tests -v
python3 -m ice_unsteady_seal demo-case --output /tmp/ice_seal_case_01.json
python3 -m ice_unsteady_seal predict --engine matched --case /tmp/ice_seal_case_01.json --output /tmp/ice_seal_matched_01.json
python3 -m ice_unsteady_seal predict --engine linear --case /tmp/ice_seal_case_01.json --output /tmp/ice_seal_linear_01.json
python3 -m ice_unsteady_seal predict --engine nonlinear --case /tmp/ice_seal_case_01.json --output /tmp/ice_seal_nonlinear_01.json
```

Output paths are created exclusively; existing results are not overwritten. Known future forcing is an explicit supplied scenario, not a future observation silently treated as a forecast input. Negative water inputs, incompatible shapes, dry states and invalid hydraulic branches are rejected or flagged. The linear screen is an amplitude/criticality diagnostic, not an error bound.

The fixed-frame finite-cutoff variable is **total head**. Its exact surrogate bookkeeping is not physical lake-water closure. The channel surface-stage correction is reported separately. In the matched infinite-pond theory total head and pond stage coincide; finite-cutoff, interpolation and dynamic matching errors still matter. Neither a tiny ledger residual nor linear/nonlinear agreement establishes environmental accuracy.

## Verification and figures

The verification suite checks independent quadrature, nonlinear spatial shooting, adaptive nonlinear BVP derivatives, independent time integration, the uniform moving-gauge cancellation, positive shape-generator structure, refined translation covariance, criticality rejection and declared storage bookkeeping. Unit records, manufactured CLI cases and versioned verification data are in `records/`.

Three figure sets have PNG and vector SVG exports with underlying JSON and hashes: base/spatial/pond-limit checks, total-head versus stage bookkeeping, and the matched shape kernel/control-height distinction. All are **verification figures**, not withheld forecast-skill results.

The first source record is retained as `verification.json`; version 2 added the matched engine, and `verification_v3.json` incorporates its reviewed diagnostic correction. The current `figures/` use version 3; earlier figure exports are retained in `figures_v2/`. This is numerical-verification versioning, not a preregistered scientific performance experiment. See FAILURE_LOG.md for failed checks and fixes.

## What is not yet established

The full nonlinear moving-seal/outer matching problem, nonlinear branch basins and large breach/reseal transitions remain open here. The code does not resolve water-wave propagation, hydraulic jumps, cold-ice conduction, solar heating, an independent turbulent heat-transfer law, hydrofracture, sediment or whole-glacier mechanics. The retained incision closure is dissipation-driven; this work does not claim a newly discovered heat-transfer mechanism.

The previous 120-profile heat benchmark is unchanged and still favors its classical two-mode baseline. This model addresses a different physical prediction. No superiority over existing coupled glacier models, Himalayan transfer, satellite correction, agency readiness, Nature suitability or award outcome has been demonstrated.

## Reading map

- [MATHEMATICS.md](MATHEMATICS.md): scaling, governing equations, hydraulic operator and gauges.
- [MATCHING_THEOREM.md](MATCHING_THEOREM.md): infinite-pond hydraulic matching and full-line linear well-posedness.
- [MATCHED_LINEAR_CLOSURE.md](MATCHED_LINEAR_CLOSURE.md): outer-matched shape/discharge law and nonzero kernel.
- [STABILITY_AND_SLOW_LIMIT.md](STABILITY_AND_SLOW_LIMIT.md): the canonical static family and stability qualifications.
- [MEMORY_TAIL.md](MEMORY_TAIL.md): the algebraic lower bound, finite exponential-closure limitation and physical ordering caveat.
- [NONLINEAR_MODEL.md](NONLINEAR_MODEL.md): finite nonlinear equations and local solution criterion.
- [PRIOR_ART.md](PRIOR_ART.md): close literature, actual reading scopes and incomplete priority coverage.
- [RED_TEAM.md](RED_TEAM.md): independent rejection tests and physical limits.
- [HANDOFF.md](HANDOFF.md): concrete review instructions without disturbing Claude's current work.
