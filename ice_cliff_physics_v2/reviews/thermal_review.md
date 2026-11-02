# Independent thermal integration review — post-fix result

**Pass for the reviewed implementation scope. No environmental accuracy, historical novelty or publication-level claim follows.**

## Scope and reproducibility

The reviewer did not edit parent integration files. Twenty seeded synthetic warm configurations were defined before execution; exactly the same configuration objects were rerun after the parent fixes. The original failing run, source hashes, frozen cases and exact failure traces are retained in the *_pre_fix artifacts. Source hashes are unchanged during each individual run. The parent implements any code changes; their effects are explicitly versioned here.

## Post-fix execution

- 17 of 20 completed; the other 3 are declared physical applicability refusals, not silently dropped runs.
- No unexpected failure remains. Maximum normalized global energy residual: 4.42e-16 of gross incoming/atmospheric heat.
- Maximum completed geometric phase-volume residual: 3.86e-15 m³; water bookkeeping residual: 1.2e-09 m³. Water bookkeeping is not an independent validation of the mass model.
- Cold-reservoir wet-debris input is refused. Zero external inflow with atmospheric melt is refused without a declared drainage contract and accepted under its stated synthetic free-drainage assumptions.

Explicit refusals:

- T02: local route meltwater feedback requires a declared downstream free-drainage contract.
- T12: local route meltwater feedback requires a declared downstream free-drainage contract.
- T18: local route meltwater feedback requires a declared downstream free-drainage contract.

## Regression fixtures for the four original numerical failures

The new local-runoff applicability guard can stop a configuration before its earlier numerical failure. To test fixes directly, the four formerly failing configurations were also frozen with an explicit separated downstream runoff contract and a prescribed 0.001 m³/s drainage capacity. This is an implementation regression assumption, not new environmental evidence. All four complete:

| Case | Normalized global energy residual | Max post-geometry profile/scalar enthalpy mismatch | Actual maximum height bands |
|---|---:|---:|---:|
| T11 | 4.28e-16 | 3.38e-16 | 0 |
| T12 | 2.54e-16 | 3.62e-16 | 0 |
| T16 | 4.57e-17 | 1.93e-16 | 128 |
| T19 | 1.75e-17 | 2.19e-16 | 70 |

Every exported canonical column area exactly agrees with its corresponding post-geometry scalar reservoir area. The largest normalized post-geometry enthalpy mismatch is3.62e-16.

## Corrected issues

1. Canonical total area now survives wet/dry class reclassification. A one-ULP sum of rounded class areas no longer causes a false contraction refusal. True requested contraction remains refused; no outgoing enthalpy is silently discarded.
2. Heat conductance areas use the authoritative geometric class partition. A mathematically zero bare-water class has exactly zero water conductance, avoiding an unallocatable phantom heat source.
3. Temperate phase accounting preserves supplied latent energy exactly, and resolved dry-column heat is used directly instead of subtracting nearly equal wet/total heats. A zero dry class no longer receives spurious ULP-scale energy.
4. The actual final normal profiles are exported separately as debris_columns_post_geometry, with sampling times, instead of exposing only the post-heat/pre-geometry profile.
5. The meltwater applicability policy now checks toe and bypass route ratios, not only global inflow. High local feedback requires a declared drainage/separation assumption. This necessary ratio guard is not a proof that threshold/contact predictions are insensitive to a five-percent perturbation.

## Physics accounting inspected

The upstream floor removes heat once; bare and covered face patches share the remaining stream. Bypass water receives its independent routed share. Covered-interface heat enters phase change once. Distinct wet/dry temperature profiles persist between reclassification events; the stored global mean is only an energy projection. New debris area imports the prescribed normal-profile enthalpy. These checks establish implementation behavior for the declared closures, not that real boulders or cliffs obey them.

prescribed_all_cover supplies the declared prepared external flux to both wet and dry covered areas; that input must exclude the separately modelled water-convective exchange. fixed_cover_interval_m is an Eulerian height mask. Relative-cover modes instead prescribe cover redistribution as the floor lowers; they do not implement mobile debris or material-mask advection.

## Remaining scientific limits

- The pressure/contact rating and Manning slope remain imposed reduced closures with only a necessary head check, not resolved environmental flow.
- Normal slab diffusion resolves thickness, not full lateral wet/dry heat conduction. Representative donor averaging on reclassification cannot recover a face-height temperature correlation.
- Height-band face retreat uses an explicit adiabatic-ledge approximation and does not predict structural collapse.
- Adaptive class alignment changes the actual number of height bands as wetting history and time step change. Initial band count alone is not a mesh-resolution certificate, and time refinement may also change spatial representation.
- Source heat-transfer correlations and imposed temperatures still need independent environmental validation; exact conservation does not establish field accuracy.
- The present review is synthetic development/regression testing. It is not a blind test and must not be counted as independent glacier observations.

## Decision

The reviewed fixes remove the four identified numerical blockers and close the missing output-state audit. A bounded frozen campaign may proceed under the declared applicability envelope, preserving every refused case and comparing all treatments over matching supported cases. Numerical refinement and independent physical comparison remain separate gates.
