# Passive handoff — do not interrupt the active CFD milestone

This contribution is confined to `ice_cfd_checks/` on a separate Codex branch. It changes no existing file, build environment, OpenFOAM case, dependency list, job runner or result. It does not request that Claude restart or redirect any task in the user's pasted prompt. There are no solver installation steps here.

The active prompt and current source layout were read before writing this package. The source inspection found existing Graetz, step/ramp and budget-check scripts, so those were not edited or reimplemented as the main task. This contribution supplies independently coded analytical references, exported-data scorers and a separate future response-analysis study.

## Useful when current outputs become available

1. Review the additive commit after the current work reaches a safe boundary; no immediate cherry-pick is required for a running solver.
2. Run `python3 -m unittest discover -s ice_cfd_checks/tests -v`.
3. Generate the known examples with `python3 -m ice_cfd_checks reference --output <new directory>`.
4. Export the actual solver fields into the explicit schemas in README.md. Record any spatial/time averaging and source field definitions. This package does not yet parse arbitrary OpenFOAM files or function-object formats.
5. Score remote-section jump measurements, hydrostatic rest or thermal ledgers using the corresponding CLI. Define acceptance tolerances from the benchmark and uncertainty before interpreting actual results. Never label self-generated exact examples as independent CFD validation.

Do not compare raw p_rgh with physical hydrostatic pressure without reconstructing the correct pressure and reference convention. Do not substitute mixture flux for liquid flux in a VOF heat ledger. Prefer numerical face fluxes to alpha*U postprocessing where the solver transports them differently. For reversing flow, export the correct donor/face temperature for each sign.

The jump scorer uses unit profile factors and hydrostatic measurement sections. When velocity profiles, turbulence stress, friction or nonhydrostatic pressures matter, inspect the integrated momentum terms in the actual solver instead of treating the ideal residual as a coding failure.

The dam-break benchmark supplies a transient check, but it is an exact shallow-water reference. A short, deep VOF gate release includes vertical acceleration and air/free-surface dynamics outside that reference. Choose the appropriate asymptotic comparison or an independent VOF benchmark; the package cannot remove that modeling distinction.

The Stefan reference supplies target fronts, temperatures and energy budgets for equal-density, conduction-only melting. It does not exercise the unequal-density interface condition. Maintain the separate mass/interface test required by the user's prompt.

## Future scientific work

SCIENTIFIC_STUDY.md is preparation for the next research decision, not an instruction to begin another campaign during verification. The small linear-response analyzer evaluates a supplied physical linearization; it does not generate one from an unfinished solver. Its manufactured non-normal example is only a numerical test.

If the verified model later exports a validated reduced tangent operator, the tool can test whether finite-time thermal–morphological responses affect a declared incision observable beyond modal stability. That requires operator validation, careful scaling, conserved inputs, independent nonlinear checks, regime relevance and a full prior-art comparison.

No model improvement, new ice mechanism, field validity or novelty is established by the current package.
