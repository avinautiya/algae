# Passive future handoff

Continue Claude's existing prompt to completion. This package is offline research/code on a separate branch. It changes none of the active build, case, thermal solver, geometry or job files. It launches no CFD campaign and needs no environment changes.

The sole code dependency is the project's independent `ice_cfd_checks` contribution (commit `65819f5cf101f485e061865b86a13f9a56d8c0b9`). Both packages are included on the branch and in the delivery ZIP. The new commit can be reviewed after the current milestone reaches a safe boundary. Do not import unrelated history or restart running jobs just to use it.

When ready:

1. Run the new package's tests and synthetic commands in README.md.
2. Read LITERATURE_AUDIT.md before describing the candidate mechanism as new.
3. Use PHYSICS_SCOPE.md to assess the intended environmental regime. A small bulk Richardson ratio alone is not a proof of negligible buoyancy in a recirculation zone. Conversely, a synthetic local ratio is not evidence that such a zone exists in the real case.
4. If a physical reduced operator later becomes available, record its basis, units, normalization, validity, observation map and conservation constraints. Validate derivatives and nonlinear response before replacing the manufactured optimizer input.
5. Use the optimizer to select a bounded, physically feasible perturbation experiment. Check exact nonlinear water/heat budgets, temperature/discharge bounds and time-grid convergence. Do not interpret an unconstrained signed perturbation as a realizable meteorological history.
6. Predeclare a forecast observable and held-out cases. Compare appropriate physical/empirical baselines at matched input information. The old event-averaging findings remain part of the record.

The active solver's success/failure logs must remain the source for its runtime status. This package has neither diagnosed nor repaired the current interFoam execution failure. It provides theory and later analysis support without entering that workstream.
