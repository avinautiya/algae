# Independent proof and numerical evidence review

Pinned repository commit: `455bf85856bcf11771f0b37bd63cc58e89f05f32` on `claude/recovery-m5`. Seven source files were fetched and every snapshot's Git blob SHA verified. No Claude file or frozen result changed. Six independent tests pass; no expensive campaign was run.

## Verdict

The rigorous21-day no-go survives. The review found correctable reporting and parameter-path errors, plus a spatial-refinement design problem; none supplies evidence of sustained recovery or superior prediction.

## Corrections required

1. **Not every area-sweep case has zero repair.** `viability_check.json` already records positive repair for narrow toe footprints. AtQ=0.14m³/s andAtoe=0.001m², independent RK2 integration reproduces0.1034348mm repair in21days, final gap0.249896565m. AtAtoe=0.01m², the saved result has0.08495mm repair. The statement that the gap always grows throughout the adversarial sweep is false. The required repair is61.076mm, so the sustained-recovery exclusion remains valid. The declared12base cases can still all show zero repair; distinguish those from the expanded sweep.
2. **The164×NTU statement holds only on a fixed-film-limit parameter path.** KeepinghT fixed gives requiredNTU0.221714 versus declared0.00135506, ratio163.619. Ifh itself varies, bothNTU and the film-limit prefactor change. Solving the full bound givesh≈2474.027Wm⁻²K⁻¹, ratio12.370,NTU≈0.0167623 as the point removing this exclusion. That does not establish actual recovery, independently justify thath, or permit tuning it for a win. Equal-area symmetric models still forbid reset by their own monotonicity argument.
3. **The spatial series changes physical geometry.** The actual Sim constructor usesny=2round(Wl/dx)+1 and a one-cell divider. WithWl0.25m anddx0.25/0.125/0.0625m, divider widths are0.25/0.125/0.0625m and total upstream pool widths0.75/0.625/0.5625m. Thus toe-flow differences conflate discretization with a changing pool/divider. They are not a fixed-domain grid error estimate or valid convergence order. A future study must hold physical domain, divider, boundary injection and ROI fixed with consistently rasterized masks. Merely adding finer meshes or prescribing four cells across a lane does not resolve this confounder. `require_fixed_geometry` provides a fail-closed descriptor gate; it does not modify the solver or pretend that geometry has been corrected.
4. **Global defect headlines understate saved maxima.** Across all controls, per-step mass maximum is7.23461e-13kg, sensible maximum7.60701e-9J and sensible-plus-latent maximum7.60519e-9J. The6.2e-13/5.1e-9 figures describe the spatial subset, not every tested control. These are still small numerical residuals; changing the headline does not establish prediction accuracy.
5. **Baseline equivalence remains a hard claim ceiling.** The candidate and same-equations evolving-topography comparator are identical in this ideal lane. Agreement cannot support improvement. Frozen class-label mismatches must remain visible; another agent reviews their exact endpoint semantics.

## What the proof does establish

For g>=0 and equal hydraulic rating coefficients, bypass discharge is at leastQ/2. With common fixedh,T, independent route supplies, floor-only heat and no other geometry terms, toe incision is no larger than the film limit and bypass incision is at least its finite-water value atQ/2. Integration gives the conservative repair upper bound. The four pulse multiples permit only0.40119/0.20064/0.10033/0.05017mm repair in21days, all below61.076mm. This necessary-bound result is not a theorem about arbitrary natural glaciers or the resolved SWE with different thermal paths, sources, storage, temperature or geometry.

## Reproduce

Run `python3 check_proof.py` and `python3 -m unittest -v test_proof` in this folder. `CHECK_RESULTS.json` records arithmetic, the narrow-area counterexample, path-dependentNTU results and saved numerical maxima. No new parameters were fit to any test data.
