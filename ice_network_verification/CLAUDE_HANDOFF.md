# Independent numerical contribution and handoff

## Scope and provenance

Created independently of the production solver, in an isolated worktree. Only the
new top-level `ice_network_verification/` directory needs to be cherry-picked.
No solver files, running jobs, chemistry outputs, or existing tests were changed.

The local branch starts at older repository commit
`0fc0096782badabfdbb0aab4392d30596acb7417`; it is intentionally additive. Do not
reset a current branch to that base. Cherry-pick the contribution or copy this
directory into the current worktree.

Inspected production snapshot:

`7b1d59ad9e2b685270cf5b028cb77b3f69fb07d3` on `claude/ice-network-m1`.

Inspected the mathematical contract and exact thermal/constants files from that
pinned commit via the GitHub connector. The separate source checkout used to run
the audit is not vendored into this package. `records/claude_heat_kernel_audit.json`
records SHA-256 hashes and all 21 scenario/resolution results.

## Completed work

- Seven exact limiting-case specifications, SI units and requested times.
- Solver-export checker with explicit missing-data and malformed-data failures.
- Independently recomputed kg and thermal-energy ledgers.
- Reference-error convergence-order utility.
- An independent heat-kernel numerical audit, with 60-digit Decimal references.
- A stable alternative segment utility (`thermal.temperate_segment`).
- 26 passing tests, including deliberately wrong density, wrong outlet split,
  zero-discharge, non-finite, missing, misaligned and double-counted-energy inputs.

Identity-fixture tests exercise the checker. They are not results from Claude's
hydraulic solver. The full production solver has **not** been run against the
seven limiting cases; it needs an adapter where compatible.

## Actual findings on the production thermal kernel

Five ordinary scenarios were checked at 1, 10 and 100 uniform segments. Maximum
absolute heat-to-ice discrepancy from the high-precision ODE solution was
`1.3969838619232178e-9 W`. This supports numerical correctness of the inspected
segment law in these cases. It does not validate Dittus-Boelter for open channels.

Two additional numerical stress scenarios set `nu_factor=1e-12`. They are NOT
plausible calibrated heat-transfer settings. The original equilibrium-temperature
formula suffers catastrophic cancellation at thermal exposure near `1e-13`:

- One warm one-segment case yields `-0.0117034912109375 W` heat to ice instead of
  approximately `8.76637e-10 W`.
- One cold one-segment case yields the same negative heat instead of approximately
  `1.03818e-12 W`.
- The largest absolute error over these stress cases/resolutions is about
  `0.696106 W`.
- The algebraic thermal budget can still close exactly while the solution is wrong.

These stress results do **not** show that the default coupled simulations are
wrong. They identify a small-exposure stability defect relevant if the code later
explores zero-transfer limits or very small per-segment exposure.

## Stable algebra (no clipping)

For capacity rate `C=rho cp Q`, conductance `UA=h_T P dx`, source `Phi`,
`a=UA/C`, and excess inlet temperature `theta=Tin-Tm`:

`F = 1-exp(-a)`

`Tout = Tm + theta*exp(-a) + (Phi/C)*(F/a)`

`q_ice = C*theta*F + Phi*(1-F/a)`.

Use `-expm1(-a)` for F. For small a evaluate
`1-F/a = a/2-a^2/6+a^3/24-a^4/120+...` without subtracting near-equal numbers.
At a=0, Tout=Tin+Phi/C and q_ice=0. The independent utility implements these limits
and rejects refreezing/non-finite inputs rather than clipping negative outputs.

Against all 21 high-precision references, the stable utility's maximum heat error
was `4.440892098500626e-16 W`. This comparison evaluates the whole uniform reach;
production integration still needs segment-by-segment tests, including nonuniform
geometry and localized entrance dissipation. Do not swap it in blindly.

## Recommended integration

1. Run the 26 tests. The core suite has no third-party dependencies.
2. Review and optionally port the stable algebra into the production thermal kernel
   with the current constants and source-placement convention. Keep source and
   conductance definitions identical. Check Numba compilation separately.
3. Add test adapters for analytical cases compatible with the public solver API.
   Linear edge laws are test fixtures, not replacements for Manning/GVF physics.
4. Independently compare geometry-derived ice mass with latent-heat-derived mass;
   the global water ledger alone cannot validate geometry updates.
5. Keep nominal-pipe Reynolds-range flags descriptive: even high Reynolds numbers
   do not validate a smooth-pipe correlation for an open melting channel.

## Remaining physical issues from the reported milestone

- The open-channel heat-transfer closure is unvalidated. Sensitivity to a multiplier
  is useful but is not a substitute for experimental validation.
- Quasi-steady channel-storage neglect should receive a diagnostic, not only an
  assertion: report transit-time/stage-change scales and meltwater-to-flow ratios.
- Multiple entrance roots at very small Q need a stated branch-selection rule.
  A small instantaneous difference may still alter future topology activation;
  compare integrated effects and activation thresholds before calling it negligible.
- First-order geometry updates have reported non-negligible incision errors. Use
  refined configurations or quantified uncertainty for scientific comparisons.
- SaDS/SHED differences remain hypotheses until matched comparisons run. This
  contribution establishes neither field predictive improvement nor novelty.

## Short message for the user to pass to Claude

Review `ice_network_verification/` from `codex/ice-network-verification`. It contains
independent analytical cases and a pinned thermal-kernel audit. Ordinary heat cases
agree closely; an extreme weak-transfer stress test reveals cancellation, with a
stable alternative provided. Integrate only after checking source placement and
Numba behavior. Continue the existing physical-validation and prior-art work.
