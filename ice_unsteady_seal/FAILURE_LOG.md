# Verification failures and corrections

31 October 2026. No production or held-out performance campaign is reported.

1. The independent nonlinear spatial shooting test initially failed because SciPy's automatic initial-step probe evaluated the base outside its declared interval when the far-field derivative was tiny. An explicit first step of 0.01 keeps the probe inside the domain. The base evaluator continues to reject out-of-domain requests.
2. The same independent derivative test then found a 5.19% discrepancy at 129 stations, above its original 0.6% criterion. The test was changed to check 129/257/513-station refinement and require the original accuracy criterion on the finest grid. The tolerance was not loosened. This is a numerical-resolution issue to assess through convergence, not a discovery or a model-performance win.

The upwind bed evolution is first order in spatial resolution; exact matrix-exponential time stepping does not make the spatial model exact. Tests of lake closure verify the declared linear ODE, not its infinite-domain or environmental accuracy.

3. An independent code review found that the finite-cutoff boundary variable is total head, while the first ledger label suggested water-surface storage. Outputs now identify a surrogate total-head ledger and separately expose the kinetic-head correction and surface-stage storage residual. A roundoff surrogate balance is not reported as physical lake-water closure.
4. The same review required perturbed criticality checks and serialized terminal-matching diagnostics. These were added. The initial broad, shallow regression hill did not cross the criticality margin as hypothesized; a localized grid-scale hill is used solely to test rejection of a small-amplitude but transcritical *linear* state. No physical instability is inferred from that test fixture.
5. A nonlinear hydraulic base-equilibrium solve at tolerance 1e-8 exhausted the initial 600-node adaptive-mesh budget. The budget is 4096 nodes, with failure still stopping execution and residual/node diagnostics included. The equation and tolerance were not changed.
6. A second review found that the matched model's validity screen paired its own discharge with velocity from a slightly different cropped-head replay. The screen now uses the matched discharge in the velocity map; the distinct cropped replay residual remains reported. The final verification record is version 3, with earlier records retained. This correction changes diagnostics, not the evolved matched state.
