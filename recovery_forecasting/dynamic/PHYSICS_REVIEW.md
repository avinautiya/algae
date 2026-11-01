# Independent physical review of the dynamic baseline

Read-only review of `model.py` and `README.md`, 1 November 2026. No active
solver source, input or result was changed. Numerical state/event handling was
reviewed separately in `INDEPENDENT_REVIEW.md`.

## Verdict

The equations form an internally consistent **conditional reduced model** under
their declared contract: a fixed-footprint well-mixed reservoir feeds two ideal
outlet controls; sensible heat advects out of the reservoir, some becomes
latent ice melt downstream, and the resulting meltwater is inserted downstream
at the melting point without feeding back on the controls. The model's water
and heat ledgers are compatible with that contract. They are not evidence of
environmental accuracy, closure validity, or independent hydraulic validation.

## Heat and phase accounting

The reservoir loses `rho_w*cp*sum(Q)*T`, while each downstream outlet carries
that route's sensible energy minus its extracted power. Summing downstream
sensible energy, latent energy and reservoir storage therefore equals upstream
energy by the same RHS algebra. Ice lowering and latent energy also share the
same extracted power. This is the correct verification identity for the chosen
approximation; integrating its components in the same solver will preserve a
linear invariant even if the trajectory itself is inaccurate.

The downstream meltwater flow is `P/(rho_w*L)`, not `P/(rho_i*L)`. The existing
formula uses the correct liquid density. Ice volume lowering separately uses
`rho_i`. Adding zero-superheat water downstream dilutes the reported outlet
temperature; no extra sensible-heat source is introduced.

## Physical applicability gates to retain

1. **Neglected hydraulic/thermal feedback of meltwater.** The reservoir-route
   control and plug-flow heat calculation use the pre-melt route discharge.
   Meltwater is inserted afterward. Report the maximum dimensionless added
   fraction `P/(rho_w*L*Q)` over wet routes. The energy limit gives an optimistic
   ceiling `cp*T/L`, but actual accepted runs should report their actual ratio.
   If that ratio is not small relative to the required accuracy, this
   approximation cannot validate a quantitative recovery boundary. Real
   distributed phase addition also changes downstream residence and cooling.

2. **Time-scale separation.** A dynamic reservoir does not make the route flow
   or route temperature dynamic: both are instantaneous quasi-steady closures.
   Reservoir mixing time and route hydraulic/thermal residence time must be
   short compared with input variation and geometric evolution. The current
   configuration lacks enough information to establish those bounds. A sharp
   input step is a controlled reduced-model protocol, not evidence that actual
   water mixes or propagates instantaneously.

3. **Critical-control applicability.** Ideal broad-crested discharges omit
   entrance loss, tailwater/backwater coupling, lateral exchange and inertia.
   The spatial solver's unresolved partition/reference discrepancy cannot be
   settled by treating this same law as an independent exact oracle for abrupt
   obstacles. Use matched smooth-control verification or a better reference.

4. **Deposit geometry.** Raising a crest by one uniform deposit thickness and
   adding its thermal series resistance represents fully covering material.
   It does not resolve open pores, coarse clasts, under-rock passage or bare
   corridors. Do not use it as proof that complete coarse-debris clearing is
   the only way to restore flow. The external spatial campaign must test those
   alternatives. Prescribed removal is not a sediment-transport prediction.

5. **Fixed hypsometry and external controls.** The reservoir contains no moving
   ice in this model. Its outlet ice areas are separate downstream footprints.
   If actual incision changes pond storage area, if reservoir displacement
   occurs during removal, or if pond-floor melt contributes water, the current
   contract no longer represents that geometry without additional terms.

6. **Liquid water and temperate ice regime.** Temperatures are superheat above
   freezing, not absolute kelvin. Validation currently enforces nonnegative
   superheat but does not enforce a physically calibrated upper range. Constant
   water properties, no boiling and no cold-ice sensible warming must remain
   explicit assumptions. Large supplied temperatures can be mathematically
   accepted while physically inappropriate.

7. **Tolerance interpretation.** One scalar absolute tolerance covers volume,
   energy and elevation states of different dimensions. Its default is a
   numerical choice, not a universal physical error bound. Report geometry,
   contact timing and local heat convergence under independently varied solver
   tolerance; ledger closure alone cannot replace that study. Near a dry-route
   activation threshold, uncertainty in geometry and control coefficients may
   dominate a very small ODE integration error.

## What this can establish now

The baseline can compute a response to a supplied hydrograph inside a stated
reduced model. It improves the *information content* of the earlier static
feasibility test by retaining water/heat storage and the actual forcing order.
That does not by itself establish improved prediction over a spatial model or
measurements. Use it as a strong competing forecast, a causal-control tool and
an inexpensive debugging target.

Recovery requires a separate operational definition after return to ordinary
flow. `toe_fraction` is the toe share of reservoir outflow, not the fraction of
the added-meltwater downstream flow and not a cliff-survival probability. A
small nonzero toe fraction does not imply useful undercutting or persistence.

The most defensible next test remains post-clearance recovery under unchanged
versus previously incised substrate, with porous and lateral pathways retained,
and with numerical/closure uncertainty smaller than the proposed improvement.
