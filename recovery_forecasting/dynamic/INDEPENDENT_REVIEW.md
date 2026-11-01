# Independent reduced-model review

1 November 2026. Reviewed `model.py`, `test_model.py` and `cli.py` read-only.
No solver source, active process or existing result was changed by this review.
This reviews the version present during the audit; later fixes need rerunning
their own tests. The reviewer does not approve environmental or publication
claims.

## Verdict

The declared reservoir/parallel-control/plug-flow model is internally coherent
at the equation level. Its global sensible-plus-latent accounting is consistent
with **meltwater joining downstream**, rather than entering the reservoir.
It is a useful conditional time-dependent predictor and mechanistic baseline.
It is not a spatial flow solver, a moving-boundary CFD verification, or a
validated model of a full ice cliff.

## Mass and heat derivation

The reservoir states are liquid volume `V` and sensible energy `U` above the
melting point. With `T=U/(rho_w cp V)` the implemented equations are

    V' = I - q1 - q2,
    U' = rho_w cp I Tin - rho_w cp (q1+q2) T.

Each route receives sensible power `rho_w cp qj T`. Its passive plug-flow heat
sink is

    Pj = rho_w cp qj T [1-exp(-Gj Aj/(rho_w cp qj))],
    z_j' = -Pj/(rho_i L Aj).

For nonnegative inputs the mathematical sink satisfies
`0 <= Pj <= rho_w cp qj T`; downstream sensible energy is the remainder.
Adding reservoir sensible energy, downstream sensible energy and the route
latent sinks therefore recovers inlet sensible energy plus initial storage.
The ledger integrates this identity consistently.

Melting ice volume is `Aj*(z_j(0)-z_j(t))`; its mass is `rho_i` times this
volume, and downstream meltwater volume is `E_latent/(L*rho_w)`. Those are
consistent conversions for a temperate-ice latent model. The reservoir mass
ledger should **not** gain that meltwater while the physical topology explicitly
puts melting downstream. A separate domain-wide liquid-plus-ice mass identity
can make this scope clearer:

    Vfinal + Vout_reservoir + Vmelt_downstream
       - Vinitial - Vin - (rho_i/rho_w)*Vice_removed = 0.

This extra identity would be an accounting check, not independent validation of
the transport topology. There is no route water storage or density-jump velocity
solution because the outlets are quasi-steady ideal controls, not resolved
moving-interface flow domains.

## Recommended fixes/checks before calling outputs final reduced-model data

1. **Intermediate positivity — addressed during review.** The initial `rates()` evaluated
   constitutive laws using `max(U,0)` and `max(V,0)`, while only the final
   integrated state is checked. A negative intermediate or dense-output sensible
   state can therefore be hidden by an emitted zero temperature. The exact ODE
   should preserve nonnegative volume and energy. Check solver mesh and emitted
   states against declared tolerances; record tolerated roundoff and refuse
   significant violations. Do not silently label a clipped state physical.
   The revised source checks accepted solver states and emitted dense-output
   states, rejects `V<-1e-10 m3` or `U<-1e-6 J`, and reports raw U plus a
   roundoff-clipping flag. This resolves the reported silent-output issue for
   those checked states; it is not a proof of positivity of every internal RK
   stage or of a continuous interpolant between samples.
2. **Event physics — strengthened during review.** The initial test
   confirms exact inclusion of removal/forcing times. Also check continuity of
   `V,U,z` across each prescribed event and the expected right-limit change in
   flow/conductance. A terminal-time removal changes the reported final flux but
   has zero duration: it cannot establish sustained recovery. An explicit
   event ledger distinguishing left/right limits would help consumers. The
   revised tests compare state continuity and the right-limit flow jump against
   a no-removal control with an identical integration breakpoint. This is a
   useful direct check beyond timestamp inclusion.
3. **Add an independent reservoir-temperature limit.** At fixed equal outlets,
   zero transfer and a stage chosen so inflow equals total outflow, volume stays
   constant and `T(t)=Tin+(T0-Tin)*exp(-I*t/V)`. This tests the sensible mixing
   equation independently of its integrated algebraic ledger.
4. **Use route-relative latent tolerances.** The current phase mismatch is
   compared with a scale containing all initial/inlet sensible heat. When the
   latent sink is tiny, that scale can be large relative to actual melting.
   Report absolute residuals and normalize each route by its own latent budget
   plus a declared absolute floor. Similarly document units for scalar `atol`,
   which currently applies to both metres and joules.
5. **Report event states and output sampling limits.** Exact protocol boundaries
   do not mean contact thresholds or persistence start times are exactly
   resolved. Operational recovery scoring needs a common ordinary-flow return
   time, persistence window, and bracketing/refinement of wet/dry observations.
6. **Maintain the fixed-floor stop.** Termination when a moving control reaches
   the immutable pond floor is a scope refusal, not a successful forecast to
   the requested horizon. Add a test covering that refusal. Do not infer a
   natural incisional limit from this computational floor.

## Physical applicability limits

- The two ideal broad-crested controls assume a common well-mixed upstream
  stage and negligible route storage/inertia. They do not reproduce the current
  spatial entrance/partition gate by construction. Agreement with their own
  formula cannot independently pass that gate.
- Bed lowering changes the route crest while pond hypsometry remains fixed.
  This is coherent only when controls lie outside the fixed pond storage region
  and are represented as idealized outlets. It omits evolving route cross
  sections, lengths, contraction/entrance losses, backwater and channel widening.
- The same uniform deposit raises the hydraulic crest and adds a thermal series
  resistance. That is an explicit closure, not a treatment of under-rock gaps,
  porous flow, heterogeneous roughness or mobile grains. Complete removal is
  prescribed; debris mass, export work and transport rates are not predicted.
- No atmospheric/radiative input, frictional heating, cold-ice warming,
  refreezing, ice mechanics, fracture, lateral avulsion or moving cliff-face
  geometry is solved. There is no whole-cliff survival or glacier-wide mass-loss
  prediction from these states alone.
- A well-mixed pond storage/temperature model can add a genuine independent
  physical clock. Its resulting history response should be separated from
  quasi-steady thermal-dose behavior and from arbitrary debris-removal timing.

## What the current tests substantiate

The analytic isothermal drainage test is an independent hydraulic limit.
Symmetry, zero-transfer and tolerance refinement are useful conditional checks.
The global ledgers confirm compatible accounting; because their integrals use
the same RHS powers, they do not independently validate the heat-transfer law.
Input/hashing checks support reproducibility, not predictive accuracy.

The revised nine-test suite passed during this audit. Downstream meltwater
dilution is now explicitly exported and tested. These changes strengthen
physical bookkeeping without promoting the model into spatial CFD.

The current test suite supports **a conditional reduced-model implementation**.
It does not pass the failed spatial hydraulic partition gate, demonstrate a
new physical discovery, establish historical originality, or validate glacier
forecast improvement. The next useful comparison is frozen, held-out spatial
cases, with meaningful geometry-aware baselines and explicit scope matching.
