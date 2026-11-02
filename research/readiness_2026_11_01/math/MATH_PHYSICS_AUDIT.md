# Mathematical and physical readiness audit

Pinned source: `avinautiya/algae`, branch `claude/integration-m4`, commit `05e7312e3bbf248acf9f646267a5b8ea06a8b32e`. Read-only audit; no Claude worktree, solver, checkpoint or frozen result edited.

## Executive verdict

**No-go for a new mathematical discovery, unconditional irreversibility, full coupled mass conservation, or field superiority.** A bounded computational mechanism remains feasible, subject to prior-art and stronger comparisons. The actual equations already explain the principal fixed-forcing result through an elementary hydraulic threshold. A sharper future question is what forcing can produce *durable* rather than temporary recovery, but this audit does not establish its historical originality or usefulness in nature.

## Exact ideal-control result and proof

The imported `compute/topology.py` and wrapper `reactivation/reduced.py` use

    q_t = K_t (eta - b_t - d)_+^(3/2)
    q_b = K_b (eta - b_b)_+^(3/2),    q_t+q_b=Q.

Here b denotes ice-control elevation, d is a prescribed obstacle height, K has units m^(1/2)/s, and Q is discharge. Let g=b_t-b_b and H_c(Q)=(Q/K_b)^(2/3). If q_t=0, q_b=Q, so eta=b_b+H_c; hence

    toe dry iff g+d >= H_c(Q).

At equality the ideal toe discharge is zero. Monotonicity of total flow in eta gives the converse. Removing the obstacle changes d to zero: removal alone cannot restore flow when g>=H_c. This is an algebraic consequence of standard critical controls, not a new theorem requiring a sophisticated model.

Stefan lowering in the same code gives db_i/dt=-E_i, where

    E_i = P_i/(rho_i L A_i),
    P_i = rho_w c_w q_i T_in [1-exp(-h_i A_i/(rho_w c_w q_i))].

Therefore gdot=E_b-E_t. Under fixed Q and the present nonnegative-heat, no-other-geometry assumptions, a dry toe has E_t=0; if bypass melting persists, g increases. The dry state is absorbing **under these assumptions**. No ice deformation, lateral capture, channel migration, collapse, seepage, atmospheric melting or water-level changes are permitted to undo it. This is geometric memory of differential thermally driven incision. The current SWE mechanism does not include debris thermal storage, so it should not be presented as a demonstrated thermal-storage hysteresis discovery.

For the actual ideal inputs Q=0.07 m³/s, W_b=0.5 m and K_b=(2/3)^(3/2) sqrt(g_acc) W_b, K_b=0.852447 m^(1/2)/s and H_c=0.188924 m. Taking bare bypass melt rate h*T/(rho_i L)=0.028210 m/day yields an elementary removal-time estimate (H_c+0.08)/rate=9.533 days, close to the reduced prediction near9.5days. This is an explanatory approximation, not independent validation of the final CFD boundary around10.4days.

## Executed falsifier: reopening is not durable reset

`falsify_irreversible.py` independently solves the ideal algebraic controls, without importing the project. Eight gap-threshold checks pass. For a residual gap0.25m:

- At Q0=0.07 m³/s, q_toe=0.
- The reopening threshold is Q*=K_b g^(3/2)=0.106556 m³/s, about1.52Q0.
- A pulse to2Q0 gives q_toe=0.006834 m³/s.
- Returning to Q0 with unchanged geometry gives q_toe=0 again.

A discharge pulse can therefore falsify *unconditional* no-reopening language, but this example does not demonstrate durable recovery. For a pulse ending at time T, return-to-baseline contact requires

    g(T)<H_c(Q0)
    iff integral_0^T (E_t-E_b) dt > g(0)-H_c(Q0),

assuming unchanged baseline rating and no additional geometric terms. This is an exact identity/criterion, not a newly claimed mathematical result. Instantaneous wetting is necessary to deliver water heat in this model but is not sufficient to satisfy the integral.

For the example, the gap repair required is0.061076m. With declared h=200 W m^-2 K^-1 and Tin=0.5K, even an optimistic bound E_t<=hTin/(rho_i L), E_b>=0 permits at most1.175mm repair in one hour. The minimum repair time under that optimistic bound exceeds2.165days. The bound is conditional on fixedh/T and no additional heat sources; it is not universal for floods or nature.

## Additional conditional impossibility check

For identical bare heat-transfer coefficients and inlet temperature, define s=q/A and e(s)=rho_w c_w T s/(rho_i L)*(1-exp(-h/(rho_w c_w s))). For x=h/(rho_w c_w s)>0,

    e'(s)=rho_w c_w T/(rho_i L) * [1-(1+x)exp(-x)] >0.

With equal route widths, areas and ideal weir coefficients, a positive gap g gives q_t<q_b and thus E_t<E_b. The gap cannot shrink under any positive Q(t); pulses may temporarily overtop but cannot durably erase the gap within that symmetric model. This is a conditional diagnostic, not a claim about the full SWE. **The imported strip model uses A_t=0.5m² and A_b=1m²**, so the equal-area statement does not apply to it. Any proposed reset campaign must specify thermal footprints and the physical process capable of making E_t>E_b, rather than assume a sufficiently large pulse must heal geometry.

## Concrete defects / qualifications, ordered

### P0: readiness checker is not a typed evidence validator

The exact committed `ice_integration/release/validity_contract.py` rejects only values equal to None, False, empty string or 'unknown'. `validator_counterexamples.py` demonstrates seven invalid templates that return `may_issue_environmental_prediction=True`: string'false' independence, support fraction2, negative uncertainty, NaN uncertainty, empty source dictionary, empty survey-date list, and all-True placeholders without quantitative/source evidence. The promised numerics group is absent from `REQUIRED`. See `VALIDATOR_COUNTEREXAMPLES.json`. Typed finite domains and traceable source bindings are needed; even a correct schema only checks supplied evidence, not scientific truth.

### P1: full ice+liquid mass coupling is absent

`reactivation/swe2d.py:Sim.lower_bed` updates z_b but leaves h, hu, hv and hT unchanged. No meltwater mass/momentum source is added. The fluid advection ledger can close to1e-12 while the full water+ice phase-mass budget omits rho_i*DeltaV_ice. In an isolated stationary-ice source step with no drainage, physical updates are

    z_new=z_old-dz,  h_new=h_old+(rho_i/rho_w)dz,

with unchanged depth-integrated horizontal momentum and sensible enthalpy if meltwater is born stationary at0°C. Free-surface change is -(1-rho_i/rho_w)dz, not -dz. This operator is not a full pressure-work/ALE solver. The source may be globally small at0.07m³/s; its effect must be bounded per route near disconnection before treating the one-way approximation as adequate. No claim here that this omission explains the failed field comparison.

### P1: comparisons inherit physics and unequal structural information

Both reduced and SWE calculations use the same supplied series thermal resistance and assumed h/T. Their agreement tests aspects of routing/evolution, not the true closure. The fractional baseline changes both initial crest and shielding by averaging, so large timing error is not automatically attributable solely to missing dynamic connectivity. Add an initial-head/flow-matched baseline if isolating the benefit of subsequent geometric evolution. The fixed-geometry model cannot predict incision-driven contact loss by construction.

### P1: release memo cites superseded mechanism summaries

`DECISION_MEMO.md` points to `REPORT_M3_FINAL.md`, which still contains the earlier4/4 recovery wording, duplicate pre-removal trajectories, scope-mismatched melt errors, and v1 early-removal behavior. `REPORT_M3_ADDENDUM.md` and `SOLVER_VERSIONS.md` correct those. A final evidence table must identify exact solver/run/metric versions. The0.02day repeated loss prediction acrossE04/E06/E07 is one shared history, not three independent accuracy successes. The±0.37day GCI estimate is conditional on a particular grid/geometry/operator, not a statistical interval or proof of model accuracy.

### P2: counterfactual necessity is conditional

Uniform-incision `lower_bed` redistributes each run's instantaneous total dz over all its ice cells. This preserves that step's volume **within the same run**, not necessarily the whole time-integrated melt across full versus uniform evolving runs. It introduces deliberately nonlocal heat redistribution. Exact two-control gdot=0 is a useful logical ablation, but four full-SWE cases without loss establish neither universal necessity nor that nature realizes this intervention.

### P2: calibration influence is not a causal error attribution

The weighted one-scale least-squares influence and leave-one-out changes demonstrate instability. They do not establish that calibration leverage causes the test bias. The effect of width is evaluated after per-scenario refitting and hence includes compensation, not a clean fixed-parameter causal width effect. Tested joint81-scenario bias remains positive, which rules out the *tested combined envelope*, but not untested interacting uncertainties or misspecification. Individual marginal factor shifts must not be summed as independent error bounds.

## Decisive next falsifier, only if literature supports a gap

Freeze a temporary-reopening versus durable-reset prediction from the actual equations, including Q(t), T(t), removal schedule, identical physical footprints, all local heat sources and a well-defined post-pulse baseline window. First establish whether the chosen model can satisfy the repair integral at all. Then use an independently discretized resolved comparison with thermal closures varied and matched initial hydraulics. Refusal, temporary wetting, permanent loss and durable recovery are distinct outcomes; retain failures. This could test a quantitative environmental-recovery hypothesis, but historical novelty and field significance are not established by this audit.

## Deliverables

- `falsify_irreversible.py`, `FALSIFIER_RESULT.json`: executed hydraulic counterexample and necessary repair bound.
- `original_validity_contract.py`, `validator_counterexamples.py`, `VALIDATOR_COUNTEREXAMPLES.json`: exact-source gate failures.
- Copied `reduced.py`, `swe2d.py` and versioned reports for read-only review.

No new mathematical priority, award, journal outcome or field superiority is asserted.
