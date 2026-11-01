# Thermal memory and flood reactivation of debris-interrupted glacier streams

## Specific candidate contribution

Predict when a stream can resume supplying an ice-cliff toe after a debris disturbance. The important candidate consequence is that **exporting the obstruction need not restore ordinary-flow access**: the bypass can thermally incise while the toe route is insulated or inactive, leaving a geometric barrier after the debris has gone.

This is a quantitative reactivation problem, not a claim to discover debris insulation, avulsion, differential ablation, channel competition, or stream maintenance of ice cliffs. Those mechanisms have substantial prior literature. The targeted primary-source search did not locate the exact proposed finite-disturbance, complete-clearing, thermally evolving reactivation comparison. This does not prove historical priority.

## Work actually completed

- Implemented conservative solid-stock accounting for three grain classes, a two-route shared-head hydraulic control, advective water-heat depletion, ice-bed lowering, and explicit grain export.
- Derived a finite-flood exclusion test that can reject reconnection even after perfect instantaneous debris clearing.
- Ran exploratory parameter scans, matched flood-timing experiments, causal controls, time refinement, and an independently implemented adaptive reference solver.
- Replaced a misleading uniform insulating layer with a grain-footprint heat-transfer model; retained the old result as a comparison.
- Generated executable forecasting code and reproducible figures. Claude's solver, branch, worktree and jobs were not touched.

## Supported numerical result

In the **patchy-coverage reduced model**, a finite dense coarse deposit loses ordinary stream contact after approximately 8.34 days without a flood. An identical six-hour, eightfold-discharge flood on day 1 maintains contact through the 45-day horizon. On day 14 it removes essentially all retained material but does not restore ordinary flow after the flood.

The corresponding local toe melt is about 1.256 versus 0.1333 m³. This is a difference in synthetic local outcomes, **not a ninefold accuracy improvement or a whole-glacier melt result**.

The independently implemented Radau/Brent reference gives no-flood contact loss at day 8.340573 and local melt 0.1262706 m³; the 112.5-second production step gives day 8.339844 and 0.1262444 m³. For the late eightfold flood the reference gives 0.1333466 m³ versus 0.1332996 m³. Both predict loss of ordinary-flow access.

The sparse-deposit counterexample matters: when only 20% of initial solid volume is in the retained coarse class, the footprint model retains contact for the full 45-day horizon, while a uniform insulating layer incorrectly predicts loss near day 10.14. The sparse uniform-layer result must not be used as environmental evidence.

These configurations were selected and explored during development. None is a held-out, preregistered environmental validation experiment.

## Conservation and moving geometry

Let S_j be each solid grain volume, phi bed porosity, W and L deposition-patch dimensions, and b_t and b_b the toe and bypass ice-floor elevations. The model uses

    d = sum(S_j)/[(1-phi)WL],
    z_t = b_t + d,
    dS_j/dt = source_j - export_j,
    db_t/dt = -P_t/(rho_i Lf WL).

The shared stage eta solves

    Q = K_t [eta-z_t]_+^(3/2) + K_b [eta-b_b]_+^(3/2).

For ideal rectangular critical controls, K=(2/3)^(3/2)*sqrt(g)*W, critical flow depth is 2/3 of head, and velocity satisfies u²=g*depth. Head must not be confused with the critical-section depth.

Each branch removes sensible heat consistently with its modeled conductance:

    P = rho_w cp Q T * [1-exp(-h_eff WL/(rho_w cp Q))].

For patchy cover, cube projected footprint is S_j/D_j. A parallel bare/covered heat-transfer model uses those footprints and a grain thickness resistance D_j/k. This is a specified subgrid closure, not a resolved grain-scale conjugate heat solution. Mean toe crest height is also a closure; a clear corridor between rocks may sustain flow that a mean-height calculation misses.

The routing approximation neglects the small additional cold meltwater discharge. Heat and solid ledgers close within the declared equations; this is not a claim of a fully coupled variable-discharge water-mass solution.

## Finite-flood exclusion proof

Let the ordinary return discharge be Q0. If a bare toe crest remains above

    b_b + (Q0/K_b)^(2/3),

it receives no ordinary flow. Suppose all admissible processes can lower the toe ice by no more than B_max during the flood and the bypass does not rise. Even after removing all sediment, the final toe–bypass gap is at least

    g_min = b_t(start) - b_b(start) - B_max.

Therefore

    g_min >= (Q0/K_b)^(2/3)

excludes reactivation of that route under ordinary return flow. The proof is the hydraulic activation condition plus a bound on possible moving-boundary displacement; it does not assume a hysteresis switch.

The executable estimate bounds water-driven B_max using both the available upstream sensible heat and the bare-water conductance ceiling, then admits additional delivered heat explicitly. Initial grain sensible heat, atmospheric/lateral heat and other mechanisms must be included in that additional bound if relevant. Ignoring actual heat sources would invalidate the exclusion claim.

For a dry toe with fixed stock, fixed properties, and constant upstream flow, continuing bypass incision gives

    Q_rewet(t) = K_b [g(start) + integral(m_b - m_t,dry) dt]^(3/2).

This is a conditional evolving activation threshold. It is not an assertion of permanent abandonment: sufficiently large floods, new flow paths, widening, collapse, changing runoff or a rising bypass can change the result.

## What would make this an environmental contribution

The intended improvement is a physically computed stream-access history for a melt hotspot. Static proximity or successful sediment export is insufficient to specify that history in these equations. A useful predictive product returns contact loss, temporary rewetting, post-flood reconnection, local ice loss and uncertainty from forcing and geometry.

The decisive remaining physical test is spatial. Resolve alternate water paths around an explicit coarse deposit and evolving ice. If an uncovered corridor incises and keeps feeding the cliff, the mean-crest reduced model may overstate disconnection. A contribution must survive that test rather than declare the reduced model to be truth.

The study also needs finite-width evolution, sliding versus rolling/lift, thermal storage, external heat, deposit geometry and turbulence sensitivity. These are explicit falsifiers, not reasons to invent fitted feedbacks. Computational validation against independently resolved physics can support a computational mechanics paper without claiming a calibrated Himalayan operational forecast.

## Closest prior work

- Petersen et al. (2024), https://doi.org/10.5194/esurf-12-727-2024: stream maintenance, large-clast blockage, avulsion and reburial already established observationally/conceptually.
- Buri et al. (2016), https://doi.org/10.1002/2016JF004039: dynamic three-dimensional cliff geometry, ponds and geometric debris reburial.
- Evatt et al., https://doi.org/10.1017/aog.2024.34: quasi-steady cliff geometry and debris-thickness thresholds.
- Bolla Pittaluga et al. (2003), https://doi.org/10.1029/2001WR001112: sediment-driven bifurcation instability and branch closure.
- Douglas and Lamb (2024), https://doi.org/10.1029/2023JF007452: thaw/erosion moving fronts and transport-limited buffering.
- Beaud et al. (2018), https://doi.org/10.1029/2018JF004779: coupled sediment transport and melting ice-walled subglacial channels.

See the primary-source notes in prior_art.md and PULSE_PRIORITY_CHECK.md. The potential novelty is the specific predictive reactivation problem and demonstrated additional physical consequence, not ownership of the component laws.

## Executable event-time predictor

After mobile fines clear, a fixed retained stock gives a constant mean sediment height d and patch conductance. Let g=b_t+d-b_b. The hydraulic partition is uniquely determined by g and Q. Define

    R(g) = melt_bypass(Q_b(g)) - melt_toe(Q_t(g)).

The toe share decreases and the bypass share increases with g; each branch's heat uptake increases with its flow. Thus R increases with g under the fixed-parameter closure. If R(g0)>0, the dry threshold g*= (Q/K_b)^(2/3) is reached at

    T_loss = integral[g0 to g*] dg/R(g).

The monotonicity also gives

    (g*-g0)/R(g*) <= T_loss <= (g*-g0)/R(g0).

This eliminates time marching for that conditional forecast. It is implemented in `compute/fixed_stock_forecast.py`. The scalar quadrature itself is standard mathematics; the proposed contribution is the physically specified post-disturbance access prediction, not a claim to invent separation of variables.

In the dense-deposit development case, the predictor gives 8.341198 days; the independently integrated full reduced system including the short initial fine-clearing phase gives 8.340573 days. No time coefficient was fitted. For the sparse case, it predicts loss after about 68.69 days, correctly outside the 45-day simulated horizon. These comparisons verify a reduction within the specified equations; they do not establish field forecast accuracy.
