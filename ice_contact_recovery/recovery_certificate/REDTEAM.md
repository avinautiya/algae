# Mathematical red-team of `bounds.py`

Independent review, 1 November 2026. Reviewed `opening_bound`, `exclusion`, and the seven tests in `test_bounds.py`. No code modified by this reviewer. The author subsequently corrected the main issue and this note records the revised API.

## Verdict

The shortest-path result is correct for a **declared whole-cell graph opening problem** with nonnegative additive costs. The width relaxation is correct for a declared rowwise area-removal problem. Their maximum is a lower bound only for a feasible class satisfying **both** definitions. It is not automatically a lower bound for real finite-width water paths on a continuous ice surface.

**Correction reviewed:** the default now uses only the relaxed finite-width term. The full-cell path term enters the exclusion only with explicit `resolved_cell_route=True`. This fixes the narrow-trench false exclusion below. The revised default is valid for a piecewise-constant terrain model with transverse wetted width at least the declared W at every streamwise section, with a genuine global level upper bound and melting-only opening. Cell-centre DEM sampling alone still does not establish a lower envelope of real terrain. Neither mode is a new shortest-path theorem or a universal environmental non-recovery certificate.

## 1. Whole-cell versus continuum width — important counterexample

Take N rows of identical ice elevation z=1, level upper bound0, and dx=dy=1. There are no walls. The graph path term is N cubic metres. A continuous straight trench of width0.1m across the N-m domain needs only0.1N cubic metres of ice removed, before allowing any positive depth margin. Thus for an admissible physical route narrower than a cell, `max(path_volume,width_volume)` can overestimate the true minimum tenfold, and `exclusion` can falsely exclude reconnection.

In the original implementation, setting `minimum_cross_section_width_m=0.1` did not fix this because the path term still forced full-cell costs. The revised default returns0.1N for this example. Width0 returns a vacuous0 default bound, correctly avoiding finite-energy exclusion of arbitrarily thin continuum trenches.

Even widths greater than dy do not by themselves establish the graph lower bound when a continuous corridor can cut partial cells, exploit within-cell topography, or cross corners. A cell-centre sampled DEM is not a cellwise lower envelope of the real terrain. For rigorous continuum exclusion one needs either a continuous geometric relaxation known to underestimate every admissible route, or a proof mapping every physical admissible route into a whole-cell graph path without increasing removal cost. Such a proof is generally false for full-cell node costs.

**Safe application:** the forecast model explicitly allows only complete cell-column lowering and defines contact by a four-neighbor wet-cell path. The certificate then concerns that discrete model. Report grid dependence, and avoid interpreting a positive exclusion as resolution-independent environmental impossibility.

## 2. Proof within the discrete feasible class

Let S be the set of cells opened to the stated upper level by irreversible ice removal. Assume each selected cell is lowered over its entire area and has nonnegative volume cost c_v=A_v max(z_v-H,0). If S connects first to last row in the four-neighbor graph, it contains a simple source-target path P. Nonnegative costs give sum(S)c_v >= sum(P)c_v >= shortest_path_cost. Dijkstra is appropriate because all c_v are nonnegative; zero-cost cycles do not invalidate the argument.

If the same opened set also supplies at least widthW at every row, its rowwise removal cost is at least the fractional sorted-cell fill implemented by `width_volume`. This relaxation may distribute occupied pieces across disconnected segments and therefore cannot overestimate the cheapest admissible row cost. Summing disjoint row areas preserves the inequality. Consequently max(path,width) is valid for this intersection of constraints.

The tests correctly compare the path term against exhaustive3x3 graph subsets and compare one width example against connected whole-cell subsets. They do not test or prove a continuum lower bound, hydraulic sufficiency, or robustness to unresolved topography.

## 3. Water level and velocity head

If H is a **true upper bound on the local free-surface elevation everywhere along every candidate route and throughout the forecast**, positive-depth free-surface flow requires local bed elevation less thanH. The geometric lower bound then does not need a velocity-head correction: flow cannot have its bed above its own actual free surface.

However, an upstream reservoir stage, prescribed inlet depth or mean water level is not automatically such a bound. A moving jet, supercritical inflow, pressure recovery, hydraulic jump, transient run-up or backwater can create local water levels above the chosen inlet stage. In inviscid steady flow with inlet velocityU, the available total head is eta_in+U²/(2g). A bed crest between eta_in and that total head can be traversed without first melting below eta_in, depending on channel control. Using eta_in alone could therefore falsely exclude an admissible route.

A verified upper bound on total mechanical head can be a more conservative geometricH in suitable steady passive-flow settings, because eta<=total_head. But an unsteady or externally forced domain needs a separately justified head bound including boundary work, pressure and transient storage; Bernoulli along a steady streamline is not a universal proof. Do not use a boundedQ alone as a boundedstage without solving or bounding storage and controls.

## 4. Energy accounting assumptions

Multiplying removed volume by rho_i L gives a latent requirement for melting initially temperate ice, assuming the route is opened only by local melting. Cold ice requires more energy, so omitting its sensible warming is optimistic. Perfect heat allocation and neglect of conductive/advective losses are also optimistic.

The upper budget passed to `exclusion` must be a genuine upper bound on **all heat that can reach the candidate ice** over the relevant physical time: incoming sensible heat, atmospheric/radiative input, dissipation and any debris or liquid thermal storage release. A water sensible-heat budget alone is not necessarily an upper bound. Do not integrate the hydraulic spin-up clock as physical elapsed time.

Mechanical fracture, ice-block detachment, lifting or collapse can open a route without melting its full latent volume. Such processes invalidate this particular necessary condition unless explicitly excluded or their possible removals are accounted for beforehand. Sediment removal that uncovers a pre-existing channel must be applied before defining ice elevation; debris elevation is not removable ice volume.

## 5. Walls, connectivity and targets

An infinite result from a blocked transverse row is valid only if those cells truly remain inaccessible to every admissible opening process. A porous or mobile debris deposit cannot be designated an immutable wall merely because it is initially dry. Lateral bypass outside the computational boundary, under-rock gaps and englacial routes are excluded by the graph definition, not by nature.

Four-neighbor connectivity prevents corner-only links and is a sensible finite-volume convention. It is not automatically the topology of a reconstructed continuous bed. The first/last-row target also represents geometric through-connection, not necessarily ice-cliff contact or sustained post-flood heat delivery. A more application-specific source and toe-target mask would need their own definition.

## 6. Implementation edge cases and reporting

- Validate the provenance and units of the bound supplied to `exclusion`; currently any dictionary with the expected key is trusted. NaN required energy would silently yield `False`, so imported malformed records must be gated elsewhere.
- Boolean material/geometry inputs are currently accepted as finite scalars by Python; this is an API-quality issue, not a flaw in the theorem.
- The1e-12m width feasibility slack can declare a minuscule over-wide request feasible. This is a numerical tolerance that should be documented or made relative to the geometry scale.
- The optimistic path returned minimizes path-only cost. If the width term dominates, that path is not a constructive minimizer of the width-constrained problem. Label it accordingly.
- A zero bound means no exclusion, not positive depth, positive discharge, persistence or recoverability.
- Comparing energy strictly below the mathematical bound is correct. Numerical applications need a declared error margin; nearly equal floating-point budgets cannot support a categorical physical conclusion.

## Handoff priority

1. Keep the revised distinction explicit: default fractional-width relaxation versus opt-in whole-cell-route requirement. The counterexample regression should remain in the test suite.
2. Require or document the origin of the head/level upper bound.
3. Carry omitted opening mechanisms and heat sources as explicit exclusions.
4. Treat the result as a useful screening/diagnostic tool. Validate conditional exclusions against a refined spatial model; do not use the same graph criterion on both sides as independent evidence.

## 7. Physical-cycle proposal originality check

Reviewed `../physical.md`. The thermal-dose null is rigorous under its stated autonomous quasi-steady, homogeneous temperature-response assumptions. Its usefulness is verification and campaign reduction; time reparameterization is classical. The AB/BA bracket expression is the standard second-order flow-composition expansion, with the displayed sign consistent with firstA thenB.

The four-state rectangle is a stronger matched-input experiment than comparing averaged versus pulsed forcing: the joint Q,T histogram is identical, not merely separate averages. However, reversing the loop still generally changes startup history and the state immediately before ordinary-flow return. Returning both to the same ordinary boundary state removes instantaneous hydraulic-state differences but does not prove a geometric phase invariant. A single pair proves finite-time order dependence only; repeated cycles, phase rotations and asymptotic small-duration scaling distinguish robustness from a selected transient.

For this physics, thermal erosion, diurnal phase effects, evolving-channel memory and sediment hysteresis already have direct glacier precedents: Karlstrom2013, Ogier2021, Schuler2009 and Delaney2025 (links in the two companion reports). Therefore a broad discovery of thermal morphodynamic memory would be indefensible. I did not locate an exact equal-joint-Q/T reversed-loop test of persistent stream–cliff contact in the focused searches, but that does not establish priority. The credible contribution is a quantitative closure predicting order-dependent contact or recovery better than a strong time-resolved geometry-state baseline.

There is a practical no-effect competitor: if the two discharge-conditioned vector fields are globally proportional, their flows commute and every matched effective dose collapses. A nonzero bracket may also remain far too small to change any useful environmental quantity. Near a wet/dry discontinuity the smooth bracket derivation need not apply; establish it away from thresholds before using separate nonsmooth event analysis near contact loss. Fixed debris and ideal inflow cycles isolate a mechanism, not a whole-glacier result.

This is a defensible bounded mechanism study, not sufficient evidence of Nature-level significance or historical mathematical novelty. It should not delay the hydraulic verification gate or trigger a large campaign before the small, independently computed bracket/null comparison is persuasive.
