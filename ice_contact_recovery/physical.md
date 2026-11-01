# Physical candidate and a necessary thermal-clock null

Date: 1 November 2026. Independent research note; no active solver or result was modified.

## Decision

The strongest next physical test is **whether an oriented discharge–temperature cycle produces a consequential change in the survival of glacier-cliff water contact, even when the joint distribution of discharge and temperature and their total advected heat are identical**. The target is an order-sensitive geometry response, not another demonstration that peak discharge matters. This is a candidate discovery, not established originality or a proven environmental result.

Before testing it, enforce the exact thermal-clock null below. If changing temperature timing at constant discharge changes a passive, quasi-steady, fixed-debris result after matched thermal dose, the implementation or stated assumptions require explanation. That null can remove large redundant temperature campaigns and exposes numerical pseudo-memory.

## 1. Thermal-clock reduction

Let x collect ice elevations and any geometry-derived variables. Hold discharge Q, deposited material, boundary geometry, transfer coefficients, water properties and all nonthermal forcing fixed. At every x solve the quasi-steady hydraulic problem. Assume the steady temperature equation is linear and homogeneous in temperature relative to the melting point: advection–diffusion with zero-temperature ice and a boundary inlet superheat T_in(t). The steady temperature field and every ice heat flux scale linearly with T_in. Consequently

    dx/dt = T_in(t) F_Q(x).

For nonnegative measurable T_in, define thermal dose theta(t)=integral_0^t T_in(s) ds. If the ODE is uniquely solvable, x(t)=X(theta(t)), where dX/dtheta=F_Q(X). This follows directly by the chain rule and uniqueness; zero-temperature intervals simply pause evolution. Any contact criterion C(x)=0 therefore has a fixed dose theta_star, and the calendar forecast is the first time theta(t) reaches theta_star. Cumulative melt, being a geometry functional, is equal at matched dose. This is classical autonomous-flow time reparameterization, not a new mathematical theorem.

Conditions that break this reduction must be declared: thermal storage and finite advection time, variable Q, dissipation independent of inlet temperature, atmospheric fluxes, cold ice, temperature-dependent properties, debris motion/supply on physical time, changing external geometry, and unequal inlet temperature histories between branches. Finite quasi-steady relaxation errors and morphology-update errors can also break it numerically. Heat transfer varying with geometry is allowed: it is contained in F_Q(x).

Practical test: compare constant inlet superheat with two nonnegative pulse schedules at the same Q and final theta, using event-aligned steps and progressively smaller morphology steps. Compare the entire geometry at matched theta, not only a thresholded contact statistic. A residual that vanishes with refinement is numerical, not a physical memory effect.

## 2. Why changing discharge creates a meaningful extension

At two discharges Q_A and Q_B, define dose-domain vector fields F_A and F_B. They need not commute because discharge redistribution and heat depletion depend on evolving geometry differently at the two flows. Apply dose a in state A and dose b in state B. With the convention AB means first A then B,

    x_AB - x_BA = ab [D F_B(x0) F_A(x0) - D F_A(x0) F_B(x0)]
                   + O(a^2 b + a b^2).

This is obtained by expanding the two solution maps to second order. It is a classical flow-commutator calculation. A nonzero bracket is a local, physics-derived measure of ordering sensitivity; it is not proof of a novel environmental phenomenon.

For a scalar gap g with dg/dtheta=R_Q(g), the coefficient reduces to

    R_B'(g) R_A(g) - R_A'(g) R_B(g).

If R_B is a constant multiple of R_A throughout the relevant range, the bracket vanishes and discharge effects are reducible to one effective dose. If it does not vanish, the response requires at least the ordering of discharge-conditioned thermal doses. This distinguishes genuine state dependence from merely replacing T by Q*T or changing the input heat budget.

## 3. Decisive experiment that avoids a trivial forcing comparison

Use four strictly positive boundary states (Q_low,T_low), (Q_high,T_low), (Q_high,T_high), (Q_low,T_high), each held for the same physical duration. Run that cyclic order and its reverse, starting at the same phase/state, with identical initial geometry. Both protocols have the same joint Q,T histogram, total water and total inlet sensible heat integral rho*c*Q*T dt. They differ only in temporal ordering. All boundary transitions must be handled identically and exactly.

Start with every route wet and sufficiently far from thresholds, fixed debris and no external material supply. This measures a smooth commutator effect without wetting/drying or event-detection artifacts. Then test one configuration near an independently frozen contact boundary, where a small cumulative geometry difference could change a useful forecast. A route that differs solely because the two protocols finish in different hydraulic states is not evidence: return both to the same ordinary Q,T, relax to the same declared convergence criterion, and only then compare sustained contact.

For small block durations, the antisymmetric geometry difference should scale quadratically with duration and agree in sign/magnitude with an independently computed bracket. Repeat with frozen geometry: the geometric contribution must disappear. Repeat at constant Q and matched thermal dose: the thermal-clock null must hold. A two-state AB versus BA comparison alone may mostly reflect start/end transients; the four-state loop and ordinary-flow return control are preferable.

Do not infer a secular ratchet from one cycle. Reversed periodic protocols, phase-shifted protocols and finite-time transients are distinct comparisons. Demonstrate any accumulated effect across several cycles and subtract differences already present with fixed geometry.

## 4. Environmental interpretation and limits

The useful possible output is a forecast of when cliff-side water access is lost or restored under subdaily meteorological/runoff histories, at fixed seasonal water and heat totals. The candidate is consequential only if the ordering effect exceeds discretization, hydraulic partition, geometry and heat-transfer uncertainty, and improves predictions against an otherwise strong matched model that already resolves subdaily forcing. Beating a constant-average input alone is insufficient: earlier project results already showed that simple 3 h forcing was very competitive.

This should initially be a small mechanism test after the current hydraulic gate, not another large production campaign. If the bracket is negligible over admissible geometries, or a 3 h time-resolved reduced predictor already captures it, report the simplification and retain the thermal-clock acceleration. No extra model complexity is justified merely to manufacture a headline.

## 5. Primary-source novelty constraints

1. Petersen, Hock and Loso (2024), *Stream hydrology controls on ice cliff evolution and survival on debris-covered glaciers*, Earth Surface Dynamics 12, 727–753. Full open article inspected for stream avulsion, reburial and discussion. It already observes cliff abandonment after avulsions and reburial, and explicitly discusses exceptional-flow remobilization. Thus abandonment, debris removal and stream maintenance cannot be novelty claims. The new candidate would require a quantitative equal-joint-forcing order effect and predictive improvement. https://esurf.copernicus.org/articles/12/727/2024/
2. Karlstrom et al. (2013), *Meander formation in supraglacial streams*, JGR Earth Surface, DOI 10.1002/jgrf.20135. Primary article sections retrieved in search include discussion of strong diurnal discharge and its phase relative to solar forcing. Therefore the existence of a diurnal phase lag is already known. Full equations were not independently re-audited in this pass. https://doi.org/10.1002/jgrf.20135
3. Ogier et al. (2021), *Drainage of an ice-dammed lake through a supraglacial stream: hydraulics and thermodynamics*, The Cryosphere 15, 5133–5150. Open article retrieved. It reports time-variable discharge, temperature and incision including daily melt imprints. This is a strong environmental precedent and potential forcing source, not a validation of our cliff mechanism. https://tc.copernicus.org/articles/15/5133/2021/
4. Schuler and Fischer (2009), DOI 10.1029/2008JF001238, subglacial channel transit modeling. Primary article retrieved portions show time-dependent channel area and discharge/head phase offsets. This blocks broad claims that evolving ice channels first introduce phase memory here; its pressure-driven subglacial, near-melting-temperature model differs from the proposed open-channel cliff-contact outcome. https://doi.org/10.1029/2008JF001238

The targeted searches did not identify the exact four-state, equal-joint-forcing, ice-cliff-contact test. This is limited evidence of an opening, not historical priority. Geometric phase, Lie brackets, thermal dose and periodic morphodynamic forcing are established concepts. The potential contribution is an environmentally useful closure or consequential physical response, not ownership of this mathematics.

## 6. Executed independent falsifier

`thermal_clock_check.py` is a standalone synthetic two-route ideal-control/plug-flow model. Both floors evolve, flow is repartitioned by a common stage at every RHS evaluation, and one route has a finite stock-dependent insulating resistance. Four 48 h schedules (constant, early heat, late heat and four pulses) all deliver 86,400 K s at fixed 0.07 m³/s. With fixed stock, final geometry agrees to 2.8e-17 m across these schedules, as predicted. This numerical agreement is not an environmental accuracy estimate.

Adding an explicitly assumed stock-decay clock ds/dt=-s/(1 day), while keeping the same forcing totals and identical final stock, breaks the reduction: early versus late heat differs in toe-floor lowering by 0.01253 m. This is an illustrative counterexample, not a calibrated debris-export law or a new discovery. The important logic is that an independently evolving state makes timing physically relevant; temperature timing alone does not do so under the null assumptions.

Run `python3 thermal_clock_check.py --output thermal_clock_results.json`. The script fails if fixed-stock dose collapse is lost, if the extra-state effect disappears, or if final stock depends spuriously on temperature scheduling. The JSON contains all schedules' geometry, dose and solver work. No OpenFOAM, active project source, or external jobs are touched.

## 7. Executed commutator check: mechanism works, baseline effect negligible

`discharge_order_check.py` and `discharge_order_results.json` implement the smooth two-state AB/BA diagnostic without any active solver. The primary baseline retains h=200 W m^-2 K^-1, 0.5 K inlet superheat, ordinary Q=0.07 m³/s and the earlier two-route initial geometry. The two states are Q=0.035 and 0.105 m³/s. Four block durations, 600/300/150/75 s, are compared against an independently central-differenced flow bracket at three perturbation sizes. Every route stays wet, with route fractions above 0.45, excluding wet/dry events as the explanation.

The largest-block endpoint difference is only 3.29e-10 m in geometry norm, or 7.81e-7 of cumulative incision (0.0000781%). The common-ordinary-flow toe-share difference is 1.14e-9, and remains essentially that size after 600 s of common ordinary forcing. The bracket predicts the difference within 0.277% at 600 s and improves with halving. Thus the mathematical ordering mechanism is implemented correctly but is **environmentally negligible in this smooth baseline**. This is a no-go for a large smooth-ordering campaign or a headline based on this configuration.

An explicitly secondary lower-Q sensitivity case uses 0.0002/0.003 m³/s and equal initial hydraulic bed levels, retaining the same transfer coefficient. Its largest-block geometry difference is about 1.18e-6 m, and ordinary-flow toe-share difference is 7.53e-5. Observed halving orders are 1.970/1.985/1.992, and the bracket's relative error falls from 4.67% to 0.596%. These are useful numerical checks of the expected second-order mechanism, not justification to tune environmental parameters until an impressive effect appears. Its lower flows and geometry are not the main project baseline.

The result supports prioritizing independently verified contact transitions and real debris/geometry states over adding abstract forcing-order complexity. The null and bracket are valuable verification tools and might compress computation, but no consequential new environmental discovery was obtained in this small test.
