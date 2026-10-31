# A stronger, falsifiable physical contribution

Prepared as offline research support. This document does not change the active Claude milestone or authorize production runs before its verification gates pass.

## Candidate question

Can near-freezing freshwater buoyancy in slow recirculation near an ice-channel step materially change **where incision concentrates and how a hydraulic control lowers**, even when bulk-flow screening suggests forced convection dominates?

The possible finding is an interaction between free-surface hydraulics, local thermal stratification and evolving ice geometry, with an observable consequence for a declared forecast. It is not the discovery of density anomalies, melting, recirculation or nonmodal mathematics.

This is an unresolved candidate, not a proved new direction. A 2021 ice-morphodynamics paper already discusses destabilizing buoyancy above ice. A 2025 conference abstract explicitly describes supraglacial-channel calculations with the freshwater anomaly. A 2026 experimental paper studies near-field geometry/wake effects on local ice melting. These sources narrow the gap considerably. The specific step/control interaction and a transferable predictive benefit require full priority checking and direct evidence.

## Why it could matter

Water near 0 C can be lighter than slightly warmer fresh water. An ice bed can therefore supply destabilizing thermal stratification. A fast bulk stream may pass before buoyancy strongly acts, while water retained in a slow recirculation zone may behave differently. This offers a plausible route for local heat transfer to depart from a formula based only on bulk speed and temperature.

That explanation must compete with ordinary separation/reattachment, turbulent mixing, roughness, differing wetted area and numerical thermal diffusion. Do not attribute a discrepancy to buoyancy simply because an EOS option changes a run.

The environmental value would be improved predictions of control lowering or persistent melt concentration. A short-lived velocity change without measurable ice response is not enough. Integrate the heat-flux difference over the forecast horizon with physical latent heat and report whether it affects the chosen prediction beyond uncertainty.

## Decisive comparison after current verification

Start with fixed geometries to identify the thermal mechanism:

1. Smooth geometry, density-independent liquid thermal flow.
2. Step geometry, the same density-independent physics.
3. Smooth geometry, the justified freshwater buoyancy model.
4. Step geometry, the same buoyancy model.

Keep water and heat inputs, surface conditions, material parameters, boundary compatibility and fitting information controlled. Report area-integrated and local fluxes separately. The difference of differences isolates the modeled geometry–buoyancy interaction; it is not by itself evidence that the model is physically correct.

Compare a local negative-expansion approximation with the nonlinear EOS where appropriate. If the simpler approximation is adequate, retain it. Including a more elaborate EOS does not automatically make a model better.

For a Boussinesq thermal approximation, derive its relation to the air–water mixture, ice-melting density jump, pressure and energy conventions. Do not apply a freshwater thermal-density correction to the air phase, double-count hydrostatic gravity, or confuse a small liquid thermal density variation with the much larger solid/liquid phase-density difference. The new regime screen does not implement this coupling.

Only after fixed-geometry effects survive independent thermal/hydraulic checks should geometry evolve. Repeat a small subset in 3D where convective/turbulent structures require it; a 2D vortex is not evidence for a 3D environmental turbulence mechanism.

## Predicted discriminants

If local buoyancy is responsible, the result should depend on gravitational orientation and thermal density contrast in a way consistent with the equations. A diagnostic reversal of gravity-relative interface orientation can test causality, but the changed configuration is a numerical control, not a second glacier validation.

Assess a local residence/buoyancy timescale ratio, together with thermal diffusion, shear and turbulence. A ratio near one is a motivation for testing, not a universal transition threshold. Local residence time must be estimated from trajectories, transport or resolved turnover statistics; length divided by an arbitrarily chosen small speed is only a sensitivity calculation.

If the effect is simply a separation effect already captured by constant-density CFD, the thermal-density ablation should add little. If it is purely numerical, its magnitude/location may fail mesh, interface, timestep or scheme tests. If it is merely a larger wetted area, the area-normalized and integrated comparisons will expose that distinction.

## Mathematical contribution to prepare

Use the verified coupled model's tangent response to ask which *admissible* perturbations affect the declared incision observable most strongly. The new forcing tool can evaluate time-varying reduced operators with linearized matched-water/heat constraints. The previous response tool handles initial-condition perturbations.

Both require a physical operator. The synthetic examples included here establish numerical correctness only. Derive a structured heat-exchange perturbation that preserves water/ice power cancellation, normalize states and inputs physically, and validate the resulting response against independent finite-amplitude model runs.

A useful result might be a quantitative local regime criterion or reduced correction that predicts the mechanism's influence on incision in new configurations. It must be evaluated against a strong local/equilibrium closure and a fairly calibrated alternative. It is not automatically a novel theorem or a universal closure.

## Primary endpoints

Predeclare a physical control region or spatial bed-lowering norm before results. Report:

* Control lowering and bed-lowering profile at a fixed horizon, in metres.
* Integrated water-to-ice heat, total melt and downstream heat export.
* Location and magnitude of persistent incision concentration.
* Numerical and thermal-closure uncertainty, input/geometry uncertainty and cost.

If using lake release, establish a compatible storage coupling and a causal hydraulic path. A downstream hotspot cannot automatically influence an upstream lake. The primary local thermal result may still be useful without a lake-hazard claim.

## Falsification and stopping conditions

Drop this mechanism as a predictive addition if it is negligible under defensible environmental conditions, fails 3D/numerical checks, is explained by known geometry effects, has no forecast consequence, or cannot outperform an appropriate simpler closure on held-out configurations.

Preserve earlier negative findings from the network/event studies. The new contribution must come from demonstrated newly resolved physics, not a relabeling of the already-known fact that forcing order can matter.

The strongest final paper would connect an identified mechanism, a quantitative explanatory relation and a transferable prediction improvement. Publication level is determined by those results, not the number of equations or processes added.
