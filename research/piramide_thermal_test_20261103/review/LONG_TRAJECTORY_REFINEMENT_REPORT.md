# Source-only long-trajectory numerical sensitivity

The unchanged primary proxy and matched derivative-OLS proxy were run from
Dec2 toFeb10, preserving diffusivity, capacity, geometry, forcing reconstruction
and initial profiles. There were no state resets at test boundaries. No camera
outcome or empirical gain value was imported. Two helper tests pass.

Record:`model/frozen_piramide_future_batch_final/LONG_TRAJECTORY_NUMERICAL_REFINEMENT.json`.
Reproduction script:`model/future_batch_refinement.py`.

The frozen baseline is32debris/64ice cells,300s,2m ice. Every percentage below is
the change of the normalized dose for that exact future period. Because version3
gains remain fixed, it is also the relative change of its point prediction;
the gain need not be read to calculate this diagnostic.

| Model and refined grid | Dec30→Jan13 | Jan13→Jan27 | Jan27→Feb10 |
|---|---:|---:|---:|
| Primary64/128cells,150s,2m |+0.007968%|+0.002355%|+0.000799%|
| Primary64/256cells,150s,4m |−0.003360%|−0.019149%|−0.008934%|
| Derivative64/128cells,150s,2m |+0.012341%|−0.000601%|−0.000220%|
| Derivative64/256cells,150s,4m |+0.011569%|−0.002887%|−0.001143%|

The largest absolute prediction shift in these tests is0.0192%. This is a
finite mesh/time/deep-boundary sensitivity result, not a certified numerical
error bound or a statistical confidence interval. It does not cover geometry,
thermal-property, camera, initial-state or model-form uncertainty.

Complementarity products remain zero and no negative melt substep occurs.
Relative energy residuals are at most3.62e−13. At the2m refined domain the
primary signed deep-boundary heat is−122140J/m²; at4m it is−48530J/m². The
primary interface reaches−1.819C. These cold-content responses are supported by
the explicit ice energy balance; wet-debris pore freezing is still absent.

## Geometry and omitted advection

The raw unscaled primary component corresponds to1.809m cumulative ice-equivalent
retreat over70days, roughly90% of the2m ice-domain depth. The derivative component
corresponds to0.504m. Neither displacement was applied to the mesh. Source debris
and sensors may move with the melting surface, so these totals do not by
themselves establish how actual debris thickness or the sensor-to-ice gap changes.
They do require an explicit material-frame interpretation.

Using raw mean retreat velocity and the ice diurnal diffusion skin depth gives
mean Péclet scale0.0475 for the primary and0.0132 for the derivative component.
Those are scale diagnostics only: they are not bounds on transient advective
error and do not validate the omitted ice-advection term. Empirically corrected
outputs have a separate scaling and cannot inherit the raw physical ledger.

No source parameter, primary model, future interval, gain or study criterion was
changed in response to these diagnostics. Future field-error analysis remains
the independently sealed parent task.
