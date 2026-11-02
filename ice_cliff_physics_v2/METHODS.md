# Governing contracts and checks

## Finite stream heat

For liquid heat-capacity rate F=rho_w cp_w Q and simultaneous bare/debris conductances Gb,Gd, define N=(Gb+Gd)/F, phi=(1-exp(-N))/N and a=Gd/(Gb+Gd). With a representative debris surface temperature Ts, the mean contact-water temperature is phi Tin+(1-phi)a Ts. Bare and debris exchange use that same water field. Signed debris exchange is retained: warm debris may export heat into the stream.

The upstream floor removes heat first. The face receives the remaining sensible-water power, not the original inlet power a second time. Zero-flow exchange is zero; standing-water natural convection is not inferred from this limit.

## Normal debris diffusion

For normal coordinate xi in [0,1], areal heat capacity C_A=rhoc d and conductance K_A=k/d define

    C_A dT/dt = K_A d²T/dxi².

The ice boundary is0°C. The outer boundary couples finite-water exchange and the prescribed net surface source. Cell-centred finite volumes and half-cell boundary resistances yield a tridiagonal implicit step. The same discrete boundary flux enters the water, debris storage and ice ledgers. Time stepping is first order; analytic and separate space/time tests establish the numerical orders.

The effective water Robin condition admits an independent eigenfunction reference. Eigenvalues solve lambda*cos(lambda)+B*sin(lambda)=0. Exact integrated-series sum rules avoid misleading early-time cancellation. See `thermal/ROBIN_REFERENCE.md`; no new mathematical priority is asserted for this classical boundary problem.

Each wet/dry pool has its own normal profile. Newly exposed material imports declared enthalpy; class conversion carries the donor profile. Canonical material area prevents a one-ULP sum discrepancy from being mistaken for real material loss. True contraction without an explicit outgoing-material operator remains refused.

## Geometry

The extruded ice front is x(z), piecewise constant in height. Band recession from net latent energy E is dx=E/(rho_i L area_vertical). Face removal precedes floor lowering on the updated basal footprint. Their volumes are disjoint. A polygon-area comparison independently checks the exact geometric operator.

The hydraulic routing proxy uses the basal front position. Its uniform-step ice volume is generally NOT the remaining physical volume. Global melt accounting uses the complete band profile. Horizontal ledge areas are reported; current integrated runs prescribe zero ledge heat and make no collapse prediction.

Class boundaries are inserted losslessly before sampling heat. This preserves ice and heat while avoiding unresolved basal strips. Actual mesh growth is recorded, and the maximum-band limit refuses rather than silently coarsening.

## Time and observation semantics

Hydraulics, masks and boundary areas are sampled at `surface_partition_sample_time_s`. Geometry, port reservoirs and `debris_columns_post_geometry` belong to `state_sample_time_s`. The intermediate thermal result is pre-geometry/remapping. These are distinct states, not duplicate sensor observations.

Geometric undercut, mass removed and predicted water contact are separate outputs. A fixed-z covered footprint is a prescribed Eulerian mask, not an inferred moving debris stock. Arbitrary initial profiles require explicit coordinate/reference-frame provenance.

## Conservation and validity

The global accounting is incoming water heat + prescribed atmospheric heat + signed material enthalpy transfer = outgoing water heat + thermal storage change + fresh-ice warming + latent heat. Geometry loss must match phase volume. The maximum time-series normalized residual is reported, not just an end-of-run signed cancellation.

Route-specific phase-water ratios are checked. An explicit separated downstream free-drainage contract can qualify an otherwise unsupported one-way source, but does not constitute a solved meltwater-feedback flow field. The full shallow-water/ice-density jump and momentum-source problem is not silently declared solved.

Do not use retrospective observed forcing as an issued forecast, numerical error as a statistical confidence interval, or a component-temperature comparison as full glacier-melt validation.
