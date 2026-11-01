# Alternative contribution: an exposed-cliff sediment-budget certificate

## Environmental question

Can a model predict persistent ice-cliff retreat when the stream cannot export the debris released by that retreat? The output is a number with direct forecasting meaning: the maximum exposed retreat compatible with measured debris supply, admissible toe storage, and bounded export, or equivalently the minimum export required by a candidate melt forecast.

This is a conservation-based falsification and coupling component, not by itself a new discovery or a universal glacier forecast. It is more defensible than asserting that stronger melting automatically causes less total melting.

## Why this is a real model gap worth testing

Buri et al. (2016) explicitly model reburial through slope and debris-view-factor thresholds, and describe a source-area/debris-characteristic representation as a possible improvement requiring more information. Their limitations section is directly relevant: https://doi.org/10.1002/2016JF004039 . This does **not** imply they omit reburial.

Petersen et al. (2024) directly observe stream removal, undercutting, abandonment and subsequent reburial at Kennicott Glacier and offer a conceptual geomorphic account: https://doi.org/10.5194/esurf-12-727-2024 . The general stream-survival mechanism is therefore already published.

Kneib et al. (2022) identify changing cliff–hydrology connectivity as important to expansion/reburial and future energy-balance models: https://doi.org/10.5194/tc-16-4701-2022 . Again the broad coupling is not new.

A candidate contribution would be replacing a heuristic burial criterion with an explicitly conserved solid-volume source, finite toe capacity and bounded stream export, then discovering where this materially changes predicted seasonal cliff survival and volume loss. That consequence has to be tested; a literature gap alone is insufficient.

## Exact control-volume statement

Let x(t) be horizontal retreat of an exposed cliff of constant along-front width B and vertical height H. Let a=(1-porosity)*h_m be solid debris volume per unit horizontal mantle area, and c be solid debris volume per unit glacier-ice volume (if englacial debris is represented). For this idealized geometry the source per unit retreat is

    kappa = B*(a+c*H)     [m² of grains per m retreat].

Let S be toe solid-grain volume, I external grain input rate, and E exported solid-grain volume rate. Then

    dS/dt = kappa*dx/dt + I - E.

All debris fractions must use solid volume, never a mix of bulk and solid volume. If mantle debris has already been included in kappa, do not also prescribe it in I.

Integrating yields the exact identity

    kappa*x(t) = S(t)-S0 + cumulative_export - cumulative_external_input.

For variable geometry/source kappa(t), the left side becomes integral kappa(t)*x_dot(t) dt. A conservative lower bound kappa>=kappa_min gives a corresponding retreat bound. If kappa_min=0, there is no debris-budget retreat bound.

## Toe capacity from declared geometry

For a planar face with inclination theta and a pile at repose beta, 0<beta<theta<=pi/2, choose a pile crest at vertical height z on the ice face. The cross-sectional gap area is

    A(z) = z²/2 * (cot(beta)-cot(theta)).

Multiplying B and solid packing fraction phi gives

    S(z)=phi*B*z²/2*(cot(beta)-cot(theta)).

Full toe-to-crest coverage gives Scrit=S(H). A chosen lower exposure criterion z*=eta*H gives Scrit=S(z*), but eta is an explicit application definition, not a fitted universal constant. This capacity is a restrictive planar quasistatic pile construction. Uneven roughness, undercutting, lateral leakage, slope motion or collapse require a geometry-resolved capacity. Do not silently use this formula on an arbitrary cliff.

## Necessary survival inequality and dual forecasts

While a candidate exposed state requires S<=Scrit,

    kappa*x <= Scrit-S0 + cumulative_export - cumulative_external_input.

Thus a candidate retreat trajectory x_model requires at least

    required_export = max(0, kappa*x_model - (Scrit-S0) + cumulative_external_input).

If independently computed/observed maximum export is lower, the trajectory violates the accepted exposure-storage model. This is a **necessary condition**, not sufficient evidence of survival: slope collapse or hydraulic disconnection can still kill a cliff.

Under fixed geometry the thermal upper bound is

    x <= incoming_usable_latent_energy / (rho_i*L*B*H).

An energy- and debris-budget ceiling is the minimum of the two bounds. This does not mean a buried cliff stops all melting. It limits the exposed-front retreat before the specified burial criterion; residual sub-debris melt needs another branch.

## Consequential, falsifiable campaign

1. Prescribe stream Q(t), temperature and geometry; compute an upper export capacity using a closure checked for clast-to-depth ratio. If grains are too large for a bedload law, use force-balance admission plus explicit resolved removal, not an extrapolated deep-bed formula.
2. Use the same radiation, stream heat and initial geometry for the heuristic-reburial baseline and conserved-grain model.
3. Compare seasonal exposed area-time, first loss-of-exposure time, melt volume and required versus available grain export.
4. Independently resolve 2-D/3-D ice/thermal evolution plus individual grain removal for a small set bracketing the predicted capacity threshold.
5. Vary geometry/source and closure uncertainty; classify predictions as certified impossible, conditionally feasible, or unresolved. Do not call feasibility a prediction of persistence.
6. Seek a regime where all reasonable transport/geometry assumptions give the same change in survival classification and a large integrated-melt effect. If no such regime exists, the proposed environmental improvement fails.

## A tempting reversal that cannot be asserted

A model Sdot=kappa*M(S,T)-E with fixed clearing and monotonically increasing M in T often has a scalar cumulative-response comparison. Increased warming can saturate the melt response without reversing it. Selective coarse armoring, hydraulic disconnection, stream-path migration, finite cliff geometry or time-varying transport may create additional state, but each requires independently justified constitutive physics. Adding an arbitrary nonlinear inhibition function to manufacture hysteresis is not a discovery.

## Mathematical novelty assessment

The integral identity and bound are elementary consequences of conservation. They should not be sold as a new theorem comparable to a foundational fluid-dynamics result. The potential novelty is a verified model coupling and an independently demonstrated environmental regime missed by current heuristic parameterizations. That claim is still open.
