# Specific exponential-moisture model input contract

This supplements the feature audit for the selected statistical predictor. No Tower2camera outcomes have been opened.

## Canonical inputs

`Tower1_CANONICAL_FEATURES.csv` and `Tower2_CANONICAL_FEATURES.csv` each contain3361native half-hour labels from2December03:00 through10February03:00inclusive. They provide rawSM3, interpolated native2m air temperature(TA_3), source air brackets/gaps, nominal moisture-probe-depth temperature references and operating-range QCflags. Native air streams are supplied separately. The builder requires finite bracketing air observations with gaps<=18h and forbids extrapolation. It does not alter camera data or fill missing moisture.

`CANONICAL_FEATURE_MANIFEST.json` binds source hashes, feature hashes, training-reference means and site-specific epochs. Tower1development uses30December00:00, while Tower2primary uses30December03:00. Their training means are0.0501925595 and0.0704194940respectively, computed on the exact half-open training block. Source clock labels remain naive local; hourly/10minute observation-label semantics are an explicit input assumption.

The soil-temperature reference is nominalDT5at0.50m forTower1, and linear interpolation betweenDT4at0.33m andDT5at0.44m to0.39m forTower2. These are temperature QCproxies, not an independent thermometer inside each moisture probe. Across the full3361label window, Tower1has189nominal below-zero reference values; Tower2hasnone. With±0.4°C temperature uncertainty,1337/Tower1and894/Tower2samples cannot be guaranteed above0°C. Do not erase these flags or call the moisture series independently certified within manufacturer calibration conditions.

## Exact invariance and remaining uncertainty

For fixed beta, define D_B=integral_B max(Tair,0)*exp[beta*(theta−mu_train)]dt. A local training-only gain yields predicted melt Ytrain*D_test/D_train.

If a sensor has a constant additive biasb across its own training and evaluation observations, both theta andmu_train shift byb; the centered feature is unchanged. Thus predictions and the complete Tower1beta-fitting objective are exactly invariant. Different constantbvalues atTower1andTower2also cancel separately. This is a useful engineering property of the chosen form; it is elementary algebra, not a claim of new physics.

A multiplicative calibration factor does not cancel: a fixed transferred beta multiplies a different scale of anomalies. Time-varying offsets, drift, soil-temperature effects, sensing-depth changes and ice/material contact changes also remain. If raw-reading perturbations were separately bounded bydeltaandbetawerefixed, the dose-ratio error factor is bounded by exp(±2*abs(beta)*delta). Centering cancels from the ratio, so no extra independent mean-reference factor is needed. Manufacturer±.031is typical mineral-soil accuracy, not a guaranteed rock-debris bound; report this as a conditional engineering sensitivity. Refittingbetaunder uncertain readings introduces further uncertainty.

Three synthetic tests verify additive invariance, the attainable fixed-beta drift bound and non-invariance to multiplicative calibration. They use no field target.

## Frozen experiment

Primary Tower2training:2December03:00→16December03:00. Primary evaluation:16December→30December→13January→27January→10February, always03:00. Freeze all four periods and the observation operator before any Tower2target release. Camera units/geometry are source-bound from the existing sealed loader, using exposed-stake-distance increments×917to mmwater equivalent.

Beta is selected from already inspected Tower1development outcomes, using its correctly bound four intervals, then locked before Tower2camera release. Each applicable model gets exactly one Tower2training-ablation amplitude gain. Retain air-only and the equal-budget ratio-PDDpower competitor. Any ratio feature treatment for negative Tower1moisture must be declared before fitting; a fixed positive feature offset is a statistical transform, not a corrected water-content measurement. Do not discard periods to make this baseline undefined or favorable to the candidate.

Publish every normalized model prediction, all code/source hashes and selection rules before training release; publish all issued predictions before evaluation release. Keep all model results, timing conventions and sensitivity cases. The exponential model is an empirical instrument-conditioned predictor; it does not conserve a complete ice/debris energy budget or identify an intrinsic conductivity law. Its originality and environmental advantage require evidence beyond this invariance property.
