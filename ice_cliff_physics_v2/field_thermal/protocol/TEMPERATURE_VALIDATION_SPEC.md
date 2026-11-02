# Frozen field-temperature component validation

This specification is written before model fitting/scoring. Data source ranges, completeness, sensor counts and metadata have been inspected. It is not an externally registered blind study. It tests effective normal heat transfer between measured thermistors, not full cliff physics or wet-water Robin heat exchange.

## Inputs and eligibility

Use all12 named Kennicott profiles KDT_2020a/b,2021a/b/c/d,2022a/b/c,2023a/b/c. Verify original CSV and depth metadata hashes. Source cadence is5min, no missing readings in the audited files, and the source UTC/ADT columns differ by8h. Use the explicit UTC column; do not infer local time from a filename.

Boundary probes must be real sensors. The author source identifies `Base` at22cm in KDT_2020b as an imposed0°C column, not a thermometer. Its two real sensors are5/15cm, leaving no strictly interior probe for an observed-boundary subcolumn test. Exclude it from this PRIMARY test with that exact reason; it can enter a separately labelled assumed-basal-boundary analysis. The remaining11 profiles have at least three real depth sensors. Never treat the published depth-averaged conductivity values as independent material priors: their provenance can include the same seasonal ablation targets.

For each eligible profile choose the shallowest and deepest real sensor as the two Dirichlet boundaries. Predict every real sensor strictly between them. Do not impose an unobserved basal ice temperature in this primary test. Use recorded sensor depths and preserve any installation-change or medium ambiguity flag. Fixed depth over the test is a declared assumption, not verified absence of sensor motion.

Negative temperatures are retained. The linear effective-diffusion model can transfer signed temperature, but does not resolve pore-water freezing, liquid advection or moisture-dependent material properties. Report below-zero strata and model failure rather than deleting them to improve scores.

## Fixed temporal split

Let t0 be the first source timestamp. All eligible profiles have at least28complete days. Use [t0,t0+14d) as the training window and [t0+14d,t0+28d] as test. The first24h of training supply initialization/history; primary calibration objectives use [t0+24h,t0+14d), identically for all fitted families. Leave data after28d unscored in this first campaign.

Initialize the PDE once at t0 from piecewise-linear interpolation through all contemporaneous measured probes, including interior probes. Propagate through training and test; do not reset from held-out interior temperatures at the split and do not assimilate them during test. Every baseline may use the same initial measured state. Interior readings inside the training objective may fit parameters but are not repeatedly injected as physical state.

Boundary temperatures are prescribed by piecewise-linear interpolation of their original5min observations. This is an observed-boundary retrospective hindcast, not a forecast of future boundary/weather inputs. Score at the source observations every12th record (one/hour anchored at t0); no target interpolation or temperature-dependent timestamp selection. Retain native inputs for integration. If a source differs from the audited continuous cadence or includes missing values, stop and version the eligibility audit; do not silently bridge it.

## Frozen predictors

1. **Effective diffusion:** dT/dt=alpha*d2T/dz2, one homogeneous alpha per profile, fitted ONLY to the training objective. Minimize mean across interior probes of their MSE. Search a fixed logarithmic alpha grid of61points from1e-8 to1e-5m2/s; optional bounded refinement between the winning point's adjacent grid points is allowed to1% relative parameter tolerance. Flag boundary solutions; do not expand bounds after seeing test errors. These are numerical engineering bounds, not measured material-property confidence intervals.
2. **Quasi-steady depth interpolation:** interpolate the simultaneous two boundary temperatures linearly with depth. No fit.
3. **Training-only lag filter:** for each interior probe evolve dT/dt=(Tlinear(z,t)-T)/tau, initialized from the same t0 measurement. Fit tau on the same training timestamps from {0} plus61logarithmic values60s–7d. Tau0 is the quasi-steady limit. No interior test observations are inputs.
4. **Distributed-lag ridge:** per interior probe, intercept plus both boundary temperatures at lags0,1,3,6,12,24h. Standardize features with training data only; drop exactly constant columns by a1e-12 standard-deviation rule. Choose one ridge penalty per profile shared by its probes from {1e-6,1e-3,1,1e3}, using chronological inner training [t0+24h,t0+10d) and inner validation [t0+10d,t0+14d). Refit on all training objective timestamps. If the linear system is unsupported, report that failure; do not replace it with a test-selected alternative.

All four families receive identical boundary information and the same scored target timestamps. Ridge/filter are stronger temperature-transfer baselines, not mere no-physics strawmen. Do not select the field winner by dropping a baseline or a profile.

## Numerical verification gate

Verify the diffusion solver against independent constant/linear boundary and transient diffusion limits before calibration. Resolve temperature predictions with mesh/time refinement so the remaining numerical change at scored probe locations is below0.02°C, or report a numerical failure. The0.02°C target is an engineering numerical tolerance, not a claim of physical accuracy. Do not use a perfect energy ledger as a substitute for temperature convergence.

## Metrics and interpretation

Primary: MAE in°C per interior probe, then equal-probe mean within each profile, then equal-profile mean across11profiles. Also retain per-probe/profile RMSE, signed bias, daily residual summaries and fractions within0.2,0.5,1°C. Report all profile-level wins/losses against each baseline. Thousands of hourly readings are not thousands of independent glaciers. No population confidence interval or significance claim is licensed by this split. Manufacturer accuracy±0.2°C is not a95% prediction interval; boundary and interior sensor errors both contribute.

Published conductivity, test-interior data and full-season ablation must not enter this temperature model's calibration. All truth/prediction rows, parameter values, source hashes, initial states and split membership are saved. A result supports only effective thermal transfer in the measured subcolumns; it does not identify conductivity or flux magnitude, validate a water-film coefficient, demonstrate melt-volume skill, or establish Himalayan transfer.
