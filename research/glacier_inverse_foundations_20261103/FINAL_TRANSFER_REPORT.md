# Locked cross-stake melt prediction and inverse-physics findings

## Result

The newly developed moisture-anomaly correction did not improve melt prediction at the second stake. All nine models, four declared periods and the original primary model are retained. Its mean absolute period error was **13.353 mm water equivalent**, compared with **7.265** for air positive-degree exposure and **3.435** for a calendar-modulated air baseline. These are nominal errors against the published native camera operator, not confidence bounds.

| Frozen model | Four-period MAE (mm water equivalent) | Cumulative signed bias (mm water equivalent) |
|---|---:|---:|
| `anomaly_exp_air` | 13.353 | +53.414 |
| `air_pdd` | 7.265 | +29.061 |
| `buried_pdd` | 8.359 | +28.203 |
| `gradient_index` | 7.717 | -30.867 |
| `constant_rate` | 18.048 | +72.193 |
| `shifted_ratio_air` | 10.626 | +42.503 |
| `centered_softplus_air` | 13.170 | +52.679 |
| `calendar_exp_air` | 3.435 | -1.305 |
| `moisture_only_exp` | 37.105 | +148.421 |

The primary wins against two of eight baselines and loses against six. Its four errors are positive: +2.668, +10.583, +8.676 and +31.486 mm water equivalent. Air-only exposure has lower absolute error in all four periods. No candidate or baseline was reselected after release. The calendar model's good result does not establish a seasonal physical law or an accurate forecast outside this one campaign.

## What was actually tested

After all four Tower1 periods had been inspected, they were explicitly reclassified as DEVELOPMENT. Five modulation parameters were fitted there with signed grids. This improved the primary's development MAE to 2.813 mm water equivalent, but a matched calendar baseline was already better at 2.253. The primary exponent remained beta=5.4; the softplus, ratio, calendar and moisture-only alternatives had the same development information budget.

At Tower2, every model used the same single local TRAIN ablation amplitude (79.359886 mm water equivalent), and no evaluation-derived parameter. Four subsequent 14-day windows from 16 December 2022 03:00 through 10 February 2023 03:00 were scored. Tower2 moisture increases across the campaign while Tower1's decreases; this transfer challenge was documented before Tower2 outcomes were opened.

This is a locked retrospective test at a second stake of the same glacier and season. The source paper and covariates were already inspected. It is not a prospective field trial or an independent-glacier experiment. The public MAT archive was already locally available. Its decoder materializes native arrays; role gates restrict which endpoint values are returned and used, rather than providing a cryptographic or separate-custodian holdout. The evidence is a published prediction-before-inspection workflow, not proof that outcome bytes were inaccessible. Periods share endpoints, weather and local calibration. Frame counts are not independent replications.

## Auditable sequence

1. Methods, signed-development parameters, source hashes and the role-gated scoring code were published at `b2696672d0ef08c22ed22940b507082711559bc9`.
2. All nine normalized Tower2 predictions were published at `1fc6a5780523ba7b1032e08a0ec818ef3b38b4a1`, before the TRAIN outcome was released.
3. After reading only the two TRAIN endpoints, all nine absolute predictions and their SHA256 were published at `66cc20de0b0f9317f69cbe8d851bfd74ec489487`, before evaluation release.
4. Predictions SHA256: `3e2c64f9059416e8b3da77268b5f1b70ccd17ad10072101d9592e674ef81e225`.
5. Only then were the five declared evaluation endpoints accessed. Independent source auditing confirms exposed-stake distance differences times 917 equal the supplied ablation increments times 1000, to numerical precision. The source README's unit label does not justify a different operator.
6. Independent reconstruction checks raw moisture rows, exact native air integrals, all model doses, gains, predictions and scores. See `moisture/review/ACTUAL_TOWER2_POSTSCORE_AUDIT.json`.

The current package includes the final unsealed results. Pre-release artifacts remain historical snapshots; their statements that Tower2 values were unopened describe their creation time.

## Mathematical and physical work completed

**Conservative heterogeneous conduction.** A finite-volume two-layer model, heterogeneous series-resistance tests and training diagnosis are provided. Layering reduces raw TRAIN melt from 371.8 to 264.0 mm water equivalent; the observed TRAIN melt is 85.2. Better temperature fit is insufficient evidence of correct flux.

**An exact warm-regime non-identification.** Scaling all debris conductivity and volumetric heat capacity by the same positive factor preserves the temperature equation with a melting-temperature boundary but scales the heat flux. The illustrative capacity-conditioned family yields 170.6–432.2 mm with effectively identical temperature predictions. This is an exact model symmetry, not a certified uncertainty interval for the glacier. It is related to established inverse-heat identification; no historical new theorem is claimed.

**Cold-reference identification.** A known independently constrained cold-ice admittance can break that scaling ambiguity, subject to complete nonmelting history and known ice initial/deep conditions. The supplied inversion and interval tests enforce these gates. Piramide lacks resolved cold-interface excitation. NgozumpaNG2 has a resolved cold lower-debris period, but lacks the paired ice-reference conditions and independent melt required for this inversion. Cold data availability alone is not successful melt validation.

**Weak integrated-energy inference.** Positive-conductance estimation works on independently generated synthetic trajectories. On field TRAIN data the inferred deepest conductance changes approximately sixfold with integration window. A necessary LP outer enclosure remains broad and prior-limited. Deepest measured-edge heat is not automatically heat delivered to the ice. This does not support precise field melt forecasting.

**Calibration-invariant empirical index.** Constant additive moisture offsets cancel exactly in the exponential calibrated ratio. Bounded time-varying drift gives a conditional multiplicative bound. Centered alternatives can share offset invariance; arbitrary reference independence distinguishes the exponential only under explicit functional assumptions. These results use classical functional-equation mathematics, do not cancel sensor gain/drift/depth differences, and did not ensure predictive transfer in the new test.

The separate signed-flux/cold-ice phase coupler has conservative energy tests. The empirical moisture index does not inherit energy conservation merely because that coupler exists elsewhere in the package.

## A whole-family failure result, not just a failed coefficient

Post-test diagnosis supplies a stronger statement. With positive-air-temperature exposure as probability weights, let P0 be the observed TRAIN moisture distribution and PI a test distribution. Then M_beta(I)/M_air(I)=E_PI[exp(beta theta)]/E_P0[exp(beta theta)]. Each of the four observed test distributions dominates TRAIN in increasing-convex order, verified by the finite-support stop-loss inequalities. Hence this ratio is at least one for every beta>=0. Because air PDD already overpredicts each nominal target, **no nonnegative exponential moisture correction can improve any of these four period errors** under this same calibration and observation operator. Beta=0 gives the baseline itself.

First-order stochastic dominance alone holds only in periods2–4; the first period's CDFs cross, so the stronger increasing-convex check matters. This is a tested application of established mathematics, not a claim of inventing stochastic ordering. The code does not select a new beta. See `moisture/review/TRANSFERABLE_FAILURE_EXPLANATION.md`.

The conclusion is conditional on observed inputs. Under arbitrary time-varying moisture errors, the all-period symmetric error radius preserving this ordering is only about0.00228 volumetric fraction. The gate fails at the manufacturer's typical-accuracy comparison scale0.031; that scale is not a certified field bound. Thus the result cannot be promoted to an uncertainty-free statement about true glacier moisture or melt. This distinction is computed explicitly rather than hidden in a disclaimer.

## Environmental interpretation

Moisture-aware debris physics is already represented in prior work including Collier 2014 and Giese 2020; these data's source study also uses moisture. The novelty case cannot be 'existing models ignore water or layering.' The candidate was a particular low-input estimator and a rigorous directional transfer comparison. The comparison failed; current evidence does not establish a superior new melt predictor or a Nature-level environmental discovery.

The moisture observations are contemporaneous and may respond to meltwater. This is a state-conditioned hindcast, not an advance weather-driven forecast, causal estimate of wetting effects, satellite inversion, vertical cliff-lifetime prediction or full glacier model. Different probe depths, unknown field calibration gain/drift, nonconductive exchange and latent heat remain potential model-discrepancy sources, not demonstrated explanations of this failure.

The source camera precision and reset uncertainty are not quantified here. Conditional endpoint-error thresholds are sensitivity analyses, not measured confidence intervals. Manufacturer moisture accuracy is typically ±0.031 volumetric fraction under specified conditions, not 3% of the reading and not a Gaussian field standard deviation. Repetition does not divide instrument accuracy by sqrt(n).

## Research decision

Preserve this result, the earlier failed physical models and their baselines. Do not refit beta, choose a different moisture depth, delete unfavorable windows, change density/clock/operator or promote the calendar baseline into the predeclared primary. Tower2 evaluation is now development information for any future hypothesis and cannot become its independent test again.

A subsequent physical prediction claim requires an independently justified heat-scale constraint and a new site/event with compatible geometry, native temperature/moisture, surveyed ablation and quantified measurement support. Adding a neural network or more compute does not supply those observations. Tower3 targets remain unopened; using them would require a separately justified frozen protocol, not repeating model selection until a favorable result appears.

The completed contribution is executable inverse-physics analysis, a calibration-invariance result under explicit assumptions, a new locked field transfer test and reproducible falsification. Awards, publication priority and universal superiority are not outcomes that this evidence can guarantee.
