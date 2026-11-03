# Independent Tower 2 result audit

The predeclared moisture-modulated primary **does not improve held-out melt prediction over the declared baseline family**. Its arithmetic and source reconstruction are correct. All 3,524 independent checks pass, so the observed failure is not explained by a scoring mismatch or a failure to apply the locked local calibration.

Final prediction publication: `66cc20de0b0f9317f69cbe8d851bfd74ec489487`.
Prediction SHA256: `3e2c64f9059416e8b3da77268b5f1b70ccd17ad10072101d9592e674ef81e225`.
Evaluation access was authorized by the parent only after that lock. This audit reads the already-scored results; it does not independently reopen any camera archive or select additional outcomes.

## Held-out scores

Four prespecified serial periods at Tower 2; lower MAE is better, units mm w.e.

| Locked method | Four-period MAE |
|---|---:|
| Calendar-modulated air PDD | 3.434992 |
| Air PDD | 7.265305 |
| Gradient index | 7.716761 |
| Buried-temperature PDD | 8.358887 |
| Shifted moisture-ratio air index | 10.625866 |
| Centered-softplus moisture air index | 13.169711 |
| **Primary exponential moisture air index** | **13.353386** |
| Constant rate | 18.048331 |
| Moisture-only exponential index | 37.105263 |

The candidate nominally wins 2 of its 8 predeclared comparisons and loses 6. It overpredicts every period: +2.668458, +10.582845, +8.676186 and +31.486055 mm w.e. Total overprediction is 53.413544 mm w.e. The one local TRAIN ablation was 79.359886 mm w.e. The training amplitude was fixed separately for every method; none received an evaluation-derived gain.

## What was checked independently

The auditor does not import `anomaly_index.py` or call its dose function. It verifies the prediction/seal/result lock and published commit identity supplied by the parent; exact model family, fixed primary, development-fit parameters and source hashes; all canonical SM3 values against their raw decoded native source; native irregular air integration using a separate positive-piecewise-linear integral; raw profile-based gradient/buried-temperature doses; every one of nine model training gains and 36 interval predictions; every MAE, bias and pairwise status.

Constant additive moisture offsets cancel algebraically for the primary as intended. The centered-softplus variant also satisfies its centered-offset property. This correct mathematical feature did not create field-predictive superiority. It does not remove sensor drift, multiplicative calibration, different sampling depths, missing physics or seasonal confounding.

## Interpretation

The prospective directional-transfer test does not support a universally positive, transferable moisture modulation learned from Tower 1. This does NOT prove moisture has no effect on glacier melt. The covariate is measured at different depths, uncertain calibration remains, and moisture affects heat capacity, phase and latent/advective flux as well as conductivity. The simple index has no physical energy-closure claim. Tower 1's strong development improvement did not transfer; the predeclared calendar baseline performed substantially better on Tower 2.

The independently checked result supports retaining the negative transfer finding and rejecting the primary estimator's current accuracy claim. Do not reverse beta, choose another sensor, adjust periods or promote the best alternative using these same Tower 2 outcomes and call it a new independent validation. Any later model developed using these outcomes must label both stakes as development data and obtain a genuinely new test.

## Conditional measurement sensitivity

The camera source does not establish a quantitative simultaneous error bound, reset uncertainty or probability model. Therefore no achieved uncertainty, confidence interval or statistical significance is claimed.

For fixed issued forecasts and a hypothetical per-camera-reading error bound delta in mm ice, a conservative pairwise MAE perturbation bound is 4*0.917*delta mm w.e. Because the calibration TRAIN endpoint is also shared with the first evaluation interval, a second calculation permits local TRAIN gains to change consistently under the same six bounded reading errors, while retaining all forcing ratios and fitted beta values. That conservative bound is 0.917*mean[2*r_candidate+2*r_baseline+4]*delta, with frozen positive dose ratios r.

| Nominal result | Sufficient maximum endpoint error for ranking to persist with TRAIN recalibration, mm ice |
|---|---:|
| Air PDD better than primary | 0.872430 |
| Calendar index better than primary | 1.457973 |
| Gradient index better than primary | 0.849904 |
| Buried-temperature PDD better than primary | 0.716227 |
| Primary better than constant rate | 0.649595 |
| Primary better than moisture-only | 3.097571 |

These thresholds are sufficient, conservative sensitivity limits, not exact worst-case optima and not measured camera precision. Shared readings were handled coherently in the derivation and an independent 64-corner synthetic test. Serial intervals at one stake are not four independent glaciers. A large nominal loss is not made into statistical proof by assuming unknown camera precision.

## Artifacts

- `ACTUAL_TOWER2_POSTSCORE_AUDIT.json`: all checks and every conditional ranking threshold.
- `postscoreaudit.py`: independent reconstruction and score audit.
- `test_postscoreaudit.py`: five synthetic tests, including shared TRAIN/evaluation endpoint perturbations.
- Original frozen results and unsuccessful hypotheses remain unchanged.
