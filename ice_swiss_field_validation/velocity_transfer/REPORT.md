# Retrospective conditional velocity transfer test

## Result

**The piezometric-head correction does not establish a superior transferable point predictor in this test.** It improves the three-case training fit, but the fitted constant-friction model underpredicts all four subsequent tracer velocities and performs worse than both the archived-plus-head model and the training constant-velocity baseline. The result holds at every prespecified width. No split, parameter, matching rule or objective was changed after these scores.

The protocol was saved before first prediction scores (SHA256 `713b5035c19644bcf32ff21668d4fd63454a939ddf134fddaffa48bbe1edaa34`). All target observations had already been inspected in the source audit. This is a retrospective temporal comparison, not a blinded held-out discovery or independently timestamped preregistration.

## Frozen method

The first three eligible experiments (11 July morning, 11 July midday, 16 July) train each one-parameter constant-Darcy-friction model. The next four (25 July, 30 July twice, 31 July) test it. P5/P3 floor elevations and independent depth observations/deductions are imposed externally. Hydraulic diameter is `4*w*hbar/(w+2*hbar)`, where `hbar=(h5+h3)/2`. Neither test velocity nor discharge enters this diameter. Both candidate head expressions fit the same velocity-squared-error objective on exactly the same training cases. Positive amplitude `a=1/sqrt(f)` is fitted analytically. Widths 1.5, 2.0 and 2.5 m are all reported, not selected.

For 11 July the upstream 0.8 m depth is an author deduction, not a manual observation. The midday downstream depth is the nearest same-day manual value 0.30 m rather than the earlier 0.23 m chosen by the original notebook's first-same-day rule; that matching rule was fixed in the protocol before scores. On 31 July downstream 0.38 m is a pressure-sensor observation. The remaining imposed depths are manual measurements. These spatial/temporal representation assumptions are shared between candidate models.

## Test metrics at primary width 2 m

| Predictor | MAE m/s | RMSE m/s | Bias m/s |
|---|---:|---:|---:|
| Piezometric head, training-fitted f | 0.3171 | 0.3393 | −0.3171 |
| Archived-plus head, training-fitted f | 0.1289 | 0.1567 | −0.0967 |
| Training constant velocity | 0.1750 | 0.2228 | −0.1215 |
| Piezometric head, fixed f=0.3 | 0.4154 | 0.4361 | −0.4154 |
| Archived-plus head, fixed f=0.3 | 0.1171 | 0.1418 | −0.0732 |

The f=0.3 references are literature/author-context choices, not independently estimated external priors proven uninformed by this study. The main comparison is the two equally trained one-parameter models. Piezometric fitted f is 0.1788; archived-plus fitted f is 0.3220. Training squared error is 0.0430 versus 0.1010 (m/s)² respectively, illustrating that improved training fit does not demonstrate temporal transfer.

| Test experiment | Measured m/s | Piezometric fitted | Archived-plus fitted | Training constant |
|---|---:|---:|---:|---:|
| 25 July | 1.015 | 0.5973 | 0.7710 | 0.627 |
| 30 July 14h | 0.583 | 0.3757 | 0.6475 | 0.627 |
| 30 July 18h | 0.832 | 0.3757 | 0.6475 | 0.627 |
| 31 July | 0.564 | 0.3769 | 0.5413 | 0.627 |

Across widths 1.5/2.0/2.5 m, piezometric fitted MAE is 0.3180/0.3171/0.3164 m/s; archived-plus fitted MAE is 0.1291/0.1289/0.1288 m/s. The constant baseline is 0.175 m/s throughout. The predeclared consistent-superiority rule fails at all widths. Complete per-case signed/absolute errors for all models, widths and training/test cases are preserved in `CASE_PREDICTIONS.csv`; all metrics and coefficient values are in `RESULTS.json`.

## Interpretation

An algebraically inconsistent head expression can compensate for an inadequate constant-friction or representative-geometry assumption. Better empirical scores for that expression do not make its algebra physically correct. Conversely, fixing algebra is necessary for a physically interpretable model but is insufficient evidence of improved field prediction. This small comparison does not isolate the cause of the late-case underprediction.

The two 30 July tests receive identical imposed geometry/head inputs in this protocol, while their measured mean velocities differ by 0.249 m/s. Any deterministic predictor with only these identical inputs must predict the same value for both; its summed absolute error for this pair is at least 0.249 m/s by the triangle inequality. Thus pair MAE cannot be smaller than 0.1245 m/s without additional time-resolved information. This is a simple data-resolution bound, not a new mathematical theorem. It does not justify retroactively changing this test's inputs or fitting a time-varying coefficient against its exposed targets.

A future method would need independently available dynamic geometry, hydraulic state or friction information and a new unexposed validation set. No such model is established by this test.

## Uncertainty and scope

There is one channel segment and four dependent test experiments, two sharing a day and geometry. No confidence interval or significance claim is inferred. Corrected nominal heads are small (0.23 m on 30 July, 0.28 m on 31 July) relative to original floor/depth uncertainties and unknown covariance. Same-day non-simultaneous stage, width, reference-frame motion, local versus reach-mean conditions and neglected velocity-head terms remain material. The original timezone labels remain inconsistent; this same-day comparison does not resolve that issue. Width sensitivity alone is not full uncertainty propagation.

This is a conditional quasi-steady hydraulic hindcast driven by supplied surveyed/deduced heads and depths. It neither predicts those states nor tests contact loss, thermal melt, glacier forecasts, or the new full coupled cliff model. No existing author products, previous field results, or Claude worktree were changed.

## Verification and reproduction

Six synthetic/property tests pass, including the geometry-only diameter interface, known-amplitude recovery, head scaling, invalid-input refusal and metric arithmetic. Run:

```sh
python -m unittest discover -s . -p test_velocity.py
python run_velocity_test.py
python plot_velocity_test.py
```

Execute from this folder in a copied delivery to preserve historical outputs. `run_velocity_test.py` expects the sibling original_mirror observation subset and the parent saved tracer-velocity CSV. It needs only Python standard library; plotting adds NumPy and Matplotlib. The PNG/PDF figure shows the full test predictions and every width sensitivity.
