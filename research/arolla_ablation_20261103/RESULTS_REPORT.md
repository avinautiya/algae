# Native field-ablation test and uncertainty audit

The frozen positive-degree-day hypothesis fails its predeclared primary criterion. This package now compares predictions with actual field-measured ice ablation, rather than a sensor proxy, but does not establish exceptional accuracy, superiority over professional models, or a novel field-validated mechanism.

## Freeze and independent reproduction

Pre-execution code/data/protocol were published at commit `bc7916b852b9667820ef68e7442d4759e06ff342` before computing predictions. `model/results/frozen_ablation_comparison/RESULTS.json` retains all 1,053 rows: 972 native-target comparisons plus 81 contextual professional-output comparisons. No target interval fitted a coefficient. Original target values were inspected during source auditing, so this is a retrospective temporal holdout, not a blind experiment.

`targets/POSTSCORE_REVIEW.json` independently recomputes all rows directly from native workbook/weather/model CSV sources, without importing the predictor. Maximum absolute difference is 6.94e-12 mm water equivalent. This verifies implementation and source binding, not predictive superiority.

## Primary result

H1D is the source-designated nearest AWS stake, with known 6 cm debris. Its first interval provides 846 mm w.e. for calibration; the evaluation interval contains 315 mm w.e. measured ablation.

| Method | Nominal prediction (mm w.e.) | Absolute error (mm w.e.) | Timing range (mm w.e.) |
|---|---:|---:|---:|
| Constant rate | 382.06 | 67.06 | 343.69–423.00 |
| Positive degree day, primary | 229.30 | 85.70 | 213.41–264.58 |
| Net-shortwave dose | 367.45 | 52.45 | 333.91–405.75 |
| Published DEBCF, unresolved component | 215.58 | 99.42 | 204.39–237.24 |

The predeclared primary test requires the degree-day model to beat both simple baselines at H1D under every declared timing convention. It does not. Net-shortwave dose has the smallest nominal H1D error among the three simple methods; selecting it after inspecting this test would not demonstrate a successful predeclared candidate. Across four stakes, degree day has the smallest nominal pooled MAE (91.78 mm w.e., versus 95.09 shortwave and 104.81 constant rate). That secondary comparison does not replace the primary result. Four stakes share one weather episode and are not four independent environmental replications.

DEBCF is contextual only: the exported ice-only versus ice-plus-snow component and parameter independence remain unresolved. Its larger nominal discrepancy cannot establish professional-model superiority.

## Robustness to field uncertainty

Visit dates are known but hours are not. All 27 combinations of 0/12/24 h for three visits are evaluated under three hourly-bin conventions. A shared middle visit is used consistently in calibration and evaluation. These are sensitivity scenarios, not independent samples.

The source states “Uncertainty=0.1m” without specifying its probability distribution or whether it describes readings or interval differences. With source density 900 kg/m³, this corresponds to 90 mm w.e. under either explicitly conditional bounded interpretation. No Gaussian error or confidence level is invented.

For a predicted test/train ratio f and measured calibration/evaluation totals A/B, the residual is fA−B. If the bound δ applies independently to each of the three shared readings, its error radius is 2(f+1)δ. If it instead applies to each interval total, its radius is (f+1)δ. These describe sensitivity to uncertain inputs and recalibration: f times a possible true calibration total minus a possible true evaluation total. They are not the errors of an already-issued fixed forecast, which has a separate audit in reproduction/robustness/. Calibration-input sensitivity and fixed-forecast verification must not be mixed.

`timing/FIELD_RANKING_AUDIT.json` tests all three model pairs across 81 scenarios. All 243 rankings remain unidentified under either interpretation. This does not prove all models equally accurate. It means the available source uncertainty cannot support a robust ranking under those assumptions. No interval is labelled a statistical confidence interval.

The independent exact-corner audit in `reproduction/robustness/` also distinguishes input-error sensitivity from fixed issued forecasts. The smallest worst-case input-sensitivity loss belongs to degree day; the smallest worst-case issued-forecast loss belongs to net shortwave. Neither minimax ordering proves lower actual field error or robust pairwise dominance. All pairwise rankings can reverse within the declared measurement boxes.

`timing/MEASUREMENT_DESIGN.json` quantifies the resolution obstacle. At the nominal clock, distinguishing the shortwave-versus-degree-day ranking under the conditional per-reading model requires a reading bound strictly below about 6.8 mm of ice, compared with the source's unspecified 100 mm uncertainty. Even zero measurement error does not preserve that pair's ranking across every declared clock. This is a retrospective design diagnostic, not a statistical power calculation or a request to discard unfavorable clocks.

## Environmental scope and next model

The empirical comparison uses observed air temperature, measured net shortwave and elapsed time. It does not represent debris heat storage, moisture transport, turbulent exchange, evolving snow or channel/cliff contact. It is a weather-conditioned ablation hindcast, not a forecast of future weather. Complete temperature, humidity, wind and four radiation streams are available over the 1,104-hour target envelope, but snow depth and pressure are missing throughout. Missing snow must not be replaced by bare debris.

`model/environment/` audits input availability and specifies a phase-aware transient conduction model with independent thermal parameters, a signed energy/storage ledger and explicit input refusals. Adding these established physical terms is not itself novel. Any new model must be tested on distinct support; this inspected Arolla interval is not available for repeated blind validation.

`reproduction/PROCESS_COMPARISON_AND_GATE.md` identifies a useful unscored seasonal forcing contrast at Pirámide: later/earlier dose ratios differ between degree day (0.722), shortwave (0.545) and net radiation (0.424). This makes it a process-discriminating setting, not evidence of improved predictions. Snow state, thermal parameters and observation timing must be audited before scoring.

`model/thermal_memory/` adds an executable transient finite-volume debris column. It integrates the same basal flux into ice ablation and the energy ledger, preserves debris heat storage, uses independently sourced thermal parameters and an explicit initial profile, and returns finite-scenario uncertainty envelopes. It is established conduction physics and currently a numerically verified component, not a field-superior model. Negative boundary temperatures, missing source declarations and unsupported refreezing are refused rather than silently replaced. This new predictor does not modify the frozen Arolla comparison.

## Reproduce

Run `python3 verify_ablation.py` for component checks. Reproduce the original frozen scoring with the exact command in `model/PROTOCOL.md`, retaining its original pre-execution commit. `python3 make_results_figures.py` renders figures only from frozen results. It never fits a coefficient. The figures show nominal errors, declared timing ranges and conditional ranking uncertainty; none displays an invented confidence interval.
