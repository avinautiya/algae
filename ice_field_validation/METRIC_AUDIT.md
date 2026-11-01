# Independent saved-event metric audit

PASS: saved-event metrics independently recomputed; no model rerun or retuning.

Recomputed every reported metric directly from saved observed and predicted event durations using Python standard-library arithmetic. No model was rerun, selected, or retuned. Input hashes and full recalculations are in `METRIC_AUDIT.json`.

| Test | Complete paired events | Site counts | Largest metric difference |
|---|---:|---|---:|
| 2022 | 105 | P_IC_E: 38, P_IC_W: 67 | 2.66e-15 |
| 2023 base | 71 | P_IC_E: 7, P_IC_N: 21, P_IC_W: 27, P_IC_W2: 16 | 2.66e-15 |
| 2023 SW_minus15percent | 71 | P_IC_E: 7, P_IC_N: 21, P_IC_W: 27, P_IC_W2: 16 | 1.78e-15 |
| 2023 SW_plus15percent | 71 | P_IC_E: 7, P_IC_N: 21, P_IC_W: 27, P_IC_W2: 16 | 3.55e-15 |

All four models have finite positive predictions on every eligible saved event: 105/105 in 2022 and 71/71 in every 2023 scenario. There is no model-specific censoring or complete-case selection within these eligible sets. The original eligibility exclusions remain; this audit does not imply all raw stake rows were eligible.

Pooled scores weight each event equally. Site weights and alternative equal-site MAE, bias, and square-root mean site MSE are reported separately in JSON; averaging site RMSE is not the pooled RMSE. Pooled medians likewise cannot be averaged from site medians.

Melt calibration dictionaries are identical across saved 2022 and 2023 results and the frozen 2021 calibration (94 intervals). Holdout is by year, not wholly new sites: training includes E/N/S/W, 2022 tests E/W, and 2023 tests E/N/W/W2. Frozen protocol hashes match. This verifies saved artifacts, not independent reconstruction of calibration history or preregistration. Weather transfer is fit on 2021–2022 paired weather; 2023 is a secondary held-out melt year after 2022 scores were inspected.

- Pooled metrics weight events equally, not sites equally; adjacent events share endpoints and weather.
- 2022 has two site-year histories; 2023 has four, with repeated sites. These are not 176 independent glaciers or validations.
- 2023 weather transfer uses 2021+2022 meteorology, not 2023 melt fitting; 2022 results were already inspected before secondary 2023 evaluation.
- All results are retrospective observed-weather atmospheric-component hindcasts, not issued forecasts.
- Missing hydrological heat is unknown, not measured zero. These models do not validate the new coupled hydraulic-memory/debris/stream system.
- SW +/-15 percent are engineering scenarios, not uncertainty intervals. No confidence interval or accuracy guarantee is inferred.
- Saved events retain no 2023 timestamps; 2023 temporal membership relies on protocol and reviewed evaluation code rather than event-level dates.
