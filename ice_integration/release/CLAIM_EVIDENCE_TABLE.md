# Corrected claim–evidence table

Evidence classes: **NV** numerical verification (code vs analytic/independent reference, conservation, convergence); **RF** retrospective field comparison (model vs inspected field data that were visible during development; not blind); **CD** conditional demonstration (model run with imposed/assumed inputs; shows behaviour, tests nothing against nature); **IV** independent validation (blind/untouched observations, model frozen first). **No claim in this package reaches IV.**

| # | Claim | Class | Evidence (file) | Status / limit |
|---|---|---|---|---|
| 1 | Integrated point predictor reproduces Codex V5 predictions | NV | `results/EVALUATION_RECORD.json`: max rel. diff 4.2e-16, n=5994 | Verified; equivalence only, adds no physics |
| 2 | Finite-heat accounting conserves water and energy; stable depletion; restart-equivalent | NV | 15 tests (`test_integration.py`); demo residuals ≤1e-12 | Verified for implemented fluxes only; no atmospheric exchange, moving geometry or cold-ice conduction |
| 3 | Plug-flow limit of coupled solver equals finite-heat predictor (≤1 %) | NV | integration tests | Verified |
| 4 | Calibration offsets applied once; support durations split observed/reconstructed/unsupported | NV | tests; `forcing.py` | Verified for the P3/P5 sources; clock convention unresolved |
| 5 | Candidate (scale-fitted DB0.4) beats every baseline consistently | RF | `EVALUATION_RECORD.json`, `INTEGRATION_REPORT.md` | **FALSE**: consistent-improvement verdict = False (P5 vs constant-rate, DB0.4, DB0.3, Gnielinski; P3 vs constant-h) |
| 6 | Candidate pooled MAE 0.585 m; every one of 81 scenarios over-predicts (+0.52…+0.82 m) | RF | `release/ERROR_BUDGET_NUMBERS.json` | Descriptive; 6 test intervals; ranges are engineering sensitivity, not CIs |
| 7 | Local geometry reduced one family's pooled MAE 23.6 % but did not consistently beat established baselines | RF | Codex RESULTS_V5_DEVELOPMENT.md (inspected, preserved) | Retrospective; development data |
| 8 | Strict-native scoring | — | No complete interval has full native support | **No result exists** |
| 9 | Scale parameter identifiable for P5 | — | `STATION_DIAGNOSIS.json` | **Not identifiable**: one training interval, SNR<1; P3 effective n≈1.3 |
| 10 | Pressure stage is full depth | — | `CLEARANCE_CHECK.json` | **Unsupported**: differs 0.02–0.58 m (indicative) from independent depth; clearance unavailable |
| 11 | Coupled plan-view model, imposed Q, strip deposit: toe contact lost after ≈80 h, melt volume grid-converged to ≤1.5 % | CD (with NV) | `results/demo/`, `results/convergence/` | Idealised geometry, assumed Q, no observed target; contact timing ±~10 % at dx 0.05 |
| 12 | Wet/dry threshold does not affect toe-contact result | NV | `results/convergence/` | This strip geometry only |
| 13 | Connectivity-aware reduced model beats fractional-cover and fixed-geometry models at predicting contact-loss time and partition vs the depth-averaged CFD | NV (model-to-model) | `ice_cliff_survival/reactivation/REPORT_M3_FINAL.md` | Common thermal closure and imposed inputs; tests hydraulics/connectivity, not insulation; 10 configurations; 3 of them not held out in configuration |
| 14 | Natural boulder (porous) behaviour | — | Ergun-drag closure unverified | **Unsupported**; E08–E10 are closure-dependent |
| 15 | Rochemelon independent test | — | `rochemelon/` | **Not run**; all gates unresolved; paper's coefficients were tuned to this event |
| 16 | Prediction of ice-cliff retreat, melt volume or lifetime in the field | — | none | **Unsupported** |
| 17 | Independent validation of any component against nature | IV | none | **None exists** |
