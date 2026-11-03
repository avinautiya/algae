# Locked Tower 2 moisture-index scoring

`moisture_scorer.py` is additive; it does not edit any historical result or fit any parameter. It reuses the previously tested deferred-truth scoring pattern with a new exact nine-method family, four source-local camera periods, and the fixed Tower 1 development parameters. It is not a standalone authenticated Git publisher: the caller must verify publication of the supplied full commit externally.

1. Build a JSON prediction artifact with `primary_candidate=anomaly_exp_air`, `baseline_ids=BASELINES`, all nine `models`, `primary_metric=mean_absolute_period_ablation_error_mmwe`, and `tie_tolerance_mmwe` fixed before outcomes.
2. Copy `development_parameters=FIXED_PARAMS`, `development_fit_sha256=EXPECTED_DEVELOPMENT_FIT_SHA256`, and `development_parameter_source=Tower1_released_periods_only`. These are checked exactly. The parameter values are empirically learned on the four already-inspected Tower 1 periods, not independent field validation.
3. Supply the exact source target operator: native ice ablation mm w.e., source density917, source-local naive timestamps without an invented timezone, and a documented observation-method label.
4. Call `write_prediction_seal`. Commit and publish the prediction artifact, seal, code and predeclared protocol. Verify the remote commit before passing its40-character identity and the prediction SHA256 to `score_after_freeze`.
5. The deferred truth loader is called only after local locks and schema are validated. It must return the same source/operator and four observations under `observed_mmwe`. Missing/nativeNaN values fail scoring; they are not silently dropped or interpolated.

Prediction rows use `start`, `end`, `predicted_mmwe`. Truth rows use `start`, `end`, `observed_mmwe`. All intervals must exactly match the four December16T03→December30T03→January13T03→January27T03→February10T03 periods declared in the module.

Eight synthetic tests pass. They check locked family, fixed development parameters, tamper-before-truth refusal, exact camera support/density, missing-value refusal, and the primary MAE against cancellation. No Tower 2 outcome was accessed during their creation.

Interpretation is nominal prediction comparison at one previously withheld stake. Shared camera endpoints and serial intervals do not supply four independent glaciers. No confidence interval, professional energy-balance superiority, universal superiority, Nature significance or award probability is computed.
