# Locked three-period batch scorer

`batch_scorer.py` is additive. It does not edit the earlier single-period scorer, specifications or results. No future camera targets were opened while implementing it.

## Same API, new fixed primary metric

Use `write_prediction_seal(predictions_path, seal_path)` and later `score_after_freeze(predictions_path,seal_path,expected_prediction_sha256=LOCK,freeze_commit=FULL_COMMIT,truth_loader=DEFERRED_LOADER)`. Publication must be verified externally. Hash/source/family checks occur before the future truth loader is called.

Prediction schema remains the earlier schema, except:

- `primary_metric='mean_absolute_period_ablation_error_mmwe'`.
- Exactly3 contiguous source-camera intervals must exist for every model.
- All39 frozen model IDs are required, including36 thermal geometry×hold×capacity scenarios and3 simple baselines.
- The primary is fixed to `forward__middle_gap_20mm__left__proxy`; the same4 primary baseline IDs are required.

The source/operator, density917, local-naive clock, exact common support, signed increments, code/hash lock, complete secondary family and tie-tolerance checks remain. Unknown camera precision is not converted into an invented error bar.

Primary error is **mean over3 periods of |predicted ablation−observed ablation|**. Each model also reports its3 period errors, signed cumulative bias and absolute cumulative bias as secondary context. A model whose positive/negative errors cancel cannot win merely because its summed bias is small. All comparisons/ties and all39 models are retained.

These are serial periods from one glacier with shared endpoint readings. The result explicitly supplies no confidence interval, bootstrap or independent-glacier claim. The fourth proposed future period is excluded for **all** methods before outcome access because the air-temperature source ends before its requested endpoint; it is not silently extrapolated or selectively scored for the candidate.

## Preserve the known first comparison

After batch scoring, call:

`pooled_with_prior_result(batch_result, prior_v3_results_path, expected_prior_result_sha256=PREDECLARED_HASH)`

This verifies the pinned earlier result, same dataset/operator, identical39-method family and adjoining clock support. It reports **secondary four-period MAE**, including the earlier V3 period where the candidate lost to air-temperature PDD. It does not replace the frozen3-new-period primary metric or hide that earlier result. The source ID should identify the same dataset/stake in both artifacts; calendar windows are recorded separately.

Eight synthetic tests pass, covering39-method retention, fixed primary, exactly3 periods, errors that cancel in total but lose on primary MAE, ties, hash/code locks before target loading, missing endpoint refusal, native operator identity and pinned four-period context. No actual future outcome was scored in these tests.
