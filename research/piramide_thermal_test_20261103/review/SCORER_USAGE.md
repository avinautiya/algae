# Locked native-camera scorer

This scorer does not read any real evaluation camera target by itself. It accepts a deferred loader supplied by the data owner, after checking an externally locked prediction hash, a prediction seal and its own code hash. A commit ID is recorded, but this function does not contact GitHub or prove when that commit was published. The parent must verify the actual published freeze before calling it. Declarations cannot prove a target was never inspected.

## Order

1. Freeze the protocol, train-only calibration/selection procedure, parser, camera archive CRC/source hash and code before any permitted training target access.
2. Fit using the permitted training camera block only; freeze every actual evaluation prediction and the predetermined primary/secondary model identities. Do not inspect evaluation outcomes.
3. Call `write_prediction_seal(prediction_path, seal_path)`. Publish/verify the predictions, seal, code and protocol; independently save the actual prediction SHA256 and full freeze commit ID.
4. Only then call `score_after_freeze(..., truth_loader=deferred_loader)`. The loader must independently validate raw source CRC/hash and the frozen camera observation operator. The scorer enforces matching returned source/operator/support, not raw data authenticity by itself.

## Prediction JSON

Required fields:

- `source_id`: frozen source identifier; include version/member identity rather than a generic glacier name.
- `heldout_targets_seen_before_prediction_freeze`: false, explicitly.
- `target_operator`: the same object that the frozen camera loader will return. Required entries: `quantity='native_ice_ablation_mmwe'`, `observation_method`, `ice_density_kg_m3=917`, `clock_basis='source_local_naive_no_utc_claim'`. Add the source CRC/hash, footprint/stake ID and exact operator convention here so equality is checked too.
- `primary_candidate`: one predetermined model ID.
- `baseline_ids`: nonempty distinct IDs for the entire predetermined primary baseline family.
- `secondary_model_ids`: optional remaining geometry/parameter variants. All declared IDs must remain in `models`; none may be silently dropped. Only the primary baseline family determines the nominal primary win flag.
- `primary_metric='absolute_total_ablation_error_mmwe'` and fixed `tie_tolerance_mmwe` (normally0).
- `models`: ID→list of `{start,end,predicted_mmwe}`. Every model must have exactly the same contiguous ordered calendar intervals, with no gaps or overlaps. Source local-naive timestamps remain naive; no UTC is invented.

Truth loader returns `{source_id,target_operator,intervals:[{start,end,observed_mmwe}],uncertainty_semantics:'unknown'}`. Its support must exactly match predictions. Source density917 is preserved through exact operator equality; changing to900 fails. All finite signed observed increments are retained, including camera-noise negative increments; they are counted, not clipped. Negative predictions are also retained and flagged, rather than omitting a poor model.

## Metrics and claims

Primary error is **|sum(prediction)−sum(observation)| over one entire held-out calendar block**. Secondary interval MAE is mean|prediction−observation|; signed cumulative and mean interval biases are also reported. Cumulative errors can cancel while interval MAE remains high, so both remain visible. Zero differences under the frozen tolerance are ties, not wins. An empty family is refused; `all([])` cannot manufacture success.

All model metrics and primary comparisons are returned. Nominal victories do not produce statistical confidence intervals, independent-frame bootstraps, a universal-superiority claim or a professional full-energy-balance victory. A many-frame camera record is explicitly one evaluation calendar block here. Unknown source uncertainty semantics stay unknown; the scorer does not invent measurement-error bars.

Nine synthetic tests pass, including deferred-loader ordering on corrupted artifacts, scorer/family integrity, same support/native density, metric cancellation, ties, geometry-family retention, negative observations and no invented timezone. No real evaluation camera targets were opened by these tests.
