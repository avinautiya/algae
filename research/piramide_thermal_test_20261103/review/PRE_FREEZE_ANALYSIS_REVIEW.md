# Pre-freeze analysis review

Reviewed `analysis/run_melt_comparison.py`, `ANALYSIS_SPEC.json`, `calibration.py`, role-gated camera loader and model freeze before camera outcome access. No outcomes were opened by this review.

## Correct as written

- Normalized basal heat has units K·m. Multiplying by bulk volumetric heat capacity J/(m³K) gives J/m²; dividing by latent heat J/kg yields kg/m², numerically mm water equivalent. **No additional factor1000 is required.** A gain times latent heat has the units of volumetric heat capacity, but is labelled effective predictive gain, not an independent material measurement.
- Each model uses the same cumulative training-camera endpoint difference. No evaluation outcome is passed to the gain estimator. The camera loader releases only authorized role endpoints after its receipt checks.
- Constant rate and positive-temperature baselines also use one training gain. The buried-DT1 baseline supplies a particularly relevant simpler same-input comparison. Air-temperature reconstruction is separately named and its18h maximum source gap is explicit; it is not fully native forcing or proof of superior physics under an equal information budget.
- The positive-part linear interpolation integral is correct, including a zero crossing. Surface-temperature conditioning is accurately labelled hindcasting.
- Gain diagnostics and alternative geometries are not used to choose the primary model. Primary cumulative error, source density917 and one held-out calendar block agree with the scorer.

## Fix before the first outcome access

1. Bind the actual `normalized` argument to the receipt's published prediction hash, not merely to a separately validated receipt filename. Otherwise the loader can validate one frozen artifact while `fit()` reads another.
2. Bind the actual `spec_path` to the published global manifest or an explicit expected analysis-spec hash. Merely recording its hash after training is not verification against the locked spec.
3. Before calling `load_after_freeze()`, explicitly require normalized start/split/end to equal `primary_train` and `primary_evaluation`. Authorized-label refusal usually catches mismatches, but the analysis split should be checked directly before outcomes are loaded.
4. Replace the catch-all method-name fallback with a whitelist of the exact two method names: `training_temperature_forward_fit` and `training_derivative_OLS_NOT_source_NTLS`. Refuse duplicate generated IDs before assignment to `models`; dictionary overwrites would silently discard a frozen family member.
5. Verify decoded feature-file hashes against the global freeze before integration. Hashing them only for an output record does not establish that the input was the frozen source.

These are integrity checks, not alternative choices of dates, methods or parameters. Do not inspect outcomes to decide whether to make them.

## Remaining interpretation limits

One train target calibrates one amplitude for each model. Any improvement therefore tests transfer of the temperature-derived temporal response ratio into the later block, not a measured heat capacity. A low nominal cumulative error can coexist with poor subinterval dynamics. Here only one evaluation endpoint pair is being released, so interval MAE equals the primary absolute error and cannot establish short-time dynamics.

A predeclared interpolation and finite geometry family are defensible sensitivities, not continuous uncertainty coverage. The source's camera precision/reset semantics remain unquantified. No nominal winner in this one-block comparison establishes statistical, multi-glacier, full-energy-balance or universal superiority.
