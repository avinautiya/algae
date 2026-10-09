# Supplementary ML comparison

This additive study does not edit existing physics, production jobs, observations or held-out records. Starting commit: `aaff7ca2e275d054dae902e6655eb0fcee0899a9`.

Read [PROTOCOL.md](PROTOCOL.md) before interpreting any result. The study compares a training mean, ridge, random forest and Gaussian process on both inverse abundance retrieval and forward four-band **ground HCRF** prediction. This is retrospective supplementary research, not new confirmatory evidence, an agency-ready satellite product, or a glacier-melt validation.

## Run and test

From the repository root, using a Python environment with the listed dependencies:

```sh
python -m unittest discover -s comparison_study/tests -v
python -m comparison_study.study --outdir records/my_new_ml_run
# Resume only the same inputs/code/dependency versions:
python -m comparison_study.study --outdir records/my_new_ml_run --resume
```

No PySCF/BioSNICAR or expensive chemistry run is needed. The existing field loader performs the same SRF averaging used in the project's field validation. One numerical thread and one tree worker are used. Completed task/model/fold prediction bundles are checkpointed and their signatures checked on resume.

Outputs include per-sample/band predictions, complete inner-CV scores, chosen configurations, training-range extrapolation flags, matched positive-count physics point comparisons, source-data hashes, versions, runtime and maximum RSS.

## Correct interpretation

- Positive-count target only: zero counts are inventoried and excluded explicitly, not replaced with a pseudocount. Continue using the project's Poisson-count experiment for all-count predictive scores. Do not compare these RMSEs directly with its log predictive score.
- Outer transfer is between site/year datasets; inner tuning holds out entire training sampling days. Scalers are fitted separately inside each inner fold. Already-examined sites are not untouched confirmatory evidence.
- Forward models receive observed noisy abundance, not latent exact abundance. Their outputs are four HCRF bands, not broadband hemispherical albedo. No forward RMSE superiority over physics is claimed because matching physics per-band point predictions are not supplied by existing records.
- Fewer than five evaluation days produces **no bootstrap confidence interval**. Intervals with more days describe conditional day-block sampling uncertainty, not transfer across independent sites.
- Physics comparisons use the existing reported medians versus ML squared-error point predictors. They do not compare complete predictive distributions or propagate molecular uncertainty.
- Simple min/max training-range flags are warnings, not a calibrated applicability-domain classifier. All extrapolations remain reported; they are not discarded to improve scores.
- No outer-test winner is promoted as a validated production model. Final selection would require additional independent data.

## Claude handoff

Review/cherry-pick this branch after reconciling concurrent work. It adds only `comparison_study/` and new `records/ml_comparison/` outputs. Existing pipeline files are untouched. Use the saved results as a supplementary benchmark, not a replacement for frozen physics evidence.

The next useful integration is matched forward per-band physics predictions on exactly the same sample set, followed by independent albedo/absorbed-SW validation. Do not assume these small field-trained ML models transfer directly to Sentinel-2 pixels: atmospheric, geometry and spatial-support differences remain.
