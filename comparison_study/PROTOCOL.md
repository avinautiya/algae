# Supplementary ML comparison protocol v1

Defined before running this implementation, but **after existing held-out physics results were known**. This is a retrospective supplementary analysis, not a new preregistration or untouched external validation. Do not change candidates after viewing these results and retain the same label.

## Questions and targets

1. **Inverse:** predict observed log10 abundance from four SRF-averaged ground HCRF bands B2/B3/B4/B8, conditional on positive observed counts. Primary descriptive endpoint: sample-weighted RMSE in dex.
2. **Forward:** predict the same four ground HCRF bands given observed log10 abundance. Primary descriptive endpoint: RMSE pooled across samples and bands; also report per-band RMSE and bias. No solar-angle input is used in v1 so these are direct comparisons to existing statistical forward models, not radiative-transfer substitutes.

These are field-spectrum tests, not actual satellite validation, hemispherical-albedo prediction or melt validation. Observed abundance is noisy; forward inputs are not latent exact cell counts. Between-site instrument/geometry differences may contribute to error.

## Data and split

Use the repository's `field_validation.field_band_reflectance` conversion without changing source observations. Export every retained field row and its dataset/sample/day key. Explicitly list excluded zero counts: 5 S6 samples have no finite log10 abundance, and no pseudocount is introduced. The existing Poisson-count experiment remains the all-count evaluation; scores here cannot replace its result.

Outer folds remain S6 2017 → southern Greenland 2021 and reverse. Inner selection is leave-one-sampling-day-out within the outer training site. Fit feature and target scaling inside each inner fit, never on the validation day. Select hyperparameters using mean daily MSE (each day weighted equally). Evaluate outer metrics sample-weighted and state the distinction. The reverse fold has only two training days, limiting model selection stability. These data have already informed project development; do not claim untouched site independence.

## Fixed models

- Training mean: mandatory no-feature reference.
- Ridge: alpha in [0.1, 1, 10, 100], standardized predictors and target(s).
- Random forest: 100 trees, max_depth 3, min_samples_leaf in [2, 5], all features, fixed seed, n_jobs=1.
- Gaussian process: fixed RBF length scale in [0.5, 2], WhiteKernel noise in [0.05, 0.2], all kernel parameters fixed; standardized predictors and target(s), no optimizer. Predict each forward band separately.

Do not add neural networks or expand the search based on outer results: fewer than 65 rows does not justify a large tuning search. Save every inner score, chosen configuration and training-only out-of-range flags. No uncertainty-calibration or predictive-density claim is made for RF/GPR here; GPR posterior standard deviations are not substituted for externally calibrated errors.

## Comparison and uncertainty

Join existing physics positive-count point predictions by dataset/sample, never row order. Verify exact matched targets/sample sets. Physics medians and ML squared-error estimators differ, so describe this as comparison of reported point estimates rather than a complete proper-score comparison.

Compare RMSE with paired sampling-day bootstrap of whole blocks. Fewer than five days: suppress confidence intervals and state insufficient blocks. With five or more days report intervals as conditional sampling uncertainty, not generalization across independent sites or full model uncertainty. No winner declaration, hypothesis test or melt claim.

Do not compare four-band ML outputs to physics broadband-HCRF metrics. They are different endpoints. A future full predictive-density comparison requires an explicitly validated covariance and count-observation model; v1 does not invent one.

## Reproducibility and resources

Run locally with one numerical-library thread and one tree worker. No BioSNICAR, PySCF, GPU or network assets are needed. Record code commit/hash, Python/library versions, raw empirical CSV hashes, converted data hash, runtime and available RSS measurement. Output directories must be new; don't overwrite reports. Resume each fold/task/model only if its saved input/code signature matches; every completed bundle stores predictions and inner scores.

References: [scikit-learn grouped cross-validation and leakage guidance](https://scikit-learn.org/stable/modules/cross_validation.html), [Gaussian process implementation](https://scikit-learn.org/stable/modules/gaussian_process.html).
