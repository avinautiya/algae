# Executed supplementary comparison v1

Outputs: `records/ml_comparison/`. All fixed candidates were run; none was promoted as a validated production winner. This is a retrospective positive-count ground-HCRF study, not the project's frozen all-count experiment.

## Inverse abundance prediction

RMSE in log10 cells/mL (dex), positive counts only:

| Model | Train S6 → test southern Greenland (18 samples / 2 days) | Reverse (41 samples / 8 days) |
|---|---:|---:|
| Training mean | 0.708 | 0.879 |
| Ridge | 0.451 | 0.572 |
| Random forest | 0.437 | 0.640 |
| Gaussian process | 0.459 | 0.799 |
| Existing TD-DFT D reported median | 0.447 | 0.649 |
| Existing Tier A reported median | 0.394 | 0.937 |

Physics rows are computed from existing saved predictions on exactly the same positive-count samples. They are medians, while ML predictors target squared error; this is not a full posterior/proper-score comparison.

No ML model dominates both transfer directions. In the primary direction RF is only 0.010 dex better than the TD-DFT median and remains worse than Tier A. No confidence interval is reported with only two days. In the reverse direction ridge is 0.077 dex better than TD-DFT, but the paired day-block interval is approximately −0.033 to +0.170 dex; a superiority claim is unsupported. Its improvement over Tier A is 0.365 dex (conditional day-block interval 0.322–0.407).

Inverse ML bias is about −0.33 dex on southern Greenland, and +0.39–0.52 dex on S6: transfer bias remains substantial. Twelve of 41 reverse-fold samples have at least one band outside the training min/max range. Do not hide those samples or call the models broadly transferable.

## Forward four-band HCRF prediction

Pooled sample/band RMSE in dimensionless HCRF:

| Model | Southern Greenland | S6 |
|---|---:|---:|
| Training mean | 0.126 | 0.178 |
| Ridge | 0.072 | 0.115 |
| Random forest | 0.064 | 0.129 |
| Gaussian process | 0.069 | 0.137 |

These are useful numerical baselines for the next controlled physics test. Existing physics records contain forward log density and broadband HCRF summaries, not matched per-band mean predictions, so this study **does not establish lower four-band RMSE than the physics models**. Different endpoints must not be compared as if they were identical.

Thirteen of 41 reverse-fold forward inputs fall outside the training log-abundance range. No hemispherical-albedo, satellite-footprint, absorbed-energy or melt improvement is demonstrated.

## Verification and handoff

- Nine regression tests passed.
- Run uses one numerical thread/tree worker, no BioSNICAR/PySCF, approximately 7 seconds and 166 MiB maximum RSS on the recorded macOS environment.
- Full predictions, all inner-day tuning scores, comparison intervals, exclusions and source/dependency hashes are saved.
- Five zero-count samples are explicitly inventoried and excluded from this conditional log-target analysis. Existing all-count Poisson scoring is unchanged.
- Few-day confidence intervals are suppressed rather than showing a degenerate bootstrap interval.

Claude should first obtain matched per-band physics point predictions and compare them to these forward outputs on the same rows. This additive branch does not fix or endorse the existing physics assumptions and does not require interrupting active calculations.
