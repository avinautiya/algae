# Completed inverse-physics and locked Tower2 transfer campaign

Use branch `codex/glacier-inverse-foundations-20261103`. This is additive to the existing recovery, integration and CFD work; do not overwrite those worktrees, checkpoints or frozen evaluations.

Read `FINAL_TRANSFER_REPORT.md` first. The actual new Tower2 result is negative: primary moisture-anomaly MAE13.353 mm water equivalent versus air7.265 and calendar3.435. Do not claim all-baseline superiority. All nine models and all four serial periods remain in the result. Native camera conversions and independent dose/score reconstruction were checked.

## What to import

Keep both sibling research directories:

- `research/piramide_thermal_test_20261103/`
- `research/glacier_inverse_foundations_20261103/`

The latter's `paths.py` resolves the former by relative path, or `PIRAMIDE_RESEARCH_ROOT`. Install only the pinned decode dependencies if needed; no CFD installation is needed for this package. Run `python verify.py` in the latter. Native MATLAB reading uses the dependency set in the older package.

## Evidence and freezes

- `moisture/frozen_tower2_final/FREEZE.json`: active signed-development freeze, all source/code hashes and external dependency hashes.
- `NORMALIZED_PREDICTIONS.json`: all normalized predictions published before TRAIN outcome release.
- `PREDICTIONS.json`, `PREDICTION_SEAL.json`: predictions published before evaluation release, commit66cc20de0b0f9317f69cbe8d851bfd74ec489487.
- `TRAIN_RELEASE_RECEIPT.json`, `EVALUATION_RELEASE_RECEIPT.json`: role gates.
- `RESULTS.json`: all models/periods and the failed-primary verdict.
- `moisture/review/ACTUAL_TOWER2_POSTSCORE_AUDIT.json`: independent raw-feature and score reconstruction.
- `data/moisture/review/TOWER2_*SOURCE_AUDIT.json`: native target-operator checks.

Publication is an actual Git commit with verified blobs, not a local timestamp masquerading as an external freeze. Current result files contain evaluation values. The earlier unopened-outcome statements in historical freeze documents describe their pre-release state.

## Useful components

1. `physics/layered_column.py`: conservative heterogeneous heat model and explicit warm-regime scale ambiguity. Temperature matching alone does not identify melt scale.
2. `theory/weak_graph.py`: positive integrated conductance inference with conditional LP outer bounds. Deepest observed-edge heat must not be called ice melt without intervening storage/interface corrections.
3. `theory/ice_reference.py`: conditional known-cold-ice reference inversion. Unknown history or unconstrained ice reference invalidates its identification claim.
4. `moisture/flux_ice.py`: signed-flux/cold-content coupling with conservative tests. It is separate from the empirical index and does not validate it.
5. `moisture/anomaly_index.py`: exponential calibrated-dose estimator, full sensor-offset invariance and drift bounds under explicit assumptions. Other centered functions can also cancel constant offsets.
6. `moisture/experiment.py`: freeze, normalized, calibrate and score CLI. Run fresh experiments only in a new directory with a separately frozen protocol; existing results refuse overwrite.
7. `data/`: audited Himalayan thermal sources. Ngozumpa is resolved cold-debris evidence but not paired-melt validation or a complete cold-ice reference.

## Whole-family diagnostic

Read `moisture/review/TRANSFERABLE_FAILURE_EXPLANATION.md`. Observed increasing-convex ordering proves no nonnegative exponential coefficient can improve over the already-overpredicting air baseline on these four periods. This is stronger than a finite grid search, but conditional on nominal inputs; typical sensor accuracy is too broad to guarantee it for true moisture. Do not market it as new stochastic-order mathematics or a superior model. The analysis script and six synthetic tests are executable.

## Next justified work

Do not launch a large campaign or another neural model to pursue a favorable result on these now-inspected targets. The failed transfer is a real result. Finish integrating its tables/figures into the paper and distinguish verified mathematics, retrospective tests and unverified field hypotheses. Preserve the earlier failed Swiss, Piramide and recovery verdicts.

For a new physical prediction, establish an independent heat-scale constraint and compatible interface/geometry/ablation support before spending the remaining unseen data. Tower3 targets and the unused later Tower1 endpoint remain unexamined here. If proposing them as a test, predeclare the hypotheses, baselines, target operator, shared calibration and uncertainty gates in a new externally published freeze first. They are not a license to keep testing candidate variants until a positive result appears.

Do not claim historical new mathematics, first moisture-aware glacier modeling, global superiority, an advance forecast, satellite validation, a full vertical-cliff model, or guaranteed journal/competition outcomes.
