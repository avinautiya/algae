# Independent V3 field-comparison audit

**Arithmetic and freeze integrity pass; the candidate does not beat every baseline.** The coverage-amended comparison retains all39 models and its unchanged primary. Both independent audits pass:800 numerical/source/metric checks and the coverage/parameter-prefix audit.

The original V1 physical no-go and original V2 unscorable endpoint remain preserved. V3 replaced the nonfinite2022-12-30T03:00 endpoint with the nearest finite native endpoint2022-12-30T00:00,3h earlier, under a published availability-only amendment. New final predictions were published before that finite endpoint was released. This is an amended one-block conditional hindcast, not unchanged execution of V2.

## Actual locked result

Target block:2022-12-16T03:00 to2022-12-30T00:00, source local-naive clock. Observed native ice ablation: **57.03636361mmw.e.** The fixed primary predicts **61.39355129mmw.e.**, overpredicting by **4.35718768mmw.e.(7.639%)**.

|Model|Prediction mmw.e.|Absolute error mmw.e.|
|---|---:|---:|
|forward__middle_gap_20mm__left__proxy|61.393551|4.357188|
|derivative_OLS__middle_gap_20mm__left__proxy|63.500942|6.464579|
|constant_rate|84.440916|27.404553|
|positive_buried_T_DT1|63.454591|6.418228|
|positive_air_TA3_linear_gapfill18h|59.048301|2.011937|

The primary wins nominally against derivative OLS, constant rate and buried-temperature PDD, but loses to the explicitly gap-reconstructed air-temperature PDD. The OLS implementation is not the original source NTLS routine or a full professional energy-balance model. Air forcing includes predeclared linear gap filling up to18h, so its support limitation remains visible; its better score is nevertheless retained.

## Geometry/capacity family

All18 forward geometry×hold×capacity scenarios beat their matched OLS comparator and the constant/buried-temperature baselines. **All18 lose to the air-temperature baseline;0of18 beat all four.** Forward errors range4.283467–4.397391mmw.e.; OLS errors range5.888689–6.553031. These are deterministic finite-scenario counts, not probabilities, independent replications or continuous uncertainty bounds.

## What measurement precision would identify each ranking?

These are conditional sensitivity thresholds, not achieved instrument errors, standard deviations or confidence intervals. Under the shared-three-reading input/recalibration sensitivity, each cumulative camera reading would need an independently justified absolute error strictly below:

|Comparison|Per-reading bound in metres ice|Nominal direction|
|---|---:|---|
|derivative_OLS__middle_gap_20mm__left__proxy|0.00170250|candidate_nominal_win|
|constant_rate|0.00466594|candidate_nominal_win|
|positive_buried_T_DT1|0.00169547|candidate_nominal_win|
|positive_air_TA3_linear_gapfill18h|0.00101734|baseline_nominal_win|

Thus the small wins over OLS/buried-temperature baselines require roughly1.7mm-ice per-reading precision to be insensitive to the declared error box. Source precision/reset semantics are unknown, so this robustness is **not established**. Paint-mark spacing is not an estimate of tracking precision. The auditor separately records the already-issued fixed-forecast thresholds; those differ because its training calibration is no longer recomputed.

## Verified integrity and remaining scope

The coverage audit finds identical physical parameters, source hashes, all training targets/doses/gains and exact unchanged prefixes of the raw physical traces. No model was reselected after scoring. Raw latent/sensible ledgers remain attached to the unscaled physical component; issued predictions include an empirical training gain and are not a physically conserved rescaling.

This is one site and one later calendar block. It does not establish independent-glacier generalization, statistical superiority, professional-model superiority or historical mathematical novelty. Fixed geometry, assumed initial ice temperature, finite deep boundary and omitted wet-debris phase heat remain explicit. A further fixed batch can test persistence of the predictive benefit without dropping this air-baseline failure.

Artifacts: `V3_ACTUAL_RESULT_AUDIT.json`, `V3_COVERAGE_INTEGRITY_AUDIT.json`, and the untouched V3 `melt_comparison/RESULTS.json`.
