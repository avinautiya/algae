# Final independent scientific evidence verdict

**The implemented candidate does not establish superior field melt prediction.** The locked three-period comparison is now independently audited, all39 models are retained, and572 source/arithmetic/operator checks pass. The negative superiority result must remain the headline accuracy verdict. Numerical conservation and a physically consistent interface are useful engineering results; they do not overturn the field comparison.

## Actual prospective-batch comparison

All predictions used the original Dec2–16 training gains and frozen temperature-derived parameters. No new calibration or state reset was permitted. The three later calendar periods were scored together, including unfavorable periods; the fourth proposed period was excluded for every model before target access because air-temperature forcing did not extend far enough. The outcome is a conditional hindcast using observed temperature forcing, not a forecast of future weather.

|Method|Primary three-period MAE(mm w.e.)|Secondary four-period MAE including first V3 result(mm w.e.)|
|---|---:|---:|
|forward__middle_gap_20mm__left__proxy|14.220367|11.754572|
|derivative_OLS__middle_gap_20mm__left__proxy|14.349710|12.378427|
|constant_rate|18.173057|20.480931|
|positive_buried_T_DT1|14.192518|12.248945|
|positive_air_TA3_linear_gapfill18h|9.870827|7.906105|

The candidate wins nominally against derivative OLS by only0.129342mmw.e. and constant rate by3.952690mmw.e. It loses to buried-temperature PDD by0.027849mmw.e. and reconstructed air-temperature PDD by4.349540mmw.e. The air baseline is explicitly gap-reconstructed and remains the best reported primary baseline; its input limitation does not justify discarding its better score. The derivative comparator is OLS, not the original source NTLS implementation or a full professional energy-balance model.

The first, coverage-amended V3 period is preserved in the secondary four-period table. Its earlier air-baseline loss has not been dropped, and that pooled secondary score does not replace the frozen three-new-period primary metric.

## The three native observed increments

|Period|Observed mm w.e.|Candidate signed error mm w.e.|Air-baseline signed error mm w.e.|
|---|---:|---:|---:|
|2022-12-30T00:00:00 → 2023-01-13T03:00:00|59.046117|-4.455769|-0.042259|
|2023-01-13T03:00:00 → 2023-01-27T03:00:00|61.451154|+17.539118|+12.724439|
|2023-01-27T03:00:00 → 2023-02-10T03:00:00|89.054069|+20.666215|+16.845784|

Both the candidate and air baseline increasingly overpredict the later periods. This pattern is evidence of imperfect transfer, not proof of one particular missing physical process. Changing geometry, capacity, target units or the gain after seeing it would be new development and would need a new unused test.

## Finite scenario robustness

Across18 forward geometry×hold×capacity scenarios, the candidate beats matched OLS in15, constant rate in18, buried-temperature PDD in1 and air-temperature PDD in0. **No scenario beats all four.** These are finite deterministic sensitivity counts, not confidence levels, probabilities or independent field samples.

## Measurement uncertainty remains unresolved

The parent's bounded-MAE helper was checked independently against the triangle-inequality expression for all four primary comparisons on a grid of assumed per-reading errors. Issued forecasts remain fixed, while adjacent observed increments share their endpoint errors. For common per-reading boundδ(mm w.e.), the conservative change radius for a difference of mean absolute losses is4δ. This bound does not assume independent image frames and does not claim exact extrema.

At a **hypothetical1mm-ice per-reading bound**, candidate-minus-baseline MAE intervals(mm w.e.) are:

|Baseline|Conservative conditional difference interval|Consequence under that assumption|
|---|---:|---|
|derivative_OLS__middle_gap_20mm__left__proxy|[-3.797342,3.538658]|unresolved_by_conservative_bound|
|constant_rate|[-7.620690,-0.284690]|candidate_strictly_better_throughout_box|
|positive_buried_T_DT1|[-3.640151,3.695849]|unresolved_by_conservative_bound|
|positive_air_TA3_linear_gapfill18h|[0.681540,8.017540]|baseline_strictly_better_throughout_box|

The source does not establish that1mm error bound. These are conditional mathematical certificates, not measured camera precision, Gaussian confidence intervals or achieved robustness. At an assumed2mm-ice bound this conservative method identifies none of the four rankings. The0.129mm OLS margin and0.028mm buried-temperature margin are too small to support a strong improvement claim without substantially better measurement evidence.

The3 serial periods share readings, climate and one site. They are not3 independent glaciers. No frame bootstrap, significance claim or independent-glacier generalization is justified by this audit.

## Physics and conservation versus predictive accuracy

The raw cold-interface calculation closes its sensible-plus-latent ledger to at most2.06e-13 relative across the family. It is explicitly unscaled; the issued model applies an empirical training gain. For the primary that gain is about0.235 relative to the raw physical component, so the calibrated forecast must not be presented as an independently correct raw energy budget.

Raw cumulative ice-equivalent retreat spans0.324–3.088m across scenarios; the primary is about1.81m. These are model diagnostics, not observed retreat. They reinforce the need to quantify the fixed-interface/no-moving-coordinate-advection approximation, finite deep0°C boundary and uncertain initial ice state before asserting complete environmental fidelity. Wet-debris phase change is still absent.

## Defensible project contribution and stopping decision

The supported contribution is a reproducible, source-audited comparison pipeline and a verified dry/sensible-debris–cold-ice complementarity component, followed by an honestly unsuccessful superiority test. The present evidence does not establish that this component fixes professional glacier models, supplies novel environmental physics or reaches a particular publication/competition threshold.

Preserve the failed and unscorable stages, all39 predictions and all observed periods. Further changes should be driven by independently measured discrepancy/geometry/thermal-state evidence and tested on genuinely unused data. Repeatedly selecting another period, metric or material coefficient until a win appears would not establish the requested result.

Independent artifact: `ACTUAL_FUTURE_BATCH_AUDIT.json`. This audit recomputed every fixed-gain prediction and metric from released normalized doses and results; it did not reopen raw camera files or independently rerun the full thermal trajectory. The prefix/state continuity certificate is checked against frozen source metadata and unchanged original integrals, with that scope stated explicitly.
