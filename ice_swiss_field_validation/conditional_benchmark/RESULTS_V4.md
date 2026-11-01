# V4 elapsed-time native-temperature reanalysis

The result remains negative for consistent superiority: the fitted flow-sensitive predictor improves P3 but is worse at P5 than the fair constant-h baseline, the constant-rate baseline and smooth Gnielinski under nominal forcing. This reversal persists across the gap-reconstruction alternatives.

V4 corrects the temperature sampling weights using native elapsed-time integrals plus the independently temperature-tested primary linear imputer. It does not change model families, fit loss or train/test split. Left/right persistence are sensitivity cases. Endpoint-bounded gaps up to24h are reconstructed; longer gaps and extrapolation are refused. All81 prespecified engineering runs completed. The three imputers have identical case availability for each engineering scenario.

Nominal reconstruction retains eight training intervals and six later test intervals (four at P3, two at P5). Six intervals are excluded for remaining temperature or Reynolds support gaps. These are different cases from V3's nine test intervals, so comparing their original pooled errors directly would confound the forcing update with case membership. `VERSION_COMPARISON_COMMON_CASES.json` re-scores V3 saved predictions on the exact six V4 nominal test cases without refitting V3.

| Nominal MAE (m) | P3 (4 intervals) | P5 (2 intervals) | Pooled (6 intervals) |
|---|---:|---:|---:|
| Fitted constant rate |0.77335|0.23038|0.59236|
| Fitted constant h |0.21398|0.33715|0.25503|
| Unfitted DB0.4 |0.56390|0.47910|0.53563|
| Unfitted DB0.3 |0.63118|0.58357|0.61531|
| Unfitted smooth Gnielinski |0.43262|0.26800|0.37775|
| Fitted scaled DB0.4 |0.12755|0.50042|0.25184|

The fitted scaled-DB pooled MAE is only marginally below fitted constant-h, and that average hides opposite station outcomes. At P5, scaled DB loses to constant-h in all27 engineering cases for each of the three imputers (81/81), and loses to smooth Gnielinski in all81 as well. It improves P5 constant-rate only in9/27 cases per imputer. These scenario counts are dependent engineering checks, not statistical replications or confidence statements.

The six nominal test intervals have reconstructed temperature-duration fractions between0.67% and4.14%; the remainder is native supported time. Every accepted interval carries its exact native/reconstructed/unsupported duration metadata. Training intervals also include reconstructed input and its recorded fractions. This is an explicit imputed-forcing thermal-component reanalysis, not fully observed native-forcing validation.

All three independent arithmetic audits pass to at most1.78e−15; fit-membership and transfer-length-invariance checks pass. The original strict-native and author-daily results remain unchanged. Sources retain floor-point versus bulk-temperature uncertainty, clock ambiguity, author-derived hydraulic forcing and the two-site one-event dependence. No full coupled cliff-model field accuracy or environmental superiority is established.
