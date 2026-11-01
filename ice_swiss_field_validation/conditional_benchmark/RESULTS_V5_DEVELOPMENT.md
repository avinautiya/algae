# V5 local-geometry development reanalysis

**A physically local geometry improves its own fitted thermal family on a matched cohort, but does not produce a superior predictor.** All81 frozen development cases completed. This is a target-exposed reanalysis following V3/V4 failures, not a new held-out validation.

The sole feature change derives station-specific hydraulic diameter and Reynolds number from imposed daily depth proxies, conserved discharge and the shared rectangular width. It replaces the earlier identical P3/P5 transfer-coefficient assumption. No station-specific gain, added fitted parameter or test-based width selection is introduced. `HYDRAULIC_INPUT_QUALITY.md` retains pressure-sensor-height and reconstructed-flow caveats.

## Nominal result: linear gap reconstruction, shared2m width

| MAE (m) | P3 (4 test intervals) | P5 (2 test intervals) | Pooled |
|---|---:|---:|---:|
| Fitted constant rate |0.65849|0.28000|0.53233|
| Fitted constant h |0.32962|1.18843|0.61589|
| Unfitted DB0.4 |0.57743|0.58043|0.57843|
| Unfitted DB0.3 |0.64161|0.66168|0.64830|
| Unfitted smooth Gnielinski |0.45419|0.42017|0.44285|
| Fitted scaled DB0.4 |0.37358|1.00745|0.58487|

Local stage availability removes two P5 training intervals, leaving six training intervals (five P3, one short P5) and the same six V4 nominal test intervals. A direct comparison to the original eight-interval V4 calibration would confound the geometry change with the training cohort. `V5_COMMON_COHORT_FEATURE_AUDIT.json` refits the exact old V4 families on the same six training intervals, without changing the fit method, and scores the same six test intervals.

On this matched cohort, old global-geometry scaled DB has MAE0.41342m at P3 and1.46923m at P5; local geometry lowers those to0.37358m and1.00745m. Pooled MAE decreases0.76536→0.58487m (23.6%). That is a conditional, developmental feature effect, not overall superiority: simple constant rate and unfitted smooth Gnielinski still have lower pooled errors. At P3 the fitted constant-h baseline also beats fitted local DB.

At P5 the fitted local DB improves fitted constant-h in all27 primary-imputer engineering cases, but loses to constant rate and each unfitted comparator in all27. Left/right gap-reconstruction alternatives preserve the nominal failure; their complete results remain available. The three imputers have identical availability per scenario. No favorable width, sensor shift or clock alignment is selected.

Independent saved-output audits recompute errors to≤4.45e−16 and verify the separate rearrangement `Re=4Q/[nu(w+2h)]` to≤3.64e−16 relative. Original stage/Q input hashes match the verified archive manifest. These are numerical/provenance checks, not field skill. No further model search was undertaken after this bounded reanalysis.

The scientific conclusion is specific: physically localising geometry helps one thermal family, but calibration is inadequately constrained at P5 and the available public episode does not establish consistent superiority. The surviving predictor should not be advertised as outperforming professional glacier models or as validated coupled cliff survival. Original V2/V3/V4 failures, source gaps and all V5 scenarios are retained.
