# V4 native-time forcing reconstruction: versioned retrospective reanalysis

Frozen before V4 thermal scores. Earlier direct targets and V2/V3 scores are exposed. This is a documented forcing-method correction and sensitivity analysis, not a new blind test. Preserve strict-native and author-daily-mean outputs unchanged.

## Independently selected forcing reconstruction

The temperature-only masked-data experiment compares reconstruction methods without using measured incision as a response. Its fixed primary choice is linear interpolation between the finite boundary temperatures of gaps no longer than24h. Left- and right-end persistence are sensitivity alternatives, not candidate winners selected by thermal prediction accuracy. Preserve the complete masked-temperature experiment and its outcome as provenance.

Use original native P3/P5 corrected-source temperatures with final export offsets exactly once. Native valid segments use the same piecewise-linear interpretation and60s native-gap threshold as the previous strict analysis. Any gap reconstruction must have finite observed endpoints on both sides, must span at most86400s, and must not extrapolate beyond the data. Longer gaps and missing endpoints remain unsupported. The endpoints' full separation defines the24h cap; do not turn a longer gap into several fillable daily fragments. Daily labels/sample means are not used to fill native gaps.

Temperature integration uses elapsed time and exact integration of max(T+engineering_shift,0). For every manual interval report native supported duration, reconstructed duration/fraction, remaining unsupported duration, maximum reconstructed gap and method. Reconstructed temperature is an assumption supported by an input-only mask test; it is not an independently observed heat budget. Original floor-point versus bulk-water ambiguity remains.

## Unchanged thermal experiment

Use exactly the six model families and temporal train/test split frozen in PROTOCOL_V2_ORIGINAL_DATA.md: fitted pooled constant melt rate; fitted pooled constant h; unfitted DB0.4; unfitted DB0.3; unfitted smooth Gnielinski; fitted pooled DB0.4 scale. Each fitted family has one scalar fitted only on training lowerings. No change in loss, station pooling, parameter bounds, physical constants or derived-Re treatment. No author-fitted Nu or observed-friction product is added.

Run all27 clock/temperature/diameter cases for each of the three gap reconstructions (81cases total). Clock shifts −2/0/+2h move the manual-epoch windows relative to native forcing. Temperature shifts −0.05/0/+0.05°C are engineering sensor scenarios. Transfer lengths1.08/1.30/1.52m are fixed-Re length scenarios, not full dynamically consistent geometry simulations. Fitted DB's length dependence is absorbed by its gain; do not count duplicate cases as independent evidence.

Within each forcing scenario all six models share the same accepted intervals. Missing derived Reynolds support or remaining temperature gaps reject an interval for every model. Check that the three reconstruction methods have identical availability for each clock/temperature/diameter combination; otherwise compare their common cases only and disclose every rejection. Retain negative/small bed changes. No method, site, time range, forcing case or baseline is selected after V4 scores.

Report all nominal per-station and pooled errors and all sensitivity cases, with signed bias and observed/predicted lowerings. Compare scaled DB against both fair trained simple baselines and every prespecified unfitted comparator. A station-specific reversal or a failed sensitivity scenario remains visible. Input-filling uncertainty is not a confidence interval; two locations in one stream do not establish general environmental consistency. V4 can inform whether correcting the temperature time weighting changes this conditional component test, not prove full coupled cliff-system field accuracy.
