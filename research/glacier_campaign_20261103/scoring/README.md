# Explicit-window contact and heat endpoint scorer

This independent scorer prevents maintained flow from being counted as reconnection and distinguishes temporary rewetting from post-pulse persistence. It does not run a model, invent missing records, identify natural processes or certify field accuracy.

Input trace: an array of exact records `{"time_s": ..., "toe_fraction": ..., "delivered_heat_W": ...}`. Time is seconds, fraction is0–1, heat is signed power in watts. Specify which physical heat port this is; do not mix water heat entering a route with heat actually delivered to ice. Signed heat is retained, never clipped. Native source units must be converted explicitly by an audited adapter; days are not automatically guessed from values.

Specification: exact keys `pre_start_s`, `pulse_start_s`, `pulse_end_s`, `post_end_s`, `contact_threshold`, `maximum_sampling_gap_s`. All windows and thresholds must be fixed before scoring a challenge. The trace must cover the complete window, have strictly increasing timestamps and no interval larger than the declared sampling gap. Interpolation is piecewise linear; this does not prove events finer than the sampling interval are resolved. Native physical-time traces and quasi-steady morphological-time traces must remain separate evidence classes.

Outputs distinguish maintained contact, initial disconnection, lost-without-reconnection, temporary/incomplete reconnection, reconnection followed by reloss and recovery throughout the entire post-pulse window. Raw event times, wet/dry durations, signed heat, mean flow fraction, initial wetness and genuine-loss flags remain available. Whole-window sustained contact is an engineering endpoint, not infinite-time recovery. A short dry interval matters unless it falls only within floating-point duration roundoff. This is not a statistical confidence interval.

The sampled preloss trace has a stable SHA256 to detect exact duplicate prefixes. Identical hashes show identical sampled records only—not independent flow states or independent experimental replicates. Different hashes do not prove independence either. Source configuration/state/history provenance must also be grouped by the campaign.

```sh
python3 -m unittest -v test_score_trace
python3 score_trace.py trace.json spec.json --out NEW-score.json
```

The CLI refuses overwrite, nonfinite JSON and duplicate keys. Ten independent analytic-fixture tests pass, covering exact threshold crossing, signed heat integration, duplicates, maintained/reconnected distinctions, reloss, late recovery and malformed data. This code is additive and does not change Claude's work or historical outputs. The existing fixed-ROI melt scorer was inspected separately; contact endpoints do not replace matched-footprint melt comparisons.
