# First native camera-ablation comparison

The candidate improves nominal error over three declared baselines, but loses to the air-temperature index. It does not satisfy the all-baseline superiority objective. All 39 predictions remain reported; no best scenario replaces the primary.

## Result

The source camera stake indicates **57.036364 mm water-equivalent ablation** from December16 03:00 to December30 00:00,2022, using native clock labels. All models use one gain learned only from December2–16, whose measured total is85.201645 mm w.e.

| Fixed method | Prediction, mm w.e. | Absolute error, mm w.e. |
|---|---:|---:|
| Forward thermal candidate, primary |61.393551|4.357188|
| Derivative-OLS thermal comparator |63.500942|6.464579|
| Constant rate |84.440916|27.404553|
| Positive buried-temperature index |63.454591|6.418228|
| Air-temperature index, explicit gap reconstruction |59.048301|2.011937|

The primary relative error is about7.64%. Its absolute error is about32.6% smaller than the derivative comparator, but larger than the air-temperature comparator. Across all18 forward parameter/reconstruction cases, every case beats its matched derivative comparator, constant rate and buried-temperature index; every case loses to the air-temperature index. These18 scenarios are not18 independent field validations.

## Freeze and source history

The first, melt-only model failed its own physical gate before camera observations were opened. Negative basal heat persisted under refinement. That no-go remains at `model/frozen_piramide_v1/PHYSICAL_GATE_NO_GO.md`.

The separately frozen cold-ice extension stores cooling as sensible cold content and melts only when the interface energy balance permits it. No water supply is invented for refreezing. Its raw energy ledger closes to numerical tolerance; later empirical output gains are explicitly separate and do not inherit that physical energy ledger.

The original V2 final camera value at December30 03:00 is a native NaN, so that evaluation remains unscorable. The deterministic nearest-finite rule within12h selected December30 00:00, three hours earlier, using only timestamps and finiteness. No finite replacement value or model error was inspected during selection. Every model was recomputed for that same endpoint, with exactly unchanged training gains and parameters. All39 amended predictions were published at `eb08f2b2e62629b083d0f4e0a7f0d49b3acbdba8` before the replacement observation was released. This is a transparent eligibility amendment, not an unchanged original-protocol success.

Independent audits verify all39 model arithmetic and source binding (800 checks), and unchanged training targets/gains, parameters and model IDs. The shortened physical traces are exact prefixes of the original run. Independent native endpoint checks reproduce57.036364mm w.e. The supplied ablation-column differences equal metres water equivalent; multiplying those by1000 agrees with the frozen raw exposed-stake-distance difference times917. The field operator was not switched to improve a score.

## Uncertainty and scope

Camera endpoint precision remains unknown. The published ±2cm pertains to debris/sensor geometry, not camera ablation, and cannot be borrowed as a camera confidence interval. Conditional shared-reading calculations indicate millimetre-scale precision is needed to preserve the close model rankings; this dataset has not independently established that precision. Serial periods and shared endpoints must not be treated as independent camera-frame replications.

This is one observed-temperature-conditioned tower hindcast, not a future-weather forecast, validation of the original cliff/stream solver, or superiority over full professional glacier models. The derivative OLS comparator is not a reproduction of the source NTLS or full energy-balance model. Static geometry, ice initial/deep conditions and omitted wet-debris phase/transport remain assumptions. Existing literature already models cold ice and debris heat storage; adding those terms is not a historical novelty claim.

## Next fixed test

The unchanged39 models will be assessed over all three source-eligible later calendar periods throughFebruary10, using the same December2–16 gains. The February10–24 period is excluded for every method before seeing its outcomes because air observations endFebruary21, preventing the full declared comparison. No later outcome has been used to select a model or period. The new batch's primary metric is mean absolute period error, declared before its outcomes are opened; the first failed all-baseline comparison remains separately visible and is retained in a secondary four-period summary.
