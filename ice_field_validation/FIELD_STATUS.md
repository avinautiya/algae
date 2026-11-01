# Public-data field validation status — 1 November 2026

## What has actually been tested

Four atmospheric-component baselines were calibrated or fixed using 2021 Kennicott data and tested retrospectively against later measured stake-marker emergence times. These tests do not evaluate the newly coupled water/debris/cliff model. Observed future weather drives the hindcasts; they are not forecasts issued before the weather occurred.

The principal endpoint is time to the next 10 cm horizontal stake marker. Most source increments are fixed 10 cm events, so treating their observed duration as a prescribed integration interval would condition on the outcome. The event clock instead predicts when the specified increment is reached. Annual slope converts horizontal to normal retreat, conditionally on the stake following the face fall line.

| Mean absolute event-duration error | 2022, 105 events | 2023, 71 events |
|---|---:|---:|
| 2021 constant normal rate | 8.993 h | 7.609 h |
| 2021 positive degree-day fit | 8.396 h | 7.499 h |
| Fixed energy balance | 5.862 h | 8.113 h |
| 2021 calibrated energy balance | 8.002 h | 7.814 h |

The fixed energy balance reduces pooled 2022 MAE by 34.8%, but increases pooled 2023 MAE by 6.6%. Its 2022 east-cliff error increases despite a substantial west-cliff improvement. There is no consistent superiority across years and sites. The fitted albedo reaches its allowed upper bound, another reason not to interpret the fitted parameters as independently established material properties.

2023 uses an off-glacier weather transfer estimated from paired 2021–2022 weather. A ±15% shortwave engineering sensitivity changes the fixed-model pooled 2023 MAE from 6.061 to 9.711 h. These are scenarios, not confidence bounds or a demonstrated instrument correction. Radiometer tilt is documented and unresolved; selecting a favorable radiation scaling using test errors would invalidate that comparison.

There are two 2022 and four 2023 site-year histories, including repeated sites. Adjacent events share endpoints and weather. The 176 events are not 176 independent glaciers. The 2023 secondary test was specified after inspecting 2022 results; all four melt models were retained without refitting to 2023 retreat.

## New public evidence being assembled

The Ogier et al. (2021) author repository preserves measured daily corrected water-temperature plots, tracer discharge outputs, and authored incision observations. The original ETH archive currently refuses downloads with rate limiting. Saved notebook outputs are therefore labelled as such; graph-digitized temperatures will retain pixel bounds and must not be described as downloaded raw sensor records.

Before a thermal-incision benchmark, waterline-cut changes must be separated from stage changes. Unknown stage correction, station geometry or timing cannot be silently set to zero. Independent calibration and testing need to be fixed before using any of these observations to tune a heat-transfer law.

The Langtang DEMs provide actual geometric change, but mixed motion, collapse and changing footprints prevent relabelling it as pure thermal cliff melt. Röhl's notch-formation rates likewise differ from effective thermal melt. These distinctions are part of the measurement model, not optional caveats.

## Deliverables

- `event_benchmark/`: frozen protocols, executed source snapshots, all saved predictions and scores.
- `METRIC_AUDIT.md` and JSON: independent recomputation of every saved metric; pass to roundoff.
- `make_field_figures.py`: generates both PNG and PDF figures directly from saved results, without rerunning or fitting.
- `FIELD_BENCHMARK_REVIEW.md`, `RADIATION_TILT_AUDIT.md`: forcing and observation limitations.
- `public_hydrology/`, `field_digitized/`: provenance-preserving additional public observations and matching audits.

Public source records: Petersen et al., DOI 10.1002/gdj3.70032; Kennicott archive DOI 10.5281/zenodo.14625691 (verified original archive), latest record 16389806; Ogier et al., DOI 10.5194/tc-15-5133-2021, author repository `christopheogier/PlaineMorteLake_hydro_thermodynamics.jl` at commit `9116f98dc7fef42be82dc4de50316f4a676ca8e9`.

## Present conclusion

Actual field hindcast evidence now exists for atmospheric components. Consistently high field accuracy of the complete coupled prediction system has not been established. Keep the public field tests as a separate gate from numerical conservation, mesh convergence, and idealised recovery experiments.
