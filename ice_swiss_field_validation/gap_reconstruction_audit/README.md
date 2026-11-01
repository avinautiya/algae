# Controlled public-data temperature-gap reconstruction

This bounded audit measures reconstruction of withheld **temperature time integrals**, not glacier melt, hydraulic contact, or a new independent field campaign. The selector and all three methods were fixed before scoring in `PROTOCOL_FROZEN.json` (SHA256 `1a1c5534e1733db621c3912c80dd7fcca001924b2d0ab6c7c905b06ce1545128`). This is a local pre-execution freeze, not externally registered preregistration.

The window is 16 July–8 August 2019 inclusive. For each station and duration (1/6/24 hours), consecutive calendar blocks are generated from midnight 16 July, advancing by exactly the block duration. Only blocks whose native measured-temperature reconstruction has 100% support and finite supported endpoints are retained. No block was selected using its error. Block truths are exact piecewise-linear native-time integrals; all interior values are withheld from the three endpoint reconstructions.

## Mean absolute block-mean-temperature errors

| Station | Duration | Retained blocks | Linear bridge | Previous endpoint | Next endpoint |
|---|---:|---:|---:|---:|---:|
| P3 | 1 h | 534 | 0.00116°C | 0.00365°C | 0.00342°C |
| P3 | 6 h | 83 | 0.00412°C | 0.01581°C | 0.01480°C |
| P3 | 24 h | 17 | 0.01450°C | 0.01678°C | 0.01738°C |
| P5 | 1 h | 542 | 0.00146°C | 0.00629°C | 0.00622°C |
| P5 | 6 h | 86 | 0.00961°C | 0.03022°C | 0.03336°C |
| P5 | 24 h | 18 | 0.03533°C | 0.04219°C | 0.04048°C |

All three methods and all durations are retained. Linear bridging has lower mean absolute integrated errors in these six predefined station/duration groups. That is evidence for this **retrospective temperature-reconstruction component** on observed support; it is not consistent superiority of the complete melt predictor. At 24 hours the bridge has negative signed bias (P3 −0.01176°C; P5 −0.03259°C), plausibly reflecting the selected midnight endpoints and the diurnal cycle. No phase-shift selector was tried after seeing this result.

Full MAE in integrated °C·s, bias, RMSE, median and maximum absolute errors, candidate counts and excluded-block start times are retained in `RESULTS.json`; every block prediction is in `MASKED_BLOCK_PREDICTIONS.csv`. These errors refer to the supported native floor-temperature proxy, not true bulk temperature. The source ±0.05°C instrument-accuracy bound and unknown floor-to-bulk offset are not removed by this masking exercise.

## Meaning for the actual missing data

The actual longest P5 gap in the two-week incision interval is about 4.92 hours, while P3's is about 23.44 hours. The six-hour and 24-hour masking groups are relevant stress tests. They do not give a proven upper bound or a confidence interval for actual gaps: actual missingness can coincide with sensor handling, calibration, immersion or changing channel conditions. P5 six-hour bridge maximum error was 0.03270°C and P3 24-hour maximum was 0.03712°C, but neither is a guaranteed bound on an unseen gap. Reconstructing gaps still requires a declared sensitivity assumption and a separately labelled result.

Blocks are nonoverlapping within each station/duration, but durations overlap each other and observations share station, environment and neighboring endpoints. No independent-block confidence interval is asserted. Linear bridging and the next-endpoint method use future observations and are not causal forecasting methods. The audit tests only interpolation between observed endpoints; it does not assess extending temperatures beyond sensor deployment.

## Reproduction

```sh
python -m unittest discover -s . -p 'test_*.py'
python run_audit.py /path/to/original_archive_subset --out new_run
python make_figure.py new_run/RESULTS.json --out new_run
```

The source subset must retain the verified `SOURCE_MANIFEST.json`. The runner verifies its hash and every native temperature input hash. It imports the sibling `temperature_digitization/native_temperature.py`; no incision/melt outcomes or fitted heat-transfer coefficients are read. Three tests establish exact constant/linear-signal reconstruction and exclusion of incomplete support. No forecast or melt model is retrained here.
