# Original products recovered: these supersede graph digitization

The public author-disclosed ETH mirror supplied original CSV products after the initial notebook-plot recovery. Source archive selection, member CRC32, lengths and SHA256 values are retained in `SOURCE_MANIFEST.json` at `/private/tmp/ogier-author-original-data-20261101/`. Whole-archive checksum verification was not performed. The digitized curves remain a provenance crosscheck only, not primary model forcing.

`original_temperature_adapter.py` verifies each input against that manifest, normalizes P3/P5 daily and hourly products, retains original timestamps, and records missing values without imputation. No timezone conversion is inferred: temperature headers omit a zone, while bed elevation headers specify UTC and paper figures use CEST. Reconcile this upstream or explicitly use both clock-alignment scenarios.

## Correction semantics checked against numerical originals

Author notebook export writes `temp3_corr_raw.csv` and `temp5_corr_raw.csv` **without** final +0.035°C/+0.040°C offsets. It writes daily and hourly products with those final offsets. Original raw-suffix arithmetic means plus those offsets reproduce original daily means to at most 0.000511°C. The difference is compatible with source rounding and helper boundary-sample omissions. Apply the final offset exactly once to raw-suffix inputs; never add it again to daily/hourly originals.

The 30 digitized July 16–30 values agree with original daily products within their previously declared graphical bounds (maximum absolute discrepancy 0.01559°C). They are not 30 additional independent field observations.

## Sampling-weight issue and native time integration

The author's helper computes an arithmetic average of finite samples. Sampling cadence varies between 5, 30 and 60 seconds, notably during tracer experiments. Therefore calendar-day labels do not make the product a physical 24-hour time average. No new calibration is required to correct the integration measure.

`native_temperature.py` applies the documented raw-suffix offset once and integrates piecewise-linear temperature between **adjacent finite** timestamps separated by at most 60 seconds. It reports supported duration, fraction of requested interval, and the longest unsupported interval. The 60-second support limit is an explicit analysis setting, not a guarantee of observability. Missing values and long gaps are not bridged. Signed temperatures remain available; `integrate_positive` is a separately labelled melt-only `max(T+shift,0)` scenario with exact zero-crossing integration.

API:

```python
from native_temperature import TemperatureSeries
series = TemperatureSeries.from_archive(archive_root, 5, max_gap_seconds=60)
result = series.integrate(start_naive_datetime, end_naive_datetime)
melt_scenario = series.integrate_positive(start_naive_datetime, end_naive_datetime,
                                          extra_offset_C=-0.05)
```

Do not call partial support a complete forcing record, normalize a partial temperature integral to a full interval without a declared gap assumption, or compare a partial-support melt total to the full observed incision as if the exposure were matched. P5 July 16/24/25 native daily support is approximately 90.8%/78.1%/80.8%; corresponding author-minus-supported-time-weighted means are −0.0314/+0.0277/−0.0086°C. This mixes sample weighting and partial-support differences, so those differences alone do not identify one causal bias.

P3/P5 instrument accuracy is ±0.05°C according to Ogier et al. §3.1.2. Sensors measure floor-point water temperature; bulk water is likely warmer and the unknown representativeness offset remains separate. A deterministic instrument-accuracy bound must not be presented as a confidence interval. The adapted forcing does not establish model superiority.

Ten tests cover raster provenance/calibration, exact linear partial integrals, irregular sampling, gap/NaN refusal, offset sensitivity, monotone clocks and exact zero-crossing integration. Source-data audit products and native hourly/daily diagnostic products are kept separately; no model fitting or field prediction scoring is performed in this folder.
