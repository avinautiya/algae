# Kennicott source-only boundary reconstruction

`boundary.build_boundary(profile_id,start_utc,end_utc,air_station='Debris_T2',proxy_error_multiplier=0)` reads only SHA-verified meteorological and thermistor features. It never imports the stake mass-change field. Module paths resolve relative to this packet, without external local dependencies.

Primary air is `Debris_T2` at every site; `Debris_T3` is a declared sensitivity. The first thermistor is at5cm for2020b,2cm for2023a/c and9cm for2023b. This is a **buried thermistor plane**, not the physical/radiometric debris surface. A forward column must use its corresponding vertical coordinate and separate stake thickness; the helper does not move a probe into the ice or equate profile and stake geometry.

The native initial thermistor gaps are36h(TRAIN) and48h(TEST) under the nominal ADT-noon date assumption. Parameter/proxy fitting excludes the first5days of thermistor operation, then fits an affine thermistor–air relation over the following14days(337hourly nodes). A chronological first7days fit/later7days check(168/169nodes) gives temperature-only RMSE. These are autocorrelated temperature samples, not independent melt-validation replicates.

The primary hybrid uses the proxy during the initial missing interval and5-day installation period, then measured first-thermistor temperatures. The whole target interval is retained. The source's exclusion language supports disregarding disturbed initial data; choosing exactly5days is this protocol's assumption, not an exact source-prescribed value. The full `proxy_boundary_C` series supports an explicitly separate weather-conditioned model. It still fits buried-temperature behavior using observations later in the campaign, so it is not a prospective forecast from weather available at the start of the season.

Forcing values are hourly instantaneous-label resamples, with declared linear interpolation of native5-minute sensor and hourly air data. No negative temperatures are clipped. Extrapolation and interpolation gaps exceeding10min(thermistors) or2h(air) are refused. Source UTC is explicit; output appendsZ, while the source's original stake dates remain uncertain. Nominal ADT-noon is an assumption used for feature audit, not a recovered visit timestamp.

`proxy_error_multiplier` ±1 adds/subtracts2blockedRMSE during proxy filling. For the full-proxy series it applies throughout. This produces a finite sensitivity family; it is **not** a confidence interval or guarantee that the true missing temperatures fall within the envelope. Source instrument errors, seasonal nonstationarity, depth errors, initial subsurface state and phase/geometry omissions are separate uncertainties. No model is calibrated to melt here.

Output includes source hashes, regression coefficients, reconstruction error, proxy hours, the boundary plane and original parameter-fit dates. Native AIR2 blocked RMSE is1.805C(2020b),2.912C(2023a),3.432C(2023b) and2.462C(2023c). These errors must be carried forward, not described as exact reconstructed temperatures.

From the packet root:

```
python3 -m unittest kennicott_forcing.test_boundary -v
python3 -m kennicott_forcing.run_feature_audit
```

Nine tests pass, including interpolation/extrapolation gates, exact synthetic affine recovery, UTC interpretation and a source-backed check that the uncertainty shift changes only proxy epochs of the primary hybrid. `FEATURE_ONLY_BOUNDARY_AUDIT.json` binds eight nominal traces; none contains observed or predicted melt.
