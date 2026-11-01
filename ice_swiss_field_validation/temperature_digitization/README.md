# Audited P3/P5 daily floor-water temperature digitization

This portable revision extends the preceding `ogier_daily_temperature_digitized/` work to **16–30 July 2019**, adding the first day needed by the possible two-week incision comparison. No earlier file or result is overwritten. The output contains 30 daily curve samples and a conservative subset of 23 visually identifiable vertices. Near-collinear curve samples retain a flag; they must not be counted as independently identified plotted observations.

The immutable source is the saved PNG from Ogier's `Temperature_visualization.ipynb`, zero-based cell 9, at commit `9116f98dc7fef42be82dc4de50316f4a676ca8e9`. It is not a fit. The hourly plot's legend independently identifies P3 as green and P5 as red. Source hash is enforced before reading pixels. Existing P3 +0.035°C and P5 +0.04°C processing corrections are already included; this script adds neither again.

The retained helper-function source establishes calendar-day averages labelled at the same day's noon. It excludes NaNs, can include partial days, and omits boundary samples when transitioning dates. It does not establish continuous 24-hour valid immersion. The notebook uses naive `DateTime`; timezone remains unstated. A daily noon label does not imply an instantaneous noon reading or a noon-to-noon average.

## Errors which must remain separate

- Raster bounds use ±2 neighboring horizontal columns, ±2 vertical pixels, plus 0.001°C axis allowance. They are conservative analyst bounds, not confidence intervals.
- [Ogier et al. §3.1.2](https://tc.copernicus.org/articles/15/5133/2021/) states P3/P5 instrument accuracy ±0.05°C. This source bound is included separately from graphical error. Do not divide it by the square root of the number of readings as if it were independent random noise.
- The same section explicitly identifies the observations as **point temperatures at the channel floor**, with bulk water probably warmer. The floor-to-bulk representativeness offset is unmeasured. These data cannot be called independently measured bulk temperatures in a heat-transfer correlation without that additional assumption.
- Source export rounding is 0.001°C. Calibration, immersion screening, incomplete-day and temporal/spatial representativeness issues are additional to raster precision.

## Reproduction

```sh
python extract.py --image ../temperature_cell9_output0.png --out .
python -m unittest discover -s . -p 'test_*.py'
```

Requires Pillow; Python standard library performs all remaining operations. No temporary helper-notebook path, notebook execution, network call, model fitting, or temporal gap filling is needed. `daily_averaging_source_cell14.txt` preserves the helper implementation. `selection_proof.png` overlays every extracted sample on the unmodified raster geometry (black accepted vertex, purple near-collinear sample).

The 16 July 14:00 to 30 July 14:31 incision target overlaps partial calendar days. Exact integration over that interval is **not observed** from daily averages alone. Using a piecewise calendar-day temperature reconstruction is a declared forcing approximation; sensitivity to within-day timing and the unknown bulk offset must accompany it. This extraction supplies forcing covariates, not field superiority or a completed model score.
