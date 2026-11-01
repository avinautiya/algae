# External thermal prediction adapter and Rochemelon readiness

This package adds an executable finite-water-heat predictor and independent-event data readiness audit. It does not contain a successful external field validation.

The primary model parameters and comparators were fixed from the Swiss V5 analysis before numerical extraction of the external bed-change targets. They have not been refitted to Rochemelon. This is not an externally preregistered blind study: the paper and experimental design had already been examined.

## Run the actual predictor

The numerical model uses only the Python standard library:

```sh
python -m unittest discover -s model -p 'test_*.py'
python model/predict_external.py model/SYNTHETIC_PACKET.json --out synthetic_prediction.json
```

The example is synthetic and is stamped as such. It verifies execution and accounting, not environmental accuracy. All seven outputs and their fixed Swiss parameters are preserved. Six physical families allocate finite incoming water heat across supplied rectangular reaches and separate bed and sidewall melt; the uncapped statistical constant-rate reference has its own explicit nonphysical ledger status. Cold fresh-ice sensible cost is included as an approximation; no cold-ice PDE, moving geometry, atmospheric exchange, or frictional dissipation is solved. Prescribed hydraulic flow is not updated with produced meltwater. Read model/README.md and its frozen protocol before supplying environmental inputs.

Geometry, forcing clocks, interval support and inlet-to-target route must be supplied with provenance. Unknown path/geometry and missing interval support are refused. Critical-depth closure is restricted to a sourced, verified point control; it does not supply depths throughout an arbitrary channel.

## Recovered external observations

- forcing/: 55 paired graph-digitized lake-temperature/discharge samples from the primary Rochemelon 2005 article, with source pixel provenance, uncertainty distinctions and own plots. These are digitized processed observations, not downloaded raw sensor CSVs.
- target_audit/: nine labelled A/B/C profile values from three survey dates, retained as relative elevations because the raster clips absolute-altitude prefixes. Unknown survey hours and approximately 0.82–0.96m graphical bounds are preserved.

The publisher source PDF and full-page/source images are intentionally not redistributed here. The figure extractor accepts a source image explicitly and checks its SHA256; use the source URLs and retained extraction notes to obtain the original. Our numeric fact tables and our newly generated plots are included. Source PDF: https://www.cambridge.org/core/services/aop-cambridge-core/content/view/C183D075F0E9ABCD7CE5251D6096FCA4/S0022143000206357a.pdf/outburst-flood-hazard-for-glacier-dammed-lac-de-rochemelon-france.pdf . For Figure 5, Poppler pdfimages applied to PDF page 4 extracts the original embedded JPEG; the required SHA256 is recorded in FIG5_EXTRACTION_MANIFEST.json. Optional graph reproduction uses Pillow/NumPy/Matplotlib and Poppler; it is not required to run the numerical predictor.

## Current decision

The external accuracy gate remains RED. Figure 6's single bed series is not explicitly identified as A/B/C, dated profile surveys omit exact hours, and lake inlet temperature cannot silently become downstream local temperature. The source paper's event-fitted hydraulic/thermal coefficients are not independent benchmarks on that same event. No external accuracy score is reported, no target mapping was chosen by model performance, and no missing geometry was invented to obtain a win.

The original Swiss recovery and all positive/negative comparisons remain in the sibling ice_swiss_field_validation package on the parent branch. A consistent-superiority claim is not supported by the current evidence. The model is an executable conditional component predictor, not a validated full glacier forecast.
