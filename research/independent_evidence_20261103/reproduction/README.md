# One professional baseline reproduction — conditional mismatch retained

**Executed result:** DEBCF / Changri Nup / actual10cm thickness, source-prescribed run indices3350–5748,2,399hourly outputs. The first output column sums to0.7186453837m water equivalent. The raw-calendar100daily ABL observations sum to1.0233947m surface lowering; the prescribed900kg/m³ conversion gives0.92105523m water equivalent. This yields a conditional final cumulative error of−21.975864%, versus−23.3% in published Table3. The1.324136percentage-point discrepancy exceeds the0.05-point one-decimal rounding tolerance. **The published score has NOT been reproduced.** The script exits2 and records that failure.

## New authoritative source resolution

The original experiment documents concept DOI10.5281/zenodo.3134763 resolve to version3362374 (2019-05-22). Its registered license is CC-BY4.0. Its FINAL_RUNS protocol specifies per-step output as metres water equivalent and prescribes ice density900kg/m³. The accompanying `MMIC_model_run_dates.xlsx` fixes CN start2016-04-16T01:00 and end2016-07-24T23:00, with the exact one-based forcing indices. Both forcing rows are checked against these dates before scoring. The output archive specifies column1 as the actual site thickness; no sensitivity column is selected.

Sources:
* https://doi.org/10.5281/zenodo.3362374
* https://tc.copernicus.org/articles/20/1895/2026/tc-20-1895-2026-t03.xlsx
* https://doi.org/10.5281/zenodo.3048780 (original observations)
* https://doi.org/10.5281/zenodo.15754456 (submitted outputs)

This resolves the previously unknown nominal output unit, prescribed density and run-index range. It does NOT certify the exact processed observation total used in the final2026 analysis. Source daily ABL is thresholded surface lowering from daily SR50 statistics, not an independently verified pure-ice melt observation. Exact daily interval aggregation, mixed snow treatment and final processed observations remain unresolved. The finite hourly interval and calendar-day observation support are not asserted identical merely from row counts.

The calculation used the direct inclusive calendar ABL rule declared by its code, not a search over density, offset, nearby date windows, snow corrections or output columns. The published score was read from the original XLSX, not transcribed by approximation. We did not switch to another interpretation to obtain the desired value. The discrepancy does not demonstrate an error in the paper; it identifies missing processing information for exact reproduction.

## Execute

Fetch the small original files with `python3 fetch_protocol.py --out NEW_protocol_directory`. Each Zenodo file is checked against its published MD5; all downloaded sizes are capped. The full1.22GB model archive is not required: use the CRC-verified selected output members already in the recovery-review packet. The4.3MB Changri archive is already in the prediction-extensions packet.

```sh
python3 reproduce.py --protocol-dir PROTOCOL_DIRECTORY --source-data-zip changri_nup.zip --selected-outputs SELECTED_OUTPUTS_DIRECTORY --out NEW_RESULT.json
python3 -m unittest -v test_reproduce
```

Five tests pass, including an independent synthetic end-to-end source/units fixture and preservation of the published Unicode minus sign. A synthetic numeric match still leaves observation alignment and independent validation false. Real output and input hashes are retained in `CONDITIONAL_REPRODUCTION_V2.json` and `SOURCE_MANIFEST.json`.

For a compact code publication, use `PUBLICATION_FILES.json`; original PDFs and full extracted text can remain local because the fetcher restores checksum-verified originals. No source is fabricated, no field-improvement score is claimed, and no new model is trained here.
