# Native Pirámide process tables decoded

The MATLAB MCOS blocker is resolved **for native process-array access**. The public [mat-io parser](https://github.com/foreverallama/matio), version 0.4.2, runs in the isolated `venv312/` with existing Python 3.12. MATLAB and Octave executables are absent; no installation or license is needed for this conversion. Python 3.9 could not meet parser dependencies and was not modified.

`decode_piramide.py` loads only debris temperature, moisture, air temperature/humidity and wind; no camera ablation or calibration file is read. Source: [Zenodo 11581073](https://doi.org/10.5281/zenodo.11581073). All twelve native process tables are exported to `decoded/`. Every numeric column is checked against the raw MCOS `data` property, names against `varNames`, dimensions against `numRows/numVars`, and all converted timestamps against the raw datetime epoch-millisecond property. No column order, timestamps or timezone is inferred from a temperature curve. This is an internal raw-property/converted-object check, **not a comparison against an independent MATLAB export**. The parser's undocumented format implementation remains a dependency, pinned in `requirements_decode.lock`.

## Actual source-array findings

| Family | Tower1 rows | Tower2 rows | Tower3 rows |
|---|---:|---:|---:|
| Debris temperature | 4760 | 4719 | 4907 |
| Moisture | 4709 | 4712 | 4715 |
| Air T/RH/dewpoint | 12422 | 10232 | 1568 |
| Wind | 12422 | 10232 | 1568 |

Debris temperature and moisture tables have uniformly 1800-second native cadence and no missing numeric values. Although the paper describes five-minute moisture acquisition, the released process table is half-hourly; never pretend native five-minute observations were exported. Air/wind are not uniformly half-hourly: their actual steps/gaps are enumerated in the manifest, and require an explicit aggregation/coverage operator.

Tower1 temperature labels run 2022-11-24 08:00 through 2023-03-03 11:30; other native ends differ from paper campaign summaries. Native timezone fields are empty. Preserve local labels, don't append Z/UTC. Same-source relative thermal analysis can use common-clock provenance, but cross-source/camera time alignment needs an independent clock check before scoring.

Native DT5 is not fixed at zero: Tower1 reaches 4.384°C (3 February21:00), with 2899/4760 records above0.4°C; Tower2 reaches2.152°C. Tower3 starts with all depths warm (DT5 reaches35.409°C on24November12:00). These are feature-only diagnostics, not melt outcomes. They flag unresolved effective basal geometry/setup states; neither last sensor nor nominal thickness proves a melting interface. Do not clip temperatures, force last sensor to0, or turn its gradient into field-validated melt. Some moisture readings are negative, and Tower1 SM2 exceeds the source-calibration porosity0.44. No corrections or clipping were made; physical-QC and calibration interpretation are required before moisture-based heat capacity.

## Final-paper metadata and model restrictions

The final [Melo-Velasco et al.2025 paper](https://doi.org/10.1029/2025JF008360), table1/sections2.2–2.3, supersedes the preprint. It specifies ±2cm thickness/sensor-depth uncertainty. The nominal0cm sensor is covered by approximately2cm gravel; its temperature is a near-surface boundary, not radiometric surface temperature. Tower1 stayed upright; Tower2/3 fell12February/7December. In-debris instruments are separated from the towers. Homogeneous/two-layer profile inference and melt-optimized conductivity already appear in the source work; they are not new methods. Thermal sensor accuracy is±0.4°C. Camera uncertainty is not established in this audit, and its files remain unopened.

`PROCESS_BINDING_V2.json` provides source geometry and decoded rock density/specific heat. Rock heat capacity is not bulk debris heat capacity; source calibration porosity0.44 is assumed, and wet heat storage/advection is not independently identified. Diffusivity from temperature can be trained without melt leakage, but conductivity/melt magnitude still needs an independently justified heat-capacity treatment.

## Reproduce and consume

```
/opt/anaconda3/bin/python3 -m venv venv312
venv312/bin/python -m pip install -r requirements_decode.lock
venv312/bin/python decode_piramide.py DATA_Piramide.mat decoded
venv312/bin/python -m unittest discover -s . -p 'test_*.py' -v
```

Four consumer tests pass, covering expected native dimensions, original labels/clock, tampered-CSV refusal, and the basal-phase warning. `read_verified_process.py` verifies CSV hashes before use. Native source/CSV fingerprints and every full-array check are in `decoded/DECODE_VERIFICATION.json`. Do not package virtual environments; package scripts, lockfile, source MAT, binding, decoded CSVs and verification JSON.

A next thermal-column experiment can now run with authentic observations, but its output is temperature/basal-heat prediction conditional on a defensible boundary model. Melt-field superiority and originality remain unestablished; the source-geometry/phase/heat-capacity gates cannot be bypassed by successful file decoding.
