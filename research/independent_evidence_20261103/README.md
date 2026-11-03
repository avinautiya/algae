# Professional baseline reproduction and a prospective sensor experiment

This is an additive research package. Claude's worktree, jobs and frozen results remain unchanged. It contains a conditional reproduction of one published model score, a separately frozen sensor-signal experiment, an independent identifiability counterexample and a nonperformance support audit. These are not a validated new glacier-melt model.

## Executed outcome

The frozen conditional sensor experiment has now been scored once: diffusion equal-siteMAE1.1664W/m² versus0.8465W/m² for causal ridge. The declared superiority rule FAILED. Both numerical gates passed. See `HELDOUT_RESULTS.md` and the full unchanged prediction records; no claim of melt accuracy follows.

## Professional baseline

`reproduction/` resolves the official protocol's output units (m water equivalent per step), ice density900kg/m³ and exact run indices. The actual DEB_CF Changri Nup output sums0.718645m water equivalent. A predeclared calendar sum of the raw thresholded lowering, converted using that density, gives0.921055m water equivalent. The resulting conditional final-cumulative error is−21.975864%, versus the published−23.3%. The1.324136percentage-point discrepancy is unresolved. The executable correctly exits2; no density/date choice was searched to force agreement. Final processed observations and daily aggregation remain needed. The discrepancy alone establishes neither an error in the paper nor independent accuracy.

Official protocol sources are CC-BY4.0, Working Group on Debris-covered Glaciers, DOI10.5281/zenodo.3362374. Publication includes source metadata, code and checksums. `fetch_protocol.py` restores the original documents; whole downloaded PDFs and extracted protocol text are not required in the code packet. The paper/table source is https://doi.org/10.5194/tc-20-1895-2026 .

## Prospective conditional plate-signal test

`extension/PROTOCOL.md` and `extension/FREEZE.json` fix the code, data, dates, baselines and success rule BEFORE held-out scoring. A published pre-execution commit is required by the runner. All baselines and the opposite-orientation sensitivity remain visible. The representative-point diagnostic requires explicit opt-in; native interval-flux validation is blocked.

Training is1–28December2015; evaluation is29December2015–11January2016, Nepal time. The fixed28-day/14-day windows have671/672 complete site1/site2 training records and336 complete evaluation records at each site. The availability audit has not scored held-out values. Site1 diffusivity uses an additional interior-temperature training probe; that auxiliary information is disclosed. Site2 inherits diffusivity but calibrates its own raw-signal scale, so it is not zero-shot site transfer.

The output predicts archived plate readings conditional on observed boundary temperatures. It does not predict future weather, latent melt, cliff retreat or an intrinsic thermal conductivity. Plate orientation and interval alignment are unresolved; sign is learned only on training records. Sensor migration and omitted heat processes may cause failure. The strict declared superiority rule is not a promise of a favorable result.

To reproduce the now-completed conditional diagnostic from the repository root, use a NEW output directory:

```sh
python3 research/independent_evidence_20261103/extension/plate_experiment.py research/prediction_extensions_20261103/benchmark/CHANGRI_ALIGNED_AUDIT.json --out /path/to/NEW-results --commit 2a0828dd2a4c8c42562a277eec8dfb2b991cf035 --interpretation representative_point_diagnostic
```

This reruns the disclosed test; it does not create fresh validation data. `verify_evidence.py` runs the20 synthetic/unit checks separately. Professional-score reproduction commands and unresolved source requirements are in `reproduction/README.md`.

## Identifiability and scores

`decision/` includes an executed synthetic example: two parameter pairs fit the same100W training observation but give66.43 versus36.66W under changed forcing. This diagnoses sparse calibration, not a new theorem or measured field uncertainty. `coverage/` verifies timestamps, missingness and fixed support without selecting favorable target values.

Twenty unit/synthetic tests pass. No test count certifies historical originality, model superiority, journal acceptance or award selection. If the held-out experiment fails, preserve that result and do not tune against it.
