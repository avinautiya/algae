# Causal thermal-graph melt inference: algorithm development and source gates

This is an additive development campaign following the failed locked Tower2 moisture test. It implements four distinct physical algorithm families and a transferability certificate; it does not repeat the rejected moisture-coefficient search. All Tower1/Tower2 targets used here had already been inspected and are DEVELOPMENT. Tower3 native melt outcomes remain unexamined. No failed historical result was overwritten.

## Main algorithm implemented

`adaptive/adaptive_flux.py` estimates positive layer conductances from rolling past-only integrated temperature balances, computes signed heat at a conditional ice boundary, repays cold debt, and calibrates a common heat scale with exactly one local TRAIN melt interval. Primary window7days;3/14days remain sensitivity cases. This is a particular proposed algorithm assembled from established inverse-heat, passivity and calibration ideas; historical originality and field superiority are not established by implementing it.

On the second stake it lowers development MAE from7.265 to5.339mm water equivalent relative to air-only exposure, a26.5% decrease. At the first stake it raises MAE from7.906 to25.805. It loses to the stronger calendar baseline at both. The narrow second-stake gain must not be described as general or independently validated improvement.

| Primary development method | Tower1 MAE | Tower2 MAE |
|---|---:|---:|
| Causal7day thermal graph |25.805|5.339|
| Signed storage-corrected gradient |16.479|7.676|
| Increasing-resistance surface-temperature model |8.174|9.672|
| Radiation/reservoir, reference1day |17.256|13.456|
| Air-positive-degree baseline |7.906|7.265|
| Historical calendar baseline |2.253|3.435|

All units aremm water equivalent per14day period. The historical calendar parameter used all fourTower1development observations, so its Tower1score is fitted. A separate matched leave-one-period-out calendar comparison remains strong (pooled3.663versus15.721for the radiation family); it is still retrospective, not an independent trial. Four serial periods per stake share endpoints and weather.

## Other concrete work

- `algorithm/`: signed P1 basal heat-storage correction, all24capacity/gap/site cases, exact energy tests and phase-conversion diagnostics. The storage term is too small to fix the field error.
- `geometric/`: exact moving-resistance Stefan-dose calculation and classical energy/large-depth limits. Geometry growth requires an independently justified fixed thermal boundary; settlement can invalidate it.
- `alternative/`: all four radiation/cooling-reservoir scenarios, traces, ledger/refinement checks, fair retrospective baselines and two inspected figures.
- `turbulence/`: newly recovered site-pressure/radiation source members, portable archive index, CRC/SHA/source-clock audits, four-height neutral sensible/latent heat features and tests. H/LE are conditional on flux-layer assumptions; substantial gaps/stability/accuracy limits remain.
- `alternative/neutral_flux/`: an explicit no-extra-parameter heat-exchange hypothesis and executable full-period support gate. It refuses native evaluation because complete, physically supported turbulent forcing is absent. Nothing missing was silently set tozero.
- `novelty/`: primary-source audit, feature-only Tower3 reconstruction assessment and finite-family calibration-transfer certificate with counterexample models. It returns extrema/abstention, not artificial point certainty.
- `external_lirung/`: verified Himalayan thermal source and eligibility audit. Source clock/depth/binding defects and the blocked outcome-file download are recorded; no native ablation values were accessed.

## Reproduce

Keep these sibling directories together:

```
research/piramide_thermal_test_20261103/
research/glacier_inverse_foundations_20261103/
research/adaptive_melt_development_20261103/
```

The first contains the pinned parser/numerical dependencies; `package_paths.py` resolves sibling data. Optional overrides are `PIRAMIDE_RESEARCH_ROOT` and `GLACIER_ROUND3_ROOT`. Run `python verify.py` in this directory. Tests never open new outcome files. For reruns, copy a component to a fresh output directory; existing result files refuse overwrite. Download scripts require network access and verify their source metadata. The public Lirung outcome download is not retried by any test.

`PORTABILITY_CHANGES.json` records dependency-path changes after development execution. Original execution scripts are preserved in `historical_execution_sources/`; numerical equations/outputs were not changed by those path edits. Component manifests may describe their execution snapshots; the release manifest binds the complete current package. `adaptive/CORRECTIONS.md` separately records a corrected source-depth error, preserving the invalid preliminary output and excluding it from current comparison.

## Evidence ceiling and next gate

The adaptive method is a genuine executable candidate, with a narrow development gain and a serious cross-stake failure. It is not currently an overall superior, validated melt predictor. More model complexity does not turn these inspected comparisons into independent validation. Literature already contains layered debris physics, inverse thermal estimation, time-dependent properties and goal-oriented inference; the candidate contribution must be this exact method and reliable transfer, not claims that these broad ideas are new.

Before using untouched outcomes, bind independent site geometry, sensor roles, clocks, thermal/weather support and the target operator. Resolve the Lirung stake/profile association and dates without selecting by favorable rates. Native Tower3 meteorology ends6December; any later weather-driven comparison needs a separately frozen spatial-reconstruction contract. Its pre-collapse overlap error is about1.04CMAE, not an error-free forcing source.

A Nature-level contribution or STS placement is not established or guaranteed. The new code, all trial results and all unpassed gates are included so subsequent work can build on evidence rather than a renamed failed hypothesis.
