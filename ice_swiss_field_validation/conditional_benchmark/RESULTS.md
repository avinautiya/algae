# Conditional Swiss field comparison — result

**The new fitted flow-sensitive thermal predictor does not establish consistent superiority over established heat-transfer models.** It improves the nominal trained-constant-h baseline at both stations, but that improvement fails in some temperature scenarios at P5; smooth Gnielinski performs substantially better at P5.

## What ran

Original source CSVs, including direct bed elevation and epoch uncertainties, were recovered and CRC/SHA verified. This supersedes the earlier missing-CSV limitation. Two analysis freezes and both outcomes are preserved:

1. **Strict native support:** every one of the20 adjacent manual intervals has some missing quality-approved native temperature support. Consequently there are no fully supported training/test intervals, and the strict run produces no fitted field scores. Partial native integrals were not divided by coverage or compared as if they represented the complete interval.
2. **Secondary author-daily-mean scenario:** author sample-weighted corrected daily means are imposed constant during each calendar day, with missing days refused. This is an explicit forcing reconstruction; it is not a complete measured heat integral. Train on8 intervals ending by the25July visits; test on9 later complete-support intervals (P3=5, P5=4). Three intervals are excluded because temperature or derived-Re support is absent. No observations are removed because they have negative or small bed change. The targets and paper were already inspected, so this is retrospective temporal evaluation, not blind validation.

All27 declared combinations of source-clock shift (−2/0/+2h), temperature shift (−0.05/0/+0.05°C) and transfer length (1.08/1.30/1.52m) ran. None was selected by test performance. The fitted models have one pooled nonnegative scalar each and exactly the same training intervals. Nominal fitted h=3884.56W/(m²K); DB0.4 multiplier=2.5852; constant rate=3.7592e−6m/s. Values are effective fit parameters, not independently measured material constants.

## Nominal held-later interval MAE in metres

| Comparator | P3 (5 intervals) | P5 (4 intervals) | Pooled (9 intervals) |
|---|---:|---:|---:|
| Training-only constant melt rate |1.34793|1.49981|1.41543|
| Training-only constant h |0.17356|0.50955|0.32289|
| Standard DB0.4, unfitted |0.49161|0.32061|0.41561|
| Cooling DB0.3, unfitted |0.54578|0.40812|0.48460|
| Smooth Gnielinski, unfitted |0.38659|0.16476|0.28800|
| Training-only scaled DB0.4 |0.13787|0.43781|0.27118|

Nominal scaled DB reduces MAE relative to constant h by20.6% at P3 and14.1% at P5, and16.0% pooled. This supports a limited benefit from flow-sensitive heat exchange **under this imposed forcing reconstruction**, not unique novel physics. Against smooth Gnielinski, scaled DB improves P3 but worsens P5. Pooling hides this reversal and cannot justify a consistent-superiority claim.

## Engineering sensitivity and dependence

Scaled DB beats constant h in all9 distinct clock/temperature cases at P3 but only6/9 at P5. It beats smooth Gnielinski at P3 in all27 scenarios and at P5 in0/27. The length scenarios are not independent replications: fitted DB predictions are invariant to a global transfer length because its fitted coefficient absorbs the1/D factor. Thus27 scenario counts must not be presented as27 independent demonstrations. There is no confidence interval or significance calculation.

The original near-bed temperature is not a demonstrated bulk-water temperature. Daily means have sampling-frequency bias; daily Reynolds numbers are author-derived hydraulic inputs including reconstructed flow. A fixed-Re transfer-length scenario is not a dynamically consistent change in channel geometry. Source UTC headers and narrative CEST require continued clock reconciliation. The shifted-clock runs disclose sensitivity rather than resolving the source convention. Direct bed changes remain in a reference tied to the glacier and do not independently isolate every mechanical/ice-motion contribution.

## Reproducibility and next scientific decision

`PROTOCOL_V2_ORIGINAL_DATA.md`, `FREEZE_V2.json`, `SECONDARY_DAILY_SCENARIO.md`, and `FREEZE_DAILY_SCENARIO.json` precede the respective scores. `run_comparison.py` runs all models. `standard_closures.py` has four passing checks. `audit_comparison.py` independently recomputes saved errors and checks training membership and transfer-length invariance. The JSON preserves every prediction, calibration parameter, exclusion and scenario. Existing author-fitted Nu and observed-friction products are excluded from the predictors.

The study now identifies a quantitative failure to address: a single shared flow-sensitive gain overpredicts P5 while underpredicting P3, and its apparent benefit depends on forcing reconstruction. Additional model structure needs independent physical justification and a new evaluation split; fitting station-specific gains to these test errors would not establish transfer. Whole coupled cliff-melt/contact survival and general Himalayan applicability remain unvalidated by this channel-component test.
