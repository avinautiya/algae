# Changri Nup target-resolution result

**Result: the exact Table3 observed target remains unresolved. The existing data cannot be relabelled as an independent pure-ice melt target.** The declared conditional calculation was actually rerun, not merely copied:2399 model rows,100 calendar days, model0.7186453837m water equivalent, raw thresholded lowering1.0233947m. The900kg/m³ conditional conversion gives−21.975864%, against the published−23.3%.

## New source-authority distinction

The original [FINAL_RUNS protocol](https://doi.org/10.5281/zenodo.3362374) specifies900kg/m³ for conversion of **submitted model outputs**. This does not by itself establish the density or operator used in the later observation-processing pipeline.

The final [paper supplement](https://tc.copernicus.org/articles/20/1895/2026/tc-20-1895-2026-supplement.pdf), §2.2, p11, specifies915kg/m³ in the **MCC19 model description**. That is not an observation-processing instruction for all models. Applying915 to the raw lowering happens to round the resulting DEB_CF error to−23.3%. Selecting it for that reason would be outcome-based convention fitting. The executed audit retains it only as an explicitly unsupported sensitivity, never as successful reproduction.

The model total and rounded table error imply an observation total near0.936956m water equivalent. That number is an algebraic diagnostic, not a recovered field measurement, and cannot serve as validation truth.

## Distinct raw measurement operators

The verified [original archive](https://doi.org/10.5281/zenodo.3048780) yields all of the following for related but nonidentical supports:

| Source operator | Value(m) |
|---|---:|
| Sum of positive, thresholded daily ABL |1.0233947|
| Daily HEIGHT endpoint lowering |1.034822|
| Daily DIST endpoint increase |1.034820|
| Hourly DIST increase at labelled model start/end |1.044500|

These are not interchangeable. Daily means, subthreshold changes, positive height gains and hourly endpoints affect the comparison. Snow and surface displacement further prevent automatic interpretation as pure ice melt. The metadata's HEIGHT description and the thresholded ABL/ACC sums should not be assumed to define the same operator without checking the actual analysis. No operator was selected by which score looked closest.

## What would close this specific gap

The missing authoritative fields are: exact processed observation vector or cumulative total used for Table3; native daily/hourly interval convention; density and snow/ice conversion used in that processing; missing/subthreshold filtering; and the relation between that target and calibrated debris properties. A validated model cannot be created by substituting the algebraically implied target or an apparently matching density.

For a cleaner melt target, [the paper's §4.3.1](https://tc.copernicus.org/articles/20/1895/2026/) identifies Miage and Pirámide as manual ablation-stake comparisons rather than daily surface-height records. Their [Miage archive](https://doi.org/10.5281/zenodo.3050557) and the paper's Pirámide archive link are concrete alternatives for the parallel target audit. Their forcing/target locations, survey times, density convention and calibration independence still require actual verification; this report does not certify them.

## Reproduce

`resolve_target.py` invokes the pinned earlier reproduction with explicit paths and writes its result plus the source-authority audit. Expected exit2 means unresolved validation, not a program crash. `observation_operators.py` verifies the original archive MD5 and computes the separate measurement operators. Four unit tests pass, including the adversarial density-match example. The source code/datasets were read only. Claude's files, simulations and historical outcomes were untouched.
