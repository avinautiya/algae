# Arolla H1D conditional measurement-and-clock robustness audit

**No model robustly dominates another.** Every pair's absolute-error ranking can reverse under admissible bounded measurement errors, in every one of the81 frozen clock scenarios, under both uncertainty interpretations. No observations, coefficients, clocks or targets were tuned. This is retrospective sensitivity analysis, not a new predictor or a statistical interval.

## Scope and precise estimands

The audited input is the unchanged `model/results/frozen_ablation_comparison/RESULTS.json`; its SHA256 is recorded in `RESULTS_V2.json`. Only the three native-ice baseline methods and primary H1D stake are included. All81 originally declared timing combinations must exist exactly once. Recorded calibration increment X=846mmw.e.; evaluation increment Y=315mmw.e. Source0.1m × native900kg/m³ yields90mmw.e., but the source does **not** identify whether its uncertainty is a per-reading bound, per-interval bound, standard deviation or another summary. Treating it as an absolute bound is conditional, not author-certified.

Two different questions must not be conflated:

1. **Input-error/recalibration sensitivity:** for each admissible true calibration increment X*, compare the implied prediction fX* against true evaluation increment Y*. This audits the training calibration's sensitivity. It is the primary requested minimax audit and matches the earlier shared-reading ranking audit.
2. **Already-issued fixed forecast:** keep the actual forecast fX fixed and vary only the true evaluation increment Y*. This measures the error of that issued forecast. It gives a different worst-case ranking. The audit reports this separately to prevent a misleading robustness claim.

No choice made after looking at these held-out losses constitutes a prospective validation.

## Exact bounded calculation

Let true cumulative readings equal recorded readings plus e0,e1,e2, each in[−δ,δ]. Then X*=X+e1−e0 and Y*=Y+e2−e1. With fixed driver ratio f, the input-error residual is

r=fX*−Y*=fX−Y−f e0+(f+1)e1−e2.

The **same** middle error enters both increments with opposite signs. It cannot be independently perturbed twice. Because f≥0, the exact residual radius is(2f+2)δ and exact maximum absolute residual is |fX−Y|+(2f+2)δ. Under a distinct per-interval interpretation X*=X+u,Y*=Y+v with u,v independently bounded byδ, the radius is(f+1)δ. Both are analytic vertex extrema, not sampled approximations.

For pairwise comparison,

ra²−rb²=(fa−fb)X*[(fa+fb)X*−2Y*].

X* remains positive everywhere in the declared boxes. The last bracket is affine; its bounds follow directly from the same shared errors. Consequently its sign decides strict absolute-error dominance throughout the whole box. A sign reversal at corners proves ranking non-identifiability; checking only nominal errors would miss it. The code refuses boxes not establishing physical nonnegative increments rather than clipping them silently.

For an issued forecast, residual=fX−Y* and the exact radius is2δ per reading orδ per interval. Calibration uncertainty does not change an already-issued numerical forecast.

## Results across all81 frozen clocks

All values below are mm water equivalent. These are adversarial conditional error bounds, not expected errors or probabilities.

|Uncertainty interpretation|Method|Input-error sensitivity: maximum absolute error|Fixed issued forecast: maximum absolute error|
|---|---|---:|---:|
|Per reading|Constant rate|378.000|288.000|
|Per reading|PDD|326.997|281.590|
|Per reading|Net shortwave|357.080|270.750|
|Per interval|Constant rate|243.000|198.000|
|Per interval|PDD|214.293|191.590|
|Per interval|Net shortwave|223.915|180.750|

PDD has the smallest maximum input-error-sensitivity loss, while net shortwave has the smallest fixed-issued-forecast maximum loss. This difference is an estimand distinction, **not evidence that PDD fixes the field-model problem**. Within individual clocks, net shortwave has the smallest conditional input-error maximum in42/81(per reading) and55/81(per interval) scenarios; PDD in39/81 and26/81. These counts are not probabilities.

For each interpretation there are3pairs×81clocks=243pairwise tests. **All243 permit ranking reversals.** Thus a conditional minimax ordering supplies no strict field-superiority claim and cannot establish which method has lower actual error. Unknown uncertainty semantics, one primary event and measurement support limit the inference.

## Verification

`python3 -m unittest -v` passes5tests. A separate exhaustive-corner implementation checks1000deterministic random parameter boxes, each at all8 shared-reading corners or4 per-interval corners. It computes the true increments and direct squared errors without calling the analytic bound formula. It also checks the shared-middle sign, fixed-forecast distinction, equal-model identity and invalid/nonphysical-box refusals. This verifies the mathematics implemented here, not measurement accuracy or glacier generalization.

Reproduce into a **new** output path:

`python3 audit.py ../../model/results/frozen_ablation_comparison/RESULTS.json NEW_RESULTS.json`

`RESULTS.json` is the preserved first audit; `RESULTS_V2.json` adds the explicit fixed-forecast comparison and labels both estimands. No frozen source result was overwritten.
