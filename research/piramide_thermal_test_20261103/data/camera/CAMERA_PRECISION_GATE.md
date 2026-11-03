# Camera precision and correction evidence gate

**No independently justified numerical camera-endpoint error bound or standard deviation was found in the inspected material.** Do not convert the10cm ruler spacing,24MP camera specification or Tracker display precision into a measurement uncertainty.

## Primary evidence inspected

The final [Melo-Velasco et al.2025 paper](https://doi.org/10.1029/2025JF008360), section2.2.3, establishes the exposed-stake measurement operator and917kg/m³ ice conversion. Automated top/bottom selections were manually corrected when lighting/overnight identification failed. The paper mentions end-season manual measurements but provides no camera-error distribution there. The published±2cm uncertainty in table1 concerns debris thickness/sensor depth, **not camera ablation**.

The official article lists a45.6MB [Supporting InformationPDF](https://agupubs.onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1029%2F2025JF008360&file=2025JF008360-sup-0001-Supporting+Information+SI-S01.pdf). The browser could not retrieve it (cache miss/inaccessible). It was not downloaded or read, so absence of a precision estimate inside it is **not established**.

The verified [Zenodo archive](https://doi.org/10.5281/zenodo.11581073) index contains186nonempty members, including three processed camera MAT files. No image file or separately named manual-stake/Tracker/correction/uncertainty record appears in that index. This does not prove that no manual observation is embedded in a MATfile; camera outcome cells remain sealed. The source README defines distance in metres but supplies no quantitative camera precision.

## Reset/correction status

Native timestamp uniqueness and ordering are verified. That is not evidence of physical stake-reference continuity or correction accuracy. A later stage may inspect only authorized role endpoints; a complete interior reset audit needs explicit additional access and source correction provenance. Negative target increments must be retained/flagged, not removed to improve agreement. If a reference reset is discovered, invalidate the affected comparison unless an independently supported correction was predeclared; do not choose replacement endpoints after scoring.

The March manual observation is not presently an independently bound target file with timing, precision and reference continuity. It therefore cannot be called an independent validation of processed camera ablation. No authors were contacted.

## Honest robustness statement

A nominal error comparison may be reported with camera uncertainty explicitly unknown. A sensitivity threshold delta-star is useful, but it is **not an empirically established precision or confidence interval**.

For one-gain models whose issued transfer ratiosf,g are fixed before evaluation, use separate uncertainty estimands:

- For uncertainty in already-issued predictions, hold their numerical values fixed and perturb only evaluation observations.
- For a counterfactual rerun of the full training/calibration procedure, perturb all three shared stake endpoints. With positive training ablationX and evaluation ablationY, ranking fX versusgX depends on Y−mX, m=(f+g)/2. Per-reading errors separately bounded bydelta[mICE] move this quantity by at most2*917*(m+1)*delta. Require positivity X−2*917*delta>0 as well.

The independent review agent implements these exact shared-endpoint checks. Agreement for delta below a computed breakpoint supports only a **conditional sensitivity result** until an external error bound is justified. Repeated camera samples or engineering scenarios do not create independent glacier-level replication.
