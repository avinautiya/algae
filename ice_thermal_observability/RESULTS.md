# Results and scientific decision

31 October 2026. All physical numbers below are synthetic computations for the declared fixed-flow heat equation. No new field observation or OpenFOAM simulation is reported by this package.

## Proved result

The continuum theorem establishes the exact finite-bulk-observation placement criterion for prescribed unidirectional velocity, including no-slip shear. Its proof uses a weighted diffusion spectrum, nonzero coupling of every mode to bulk temperature, linear independence of exponentials and a positivity-preserving construction. The explicit analytical plug-flow witness matches both calibration readings to 3.25e-19 K while changing regional heat by 0.0043832 W. This is a positive smooth counterexample, independent of the finite-volume and LP calculations.

Historical priority is distinct from proof. Generic inverse ambiguity, Graetz spectra, optimal recovery and LP duality already exist. The exact ice-directed regional observation criterion is the candidate contribution requiring expert/full-text comparison; no first-ever claim is certified.

## Executable regional forecast

The primary case has H=2 mm, U=0.025 m/s, W=0.1 m, Re_h=50, physical constant water/ice properties, initial excess bounded in [0,0.2] K, two exact bulk observations at 0.02 and 0.04 m, and a target at 0.06–0.12 m. A steady scenario is projected for 60 s.

| Quantity | Result |
|---|---:|
| Sharp conditional regional heat | 0.232241–0.262765 W |
| Ratio of attainable extremes | 1.13144 |
| Conditional mean ice lowering | 7.58267–8.57931 micrometres |
| Scalar minimax lowering | 8.08099 micrometres |
| Conditional worst-case absolute error | 0.498319 micrometres |
| Maximum primary evaluated primal–dual gap | 2.22e-16 W |

Thus identical calibration permits a 13.14% difference between the regional heat extremes in this class. This is not a measured glacier uncertainty or a 13% improvement over an existing model. Both extreme profiles are supplied and conservatively replayed. Their upstream profiles are admitted boundary fields, not demonstrated products of real glacier meteorology.

Adding independently bounded nonnegative surface heat of 0–100 W/m² on 0.04–0.12 m increases the upper regional heat to 0.271751 W. The energy balance includes the supplied source power. Only part of that source reaches the target bed; the rest remains advected downstream under this model.

## Every predeclared comparison

All methods use the same observations on 120 seeded smooth positive inlet profiles. No training or favorable case selection is involved.

| Predictor | Mean absolute heat error | Maximum absolute heat error |
|---|---:|---:|
| Uniform-entrance amplitude fit | 9.824 mW | 25.124 mW |
| Classical two-mode Graetz fit | **1.317 mW** | **4.446 mW** |
| Conditional interval midpoint | 1.592 mW | 5.852 mW |

The classical two-mode baseline is better on average and in maximum error over these 120 smooth draws. The midpoint's proved property is optimal worst-case loss over the broader declared information class, not best average prediction under this particular synthetic distribution. All admitted true outputs lie within their intervals by construction/verification; this is not statistical confidence coverage. Three two-mode reconstructions violate the admitted inlet-temperature bounds even though their target predictions remain inside the conditional ranges.

The benchmark therefore does **not** establish a new superior point predictor or demonstrate that ordinary modal modeling is inadequate. Its positive result is an executable information-limit calculation and a precise target-specific measurement requirement.

## Verification and refinement

* **24 tests pass**, including actual CLI execution, output preservation, positive smooth continuum witnesses and future-data rejection.
* **80 additional continuous LP cases** agree with independent exhaustive active-face enumeration; maximum discrepancy 3.33e-16 in the declared scalar units. The tests also perform 32 independent small cases.
* An independent sine-series solution gives second-order convergence of plug-flow regional heat: relative forward error 9.39e-4 at 24 cells to 1.46e-5 at 192 cells.
* A separate ODE shooting calculation checks the first two no-slip thermal eigenvalues and second-order grid behavior.
* Direct quadrature of bed heat agrees with the conservative regional functional. A separate stiff ODE integration checks surface-source transport and its energy ledger.
* The supplemental **fixed-data** grid audit holds the same observation vector on every mesh. The shear reference uses 768 cells; the 384-to-768 observation changes are 3.27e-8 and 3.12e-8 K. These are error indicators, not certified PDE error bounds. At 96–192 cells the fixed-data lower/upper power bounds move by approximately 3.25e-6/3.61e-6 W, much smaller than the primary 0.0305 W conditional interval.
* Four thermal-age ranges are flagged below numerical resolution. One tiny reversed range was corrected transparently without changing primary forecasts, test cases or baseline scores; see FAILURE_LOG.md. Small numerically unresolved ranges do not prove exact continuum identifiability.

## Figures and reproducibility

Seven scientific figures have PNG and vector PDF exports, underlying JSON data, captions and hashes. They cover the continuum counterexample, attaining shear-flow forecasts, information/noise, every comparator, independent verification, thermal development and source heating. The PNGs were visually inspected and crowded ticks/panel placement corrected. Plot layout was refined after the study; the figure manifest records the final plotting script separately from the study's source snapshot.

The protocol was committed before the first production run at 70e73fb053b0dd29a1cc98573cfa3a0a483f3aa6. The numerical-width correction was committed before the identical rerun at 16d6ad05e15b8934d4340b9e90b2bf840af47e88. The unchanged protocol SHA256 is 254fc727de897721c255c0132c280d25b3000f863ab9a980ae2668454fa80cac. Original study_v1 records are retained alongside corrected study_v1_1. CSV storage uses canonical LF in Git; this does not alter its numeric values.

The corrected study's partial resource output reports 1.90 s wall time. Peak memory was not obtained because the macOS resource wrapper was denied a system query; the model completed successfully. No speed comparison against CFD is inferred from this timing.

## Scientific use and remaining gate

This provides an actual theorem, numerical implementation, conservative forecast and a reproducible comparison. Its useful role is to distinguish aggregate thermal calibration from support for regional ice erosion, then quantify what a declared model/data class can predict.

The next physical gate is whether independently generated thermal histories in the verified CFD produce forecast-relevant ambiguity beyond established developing-flow models. The current straight-film result cannot establish that turbulent glacier flow does so. Connecting the theorem-based recovery to a verified CFD-derived response operator, justified profile/source constraints and withheld physical configurations is necessary before an environmental improvement claim. Publication significance and Nature-level suitability are not established by these synthetic results alone.
