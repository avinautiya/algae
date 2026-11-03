# Continuous calibration-transfer certificate

Let g be an unknown parameter vector in the bounded nonempty polytope P={g:Ag<=b,l<=g<=u}. Let c·g be integrated TRAIN heat and f·g integrated later heat on one fixed linear/phase branch. Require inf_P c·g>0 and inf_P f·g>=0. With a positive measured TRAIN amount Y, the calibrated later amount is Y(f·g)/(c·g).

Set t=1/(c·g), y=t*g. Every feasible g has a corresponding t>0 satisfying c·y=1, A*y<=b*t, l*t<=y<=u*t. Conversely every transformed feasible point with t>0 gives g=y/t in P. Thus minimizing/maximizing f·y is a linear program computing the exact extrema over the DECLARED continuous polytope. This is classical Charnes–Cooper linear-fractional programming, not a historically new theorem.

The code first checks the original problem's feasibility and uniformly positive denominator, then solves both transformed programs and checks recovered witnesses against original constraints. Finite bounds exclude a nonzero y solution at t=0. Floating-point tolerances are reported; there is no outward-rounded interval certificate. Original constraints themselves may be only necessary outer relaxations of noisy-temperature/phase physics. The bounds then enclose that relaxation, not a guarantee of real-glacier coverage.

A scalar point prediction's worst error over this interval is at least half its width, attained by the midpoint. This is a conditional minimax statement over the declared family, not a statistical confidence interval or proof of superior average prediction. Empty feasibility falsifies the assumptions; it never yields a zero-width confident forecast.

Unlike a finite scenario grid, the method cannot miss an admissible interior or edge extremum of this specified linear-fractional problem. It does not solve arbitrary nonlinear temperature-dependent conductance, unknown geometry, switching phase history, unobserved forcing or temporal parameter drift. Those require explicit branches, new constraints or abstention. A calibrated ratio may be unique even when parameters are not (e.g.f=3c). Conversely a TRAIN scalar may leave the future ratio broad.

Next application must fix source sensor-error inequalities, geometry, linear branch, model-discrepancy allowance and all comparator information before outcome release. Using future observed temperatures makes an observation-conditioned retrieval, not an advance-weather forecast.

## Coefficient uncertainty extension

For nonnegative g and supplied coefficient boxes c_low<=c<=c_high and f_low<=f<=f_high, with a uniformly positive c_low·g, every admissible ratio obeys

    (f_low·g)/(c_high·g) <= (f·g)/(c·g) <= (f_high·g)/(c_low·g).

Minimizing the left and maximizing the right over the same feasible P gives a conservative outer enclosure, using two more linear-fractional programs. It does not assume physical sensor errors are independent: it discards their correlation, enlarging the set. It does require valid bounds for every true coefficient. Uncertainty in A,b itself must also already be conservatively included in P; the routine does not magically convert a nominal regression into a valid feasible set. Negative future heat still requires explicit phase handling and is refused by this branch. The independent corner-enumeration test covers the implemented box formulation.
