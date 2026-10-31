# Passive contribution for review after the current CFD milestone

Keep Claude's current task running. This package edits no existing solver, case, checkpoint, build script, environment or frozen result. It has no runtime connection to OpenFOAM and launches no CFD process.

The new `ice_incision_forecast` package is standalone for prediction/proof checks. NumPy/Matplotlib are needed only for its six-figure script. It does not depend on the earlier lake, response or forecast-theory packages.

When the current milestone reaches a safe review point:

1. Inspect THEOREMS.md and the conditional model assumptions. Check the conservation and sharpness proofs before integration.
2. Run the tests and continuous-polytope validation. Reproduce the demo forecast and figures in new output directories.
3. Keep synthetic figures/results labeled as such. Do not count them as field or resolved-CFD evidence.
4. Decide whether effective bed heat-transfer conductances can be bounded from the verified model/independent calibration. Use actual water-to-bed power and appropriate bulk temperature; do not mix wall gradient, point temperature, perimeter heat and bed-only heat definitions.
5. Test whether a common forcing multiplier preserves the conductance shape. If not, do not apply the forcing-uniform theorem outside its assumptions. A different time-dependent uncertainty class needs a new bound/attainability argument.
6. Preserve the known measurement limitation: uncertain near-bed point sensors do not supply exact bulk-temperature prefix anchors. Ideal observations in the figures illustrate information structure, not achieved field precision.
7. Assess whether the bounds are physically useful and sufficiently tight for the declared forecast. Stronger physical/observational constraints may narrow them; their validity requires that the true state remain within the assumed class.

The midpoint regional forecast minimizes worst-case scalar error over the admitted class. It is not a probabilistic mean and does not promise smaller error on every actual glacier. Threshold times refer to declared mean-incision levels, not independently established flood or channel-collapse events.

This is a mathematical/predictive contribution to complement the physical project. Do not abandon the active free-surface and melting verification, replace it with this reduced tool, or claim a new glacier instability from these results.
