# Declared applicability and failure classification (written before any parameter-campaign or challenge case was analysed; committed before the challenge cases run)

Applicability flags, evaluated per case on the reduced model's own diagnostics over active-edge snapshots (route discharge > 1 % of the mean inflow):
* `quasi_steady`: 90th percentile of `r_qs` (= gravity-wave crossing time / stage-change time) <= 1.0.
* `subcritical`: fraction of active snapshots with at least one station at Froude >= 1 is <= 0.05.
* `storage_neglected`: 90th percentile of channel volume / feeding-node volume <= 0.05.
A case is **inside** the declared envelope when all three hold, otherwise **outside**. (No flag is a validity proof; they mark where the stated assumptions are strained.)

Outcome classes for reduced runs against the 200-cell reference (400 cells where available):
* `ok`: completed, route agreement (see below) and hydrograph RMSE / mean inflow <= 0.10.
* `failure_numerical`: the run failed, stalled, timed out or hit its evaluation cap while inside the envelope.
* `failure_outside_regime`: the run failed, or exceeded the 0.10 hydrograph tolerance, or chose the wrong route, while outside the envelope.
* `failure_unexplained`: inside the envelope, completed, but exceeded the tolerance or chose the wrong route.
Route agreement is judged only when the reference's route margin |V1 - V2| / (V1 + V2) exceeds 0.10; closer cases are labelled `near_tie` and are reported separately, not as failures. A reference run that fails, stalls or times out makes the point `reference_failed` (reported, excluded from error statistics, never silently dropped).
Reduced configurations compared: r9/dt900 s and r17/dt450 s. No case is removed after execution.
