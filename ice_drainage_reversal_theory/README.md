# Delayed drainage reversal in an eroding ice-lake outlet

Theory packet, 31 October 2026. This addition changes no solver, checkpoint, job, frozen test, or previous result. No new simulation, training, or parameter campaign was run for it.

**Candidate physical result:** adding a higher, non-eroding outlet can lower an ice-dammed lake initially yet leave it with more water after drainage, because less water and heat pass through the persistent eroding outlet. The added outlet can run dry naturally; no prescribed valve reversal or external pumping is required.

This packet supplies proofs for stated reduced models. It does not establish historical priority, occurrence at a real glacier, superiority over existing coupled drainage models, or suitability for Nature. Those are separate questions. In particular, changing the prediction target from regional heat to residual lake water does not repair the old heat-prediction benchmark.

## What is actually derived

1. An exact retention identity for arbitrary monotone outlet ratings and a class of nonlinear incision laws. It proves a positive final lake-level effect of a finite bypass and bounds the effect from the incision-to-drawdown ratio.
2. A corresponding result for arbitrary lake hypsometry when incision is proportional to discharge, including common finite-duration water inflow. It predicts final retained volume without time-stepping the hydraulic rating law.
3. Exact finite-time counterfactual forecasts for a linear rating, including the time at which the initially beneficial drawdown changes sign. This includes an open range of stable parameters, not just a critical balance.
4. A right-half-plane zero in the local input-output response when incision increases with discharge. This identifies a delayed response that a fixed-bed reservoir cannot reproduce. Transfer-function algebra is established control mathematics; this application is a candidate contribution.
5. A separate, precise framework for improving the original regional-heat point predictor. It explains which guarantees are possible and why no average-performance win has been established from the existing data.

Start with [THEOREMS.md](THEOREMS.md), then [PHYSICAL_MODEL.md](PHYSICAL_MODEL.md), [PREDICTOR.md](PREDICTOR.md), and [PRIOR_ART.md](PRIOR_ART.md). [DECISION.md](DECISION.md) states the publication and environmental gate.

## A conditional number, not a field result

In the constant-area model, let the incision-to-drawdown ratio be beta = 1/2, initial lake level H0 = 3 w0, and initial persistent sill b0 = 2 w0. Here w0 is the initial head above that sill. Without a bypass the final lake level is w0. A finite bypass removing C = A w0/2, while the main outlet remains wet, changes the final level to 3 w0/2.

Thus the model retains **50% more water than its own no-bypass residual**, despite removing additional water early. This is exact arithmetic under the model, not a new computed benchmark or an environmental percentage. The same experiment must not be reported as a 50% improvement in prediction accuracy or a 50% change in glacier melt.

## Relationship to prior work

Lake-level/incision competition, drainage instability, channel capture, simultaneous overspill and incising drainage, and fluid-network paradoxes already exist. This packet's possible distinction is the exact finite-bypass retention and crossing-time statements. The literature screen did not certify that distinction as historically original. An existing coupled model may reproduce the phenomenon without modification.

The current scientific decision is **retain as a specific theory candidate; do not label it Nature-level established novelty and do not launch a production campaign on that premise**. The earlier thermal results remain unchanged and negative for average superiority over their strongest modal baseline.
