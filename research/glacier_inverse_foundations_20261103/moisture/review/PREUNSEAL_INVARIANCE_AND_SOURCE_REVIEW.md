# Independent pre-unseal review

Reviewed `experiment.py`, `anomaly_index.py`, the final Tower2 freeze and source
contract without calling a target reader. No must-fix execution/data-leakage
defect was identified in this review. This is approval of the declared test's
implementation, not a prediction of its outcome.

Final reviewed freeze SHA256:
`d543df0265b6e626ce9f2f29b35cf4bf8961648f85a3061008c1fd83bb0af2fd`.
All frozen source, canonical input, external thermistor, camera-member and target
loader hashes match. Camera bytes were hashed only; no outcome was decoded.

## Equal local calibration and source gates

Nine model names and their source parameters are frozen. Five global scalar
hyperparameters come only from Tower1's already released development outcomes;
the four untuned comparators retain parameter0. The primary beta is5.4, not a
parameter refitted to Tower2. All nine methods receive precisely the same one
Tower2 TRAIN-ablation scalar and one positive local amplitude gain. Four later
periods use that gain unchanged. The source-calibrated calendar comparator is
retained even though it has lower Tower1 development MAE than the primary.

The sequence binds normalized artifact SHA, source/spec SHA and access role
before TRAIN release; it then seals all nine calibrated prediction arrays before
EVALUATION access. Calibrated predictions contain no later target values. Exact
endpoint roles and the original metres-ice times917 operator remain bound.

Exposure integrates the positive part of piecewise-linear native air segments
exactly, with explicit18h maximum gaps and no extrapolation. Moisture is held
left over complete half-hour intervals. All raw moisture rows are retained.
The ratio baseline's+.031 is a declared feature regularizer, not an observation
repair. Positive gradient is explicitly an empirical feature and carries no
signed-heat or energy-conservation claim. Source flags for nonpositive/uncertain
moisture or soil-temperature conditions are limitations, not deleted periods.

This remains a same-campaign calibrated stake-transfer test using observed
contemporaneous covariates. It is not unseen-weather or independent-glacier
validation. Source physics and sensor-depth/temperature/moisture limits are
retained in the protocol. No general manufacturer-accuracy certification follows
from mathematical offset invariance.

## Exact offset and reference results

Let `a(t)>=0` denote nonnegative temperature exposure, theta(t) the raw moisture
proxy, andr a fixed reference. For positivew define

```
D_w(I;r) = integral_I a(t)*w(theta(t)-r)dt,
M_test = Y_train * D_w(test;r)/D_w(train;r).
```

For exponentialw(x)=K exp(beta*x), the factorK exp(−beta*r) cancels from the
ratio. Consequently the **calibrated prediction** is independent ofr. Constant
additive sensor offsets also cancel, even if each site has a different offset
constant over its own calibration/evaluation record. These are exact real-number
identities; floating-point evaluation still requires representable exponents.

For a data-equivariant center such as own-TRAIN mean, replacingtheta bytheta+b
also replacesr byr+b. The anomalytheta−r is unchanged. Therefore **every fixed
function of a centered anomaly**, including centered softplus, shares constant
additive-offset invariance. Exponential uniqueness must never be claimed from
this centered-offset property alone.

Centered softplus additionally fails general reference independence: for
`w(x)=log(1+exp(beta*x))`, changingr at fixedtheta does not generally multiply
all weights by one common factor. For beta1 and point exposures theta_train0,
theta_test1, r0 gives`log(1+e)/log2`, whereasr1 gives
`log2/log(1+e^−1)`. They differ. Beta0 or constant/non-discriminating data are
degenerate exceptions. The independent synthetic review confirms the distinction
for the actual frozen nonzero beta/gamma values.

## Characterization and its necessary scope

Suppose one **fixed, positive, continuous function**w on the real line is used
asw(theta−r), and require the calibrated ratio to be independent of arbitraryr
for **every** nonzero nonnegative TRAIN/test exposure distribution. Taking point
exposures atx and0 gives, for any translationb,

```
w(x+b)/w(b) = w(x)/w(0).
```

Thusg(x)=w(x)/w(0) satisfies`g(x+b)=g(x)g(b)`. Positivity permits taking logs;
continuity makes the additive Cauchy solution linear, so
`w(x)=K exp(beta*x)`. Conversely that family has the required property. Connected
interval versions require the corresponding local-translation/domain assumptions.

This is an established functional-equation characterization, not a new historical
theorem, a glacier law or a uniqueness theorem for all offset-invariant indices.
It assumes arbitrary additive references applied through the same fixedw;
it does not apply to all possible feature normalizations. For example, a
shifted-ratio power has a reference denominator that is merely a common scalar,
so its calibrated ratio can also be reference-independent. It is not a fixed
function oftheta−r and generally lacks raw additive-offset invariance. Those
facts are compatible with the stated characterization.

## Uncertainty scope

A multiplicative sensor change theta→b*theta+c replaces the effective beta by
b*beta; fixed transferred beta is not immune to that calibration gain. Drift,
temperature response, nonlinear calibration, censoring and probe-depth changes
also remain. For fixedbeta and a separately assumed pointwise reading error
`|delta(t)|<=epsilon`, uncentered dose perturbations are between
`exp(−|beta|epsilon)` and`exp(|beta|epsilon)`; hence the calibrated ratio lies
within factors`exp(±2|beta|epsilon)`. Reference cancellation prevents adding an
extra independent mean-reference factor. This is a conditional bound, not a
field confidence interval or a manufacturer guarantee in rock debris.

No energy budget or causal wetting effect follows from the index. Moisture
storage may respond to meltwater, and existing moisture-inclusive physical
models already include competing conduction/capacity/latent effects. The frozen
positive fitted coefficient is an empirical result from Tower1 development;
its fresh Tower2 performance is an unanswered test at the time of this review.
