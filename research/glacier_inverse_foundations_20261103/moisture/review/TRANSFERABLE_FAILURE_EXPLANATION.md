# Why nonnegative exponential moisture modulation fails on the observed transfer

Post-test analysis only. No parameter was refitted or selected and no additional target was read. This applies to the issued Tower 2 temperature/moisture operators and the released camera comparisons, not certified true moisture, a universal glacier law, or a historically new theorem.

## Exact decomposition

For each period I, normalize the positive-air-temperature dose weights to a probability distribution of observed moisture, P_I. Let P_0 be TRAIN. With the same local TRAIN ablation gain,

    M_beta(I) / M_air(I) = E_(P_I)[exp(beta theta)] / E_(P_0)[exp(beta theta)].

The training reference cancels. The independent reconstruction verifies this identity for the issued beta5.4 to below1e-12 relative error.

If P_I first-order stochastically dominates P_0, every increasing response has at least as large an expectation. Therefore M_beta>=M_air for every beta>=0. Since the nominal air-PDD prediction already exceeds observed ablation in all four periods, no such correction can reduce that period's absolute error.

**The first-order condition holds only in periods2–4.** In period1 the empirical CDFs cross: maximum F_I−F_0 is0.0185448182 at theta0.0664. A higher mean is not sufficient to claim first-order dominance.

## A stronger applicable conditional gate

All four empirical distributions nevertheless dominate TRAIN in increasing convex order. For finite supports this can be checked by

    E_I[(theta−t)+] >= E_0[(theta−t)+]

at every threshold in the union of their support points. The tail-difference function is piecewise linear between knots, constant below the lowest knot and zero above the highest, so knot checks suffice. The minimum is zero at the largest support value in all periods; the minimum below that point is strictly positive in the numerical calculation.

For completeness, on the bounded common support [a,b], any twice-differentiable increasing convex function satisfies

    phi(x)=phi(a)+phi'(a)(x−a)+integral_a^b (x−t)+ phi''(t)dt.

The stop-loss inequality also implies the mean inequality by taking t=a. Every term in the expectation difference is consequently nonnegative. Setting phi(x)=exp(beta*x), beta>=0, proves that the calibrated exponential factor is at least1 **for every nonnegative beta**, not just the tested finite grid. Beta0 reproduces air PDD. This is an application of established stochastic-order mathematics, not a newly invented theorem.

At the locked beta5.4 the correction multiplies air-PDD predictions by:

| Period | Exponential factor | Air-PDD signed error, mm w.e. |
|---|---:|---:|
| Dec16–Dec30 | 1.026019793 | +1.224845539 |
| Dec30–Jan13 | 1.124287692 | +4.009285362 |
| Jan13–Jan27 | 1.095725285 | +2.136641977 |
| Jan27–Feb10 | 1.100343719 | +21.690448657 |

Thus, conditional on these observed covariates and nominal camera totals, the entire beta>=0 exponential family cannot improve over air PDD on any of the four periods. This is stronger than saying that beta5.4 happened to fail. It does not rule out physically meaningful moisture effects, negative beta, other response families, other periods or different measured depths; those are separate hypotheses and cannot be validated by retuning on this test.

## Covariate uncertainty prevents a claim about unknown true moisture

For arbitrary pointwise moisture errors |epsilon(t)|<=delta and unchanged air weights, the worst case for increasing-convex dominance puts every test moisture value down by delta and every TRAIN value up by delta. Those extremal shifts suffice because all tested responses are increasing. A common constant offset applied to both distributions preserves their ordering and cancels from the prediction ratio.

Bisection on the exact finite-support stop-loss criterion gives these conditional symmetric radii, in volumetric moisture fraction:

| Period | Largest symmetric pointwise error radius preserving the ICX gate |
|---|---:|
| Dec16–Dec30 | 0.002276892 |
| Dec30–Jan13 | 0.010776426 |
| Jan13–Jan27 | 0.008395425 |
| Jan27–Feb10 | 0.008472473 |

The all-period radius is therefore approximately0.00228 (0.228 percentage points). Every period's gate fails under the comparison delta0.031 (3.1 percentage points), the manufacturer's typical-accuracy scale recorded in the source audit. Typical accuracy is not itself a certified hard error bound. The calculation demonstrates that source uncertainty is too large to promote the observed-distribution result into a guarantee about actual field moisture. Camera error and TRAIN gain uncertainty are separately reported in the main post-score audit; the positive nominal air error is not an uncertainty-free fact.

The numerical stop-loss tests use declared tolerance1e-12, with1e-14 for the radius search, not outward-rounded interval arithmetic. The smallest nontrivial nominal tail margins are about1.05e-6,7.48e-6,1.80e-6 and3.85e-6, well above those tolerances. The reported radius digits are computational precision, not sensor resolution.

## Verification and scope

Six synthetic tests pass: a translated distribution with exact exponential factor for every tested beta; a crossing-CDF example that passes increasing-convex order; a higher-mean example that fails it; equal distributions/zero weights; and the exact0.05 symmetric uncertainty radius for a0.1 translated fixture. The descriptive beta0–40 grid is retained in the JSON but selects no new model and is not the proof of the all-beta result.

Artifacts: `transferable_failure_analysis.py`, `test_transferable_failure_analysis.py`, and `TRANSFERABLE_FAILURE_ANALYSIS_ICX_UNCERTAINTY.json`. Earlier FOSD-only and ICX-only JSONs are retained as versioned analysis records, with the uncertainty JSON the complete current explanation.
