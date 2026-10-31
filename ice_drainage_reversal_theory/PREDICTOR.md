# What can establish a superior predictor

There are two different targets here. They must not be pooled into one success claim.

## A. Original target: regional heat from upstream bulk observations

The immutable 120-profile study at commit 17505c8f17f8073a0beb2ab17944211e760f5dfc reports mean absolute heat errors of 1.317 mW for the classical two-mode Graetz fit and 1.592 mW for the conditional interval midpoint. The classical baseline won. This packet does not change those scores or rerun that study.

### A theorem cannot manufacture the missing distribution

Let X be the available observations, Y the true regional heat, and b(X) the existing predictor. The optimal population predictor for squared loss is E[Y|X]; for absolute loss it is a conditional median. The feasible interval alone does not specify either conditional distribution. A minimax midpoint is therefore not a general solution to average prediction.

If b(X) is itself an admitted target value, no nonzero change at that X can improve every admitted population: choose a population supported on an admitted state with Y=b(X). The original predictor has zero error there and the changed one does not. In a convex feasible class, every point between the attainable target extremes is admitted. An average-performance guarantee must specify a population or obtain additional information. It cannot follow from model complexity or the dimension of a neural network. When b is outside a valid interval, the projection result below is an important exception.

### Exact residual correction identity

For squared loss, put R=Y-b(X). With additional information Z available at the forecast origin, the population-optimal correction is

\[
g^*(X,Z)=E[R\mid X,Z],\qquad \hat Y=b(X)+g^*(X,Z).
\]

Conditional expectation gives

\[
E[(Y-\hat Y)^2]=E[R^2]-E[g^{*2}].                        \tag{A1}
\]

Thus the oracle strictly improves mean squared error exactly when predictable residual structure exists. This is established orthogonal-projection mathematics, not a new theorem claimed for this project. It does not say a fitted correction estimates g* accurately.

For an arbitrary fitted correction g, conditional on a training dataset independent of evaluation cases,

\[
Risk(b+g)-Risk(b)=E[g^2]-2E[Rg].                         \tag{A2}
\]

The correction needs enough independently measured residual correlation to exceed its own variance. This is the actual hurdle. A complex PINN can make the right-hand side positive and worsen the predictor. For the original MAE metric, a conditional-median correction must be evaluated; a win in MSE cannot be substituted for a win in MAE.

A meaningful physical Z could contain an independently constrained upstream thermal state, forcing history or resolved mixing information that the current two-mode inference lacks. Future downstream temperature, future heat transfer and test labels cannot be used as forecast inputs. Mode count or a learned residual alone is not a new physical discovery.

### A distribution-free safe correction, with a decisive limitation

Suppose [L(X),U(X)] is a valid feasible interval for Y, including all required observation and model errors. Set p=clip(b,L,U). For every admitted Y,

\[
(b-Y)^2-(p-Y)^2\ge(b-p)^2\ge0,\qquad |p-Y|\le|b-Y|.     \tag{A3}
\]

Proof: if b>U then p=U and (b-Y)^2-(U-Y)^2=(b-U)(b+U-2Y)>=(b-U)^2; the other side is symmetric. If b is inside, nothing changes.

The existing RESULTS.md explicitly states that all two-mode target predictions in its 120-case comparison are already inside the conditional ranges, including the three cases with inadmissible reconstructed *profiles*. Clipping therefore does not improve that reported comparison. Profile infeasibility must not be passed off as target-prediction infeasibility. Also, the computed intervals are sharp for a declared discrete thermal model; they are not currently certified environmental error bounds.

**Conclusion for A:** a rigorous framework exists, but no new superior average point predictor has been established. No paper should say otherwise. The previous two-mode predictor remains the strongest demonstrated average predictor in that experiment.

## B. New target: counterfactual lake level and final retained volume

### An analytical superiority result against a precisely defined baseline class

Consider any fixed-bed reservoir with the same prescribed common inflow and added nonnegative bypass,

\[
A(H)\dot H=I(t)-F(H-b_{fixed})-u(t),\quad A(H)>0.
\]

For the higher fixed auxiliary sill, write u=G((H-h_s)_+). Standard scalar comparison shows that adding the outlet cannot raise H above its no-bypass trajectory: at a first coincident level its derivative is lower by u/A. Thus these models predict Delta H<=0 at all times. A prescribed nonnegative bypass has the same ordering property.

Propositions 1 and 2 prove Delta H_infty>0 for the admitted eroding-outlet class. A predictor enforcing fixed-bed ordering must therefore have an error in the sign of this counterfactual response. For proportional incision and constant A, the exact predictor is

\[
\widehat{\Delta V_\infty}=\frac{\beta}{1-\beta}C.         \tag{B1}
\]

For any probability distribution of admitted nonempty cases with beta>0 and C>0, finite second moment, and exact supplied parameters, (B1) has zero model-class error. Every comparator predicting Delta V_infty<=0 has squared error at least [beta C/(1-beta)]^2 in each case, so its average squared error is strictly higher. Absolute loss is also strictly higher.

This is a conditional analytical dominance result on a new observable. It is **not** an empirical win on the original regional-heat dataset, not a forecast with unknown parameters, and not superiority over coupled eroding-outlet models. Those models can already reproduce (B1) when supplied with these equations. Including them is mandatory in any future claim about improvement over existing models.

There is also a robust point-forecast version. If independent admissible bounds give beta>=beta_min>0 and C>=C_min>0, set L=C_min beta_min/(1-beta_min) and forecast p=L for the positive retained-volume change. Proposition 1 gives Y=Delta V_infty>=L. For every comparator q<=0,

\[
|Y-p|=Y-L<|Y-q|,\quad
(Y-q)^2-(Y-p)^2\ge L^2>0.                              \tag{B1a}
\]

Thus a conservative positive forecast strictly dominates the entire nonpositive-response class in both absolute and squared loss for every admitted case, even without the exact beta. This is still a conditional model-class statement. No independent site bounds establishing beta_min and C_min have been supplied here; if the admitted bounds include zero, the strict advantage is not guaranteed. It also gives no dominance over another predictor that uses erosion feedback.

For general nonlinear E, Proposition 1 gives an exact path-weighted prediction and certified class bounds if beta bounds are independently supplied. For known hypsometry and proportional incision, solve the scalar monotone root in Proposition 2 rather than time-step the hydraulic rating. This is an efficient exact calculation inside the reduced model, not a replacement for its uncertain physical closure.

### Finite-horizon forecasts

Proposition 3 gives an exact Delta H(t) and crossing time from A,k,kappa and a *forecast* bypass schedule. For a naturally activated auxiliary outlet, its rating must also be forecast. After-the-fact pulse volume and centroid are sufficient for post-pulse state reconstruction in that linear class, but are not available in advance merely because the formula uses them.

A useful eventual output would include H(t), residual volume, sign-reversal time, predicted sill level, and the parameter/forcing conditions under which the prediction is valid. It must report a range when those parameters are uncertain. This packet supplies formulas, not a newly executed or site-validated forecasting application.

### Link to the existing thermal work through a forecast-sensitive loss

Let C_m=rho_i L_f S_c [J/m]. An additional heat-to-ice error r(t) [W] enters b_dot as -r/C_m. In the constant-coefficient response with a=k/A and lambda=k(E_Q-1/A), first-order sensitivity of the lake level at T is

\[
\delta H(T)=\int_0^T K_T(s)r(s)\,ds,
\qquad K_T(s)=-\frac{a}{C_m}
 \frac{e^{\lambda(T-s)}-1}{\lambda}.                    \tag{B2}
\]

At lambda=0, K_T(s)=-a(T-s)/C_m. This follows directly by solving delta w_dot=lambda delta w+r/C_m and delta H_dot=-a delta w.

Consequently small average heat errors need not imply small lake-level errors: early errors and different temporal covariance can be weighted very differently. The squared forecast error is E[(integral K_T r)^2], which depends on the full residual covariance, not just E[r(t)^2]. This is standard goal-oriented/adjoint reasoning applied to the stated drainage dynamics. It is not a new generic learning theory.

This provides a defensible connection between the thermal solver and a consequential forecast. It does not justify changing the old heat-comparison metric until a preferred method wins. Both old and new observables must be reported separately.

## No production authorization inferred from these formulas

The publication decision is in DECISION.md. No simulations, fitting, threshold tuning, or withheld evaluation have been requested or launched by this packet. Existing experiments are neither relabeled nor erased.
