# Threshold-aware averaging for a two-outlet lake: derivation, status of each statement, validity limits

Setting: lake area A, stage H, outlets i = 1,2 with sills s1 < s2 (gap Δ = s2 − s1), quasi-steady ratings Q_i(H − s_i) (critical control, increasing, Q_i(0)=0), inflow I(t) with heat flux I·T(t), well-mixed lake temperature T_L. `A dH/dt = I − Q1(H−s1) − Q2(H−s2)`. Outlet 2 is *wet* iff H > s2.
Status labels: **[R]** rigorous (under the stated assumptions), **[A]** asymptotic estimate, **[F]** fitted or conventional rule.

## 1. Fixed geometry, no heat
1. **Activation threshold [R].** At steady inflow I the stage is above s2 iff I > I* := Q1(Δ). Proof: steady stage solves I = Q1(u+Δ) + Q2(u), strictly increasing in u ≥ 0 and equals I* at u = 0.
2. **Comparison theorem [R].** If I(t) ≥ I* for all t and H ≥ s2 at t0, then H ≥ s2 for all t ≥ t0 (A dH/dt ≥ I* − Q1(Δ) = 0 on H = s2). If I(t) ≤ I* for all t and H(t0) ≤ s2 then H ≤ s2 afterwards. Inflow that stays on one side of I* never changes the state of outlet 2 — this is what "guaranteed" means in the code.
3. **Drying bound [R].** For H ≥ s2, Q1 + Q2 ≥ I* so `A du/dt ≤ I − I*` with u = H − s2. Hence u(t1) ≤ u(t0) − D(t0,t1)/A, D = ∫(I* − I)dt. If D(t0,t1) ≥ A u(t0) the outlet **has dried by t1**. (Sufficient condition for a dry period; the mean inflow does not enter.)
4. **Wetting bound [R].** If the total rating slope is bounded by G on 0 ≤ u ≤ u_max, `A du/dt ≥ I − I* − G u`, a linear differential inequality; solving it gives a sufficient condition for staying wet and for the first-activation time (an upper bound).
5. **First activation time, exact [R].** From a dry state H < s2 only outlet 1 drains: `t_act = inf{t : H1(t) = s2}` with `A dH1/dt = I − Q1(H1 − s1)`; computed by quadrature `∫ A dH / (I − Q1)` for constant I, by one scalar ODE otherwise. Relaxation time at the threshold: τ_c = A/g1, g1 = Q1'(Δ) — the time scale below which inflow fluctuations are filtered out **[A, linearisation]**.
6. **Mean is insufficient [R, counterexample].** Take inflow J1 ≡ Ī > I* and J2 = Ī + δ for 0 < t < T, Ī − δ(…) such that ∫J2 = ∫J1 and J2 < I* on an interval of length ℓ with D ≥ A u. Same water and heat budget, one run stays wet and the other dries (§3 bound). In the development set D03: mean inflow below I* by a factor 0.84 while the series exceeds I* for 3.8 h; D02–D05 deliver 2–10× more outlet-2 volume than the averaged run although some wet hours coincide. For a *scalar threshold* the information needed is therefore not the mean but the occupancy measure of I about I* at scale ≥ τ_c; this is the Bogolyubov/Filippov statement, not a new result.
7. **Screen R [A/F].** `R = max_t CUSUM(wrong-side inflow volume) / (A |u_ss(Ī)|)`: ratio of the largest excursion volume on the wrong side of I* (reflected cumulative sum of I* − I if the mean says wet, of I − I* if it says dry) to the water volume separating the steady stage at mean inflow from sill 2. R ≥ 1 means the fluctuations can move the stage across s2 although the mean says otherwise; by item 3 this is *rigorous only as a sufficient condition for a dry period when the initial stage is the steady one*. The cutoff 1 is a natural-units choice, not a tuned number; it is [F] as a classifier.

## 2. Evolving geometry and thermal storage (extension; weaker)
* Sill lowering changes the threshold: `I*(t) = Q1(s2(t) − s1(t))`, `dI*/dt = g1 (ds2/dt − ds1/dt) ≈ g1 (E1 − E2)`-type bound **[R, kinematics]**; the frozen fixed-geometry model SE therefore errs wherever the sills move appreciably during the event (dev: SE wet hours 13.0 h vs 3.7 h in D03, 23.1 h vs 7.8 h in D05).
* Lake temperature is an exponential filter of the inflow temperature with flushing time V/Q **[R]**; sill lowering rates are ∝ T_L only while the outlet is wet, so heat reaches outlet 2 only during wet periods: the wetness–temperature covariance is not recoverable from mean enthalpy **[A]**.
* No smooth theory is applied across wet/dry transitions: the quantities above are defined piecewise (each outlet either rated or zero), and the reduced model is integrated with event-bounded blocks, not with a smoothed Q2.

## 3. What must be retained beyond mean water and heat
(i) inflow volume and enthalpy per run above/below I* at resolution ≥ τ_c (blocks), (ii) their order (the stage is a reflected quantity: the initial stage and the excursion order matter), (iii) the geometry-drift rate on the same blocks. Mean water and mean heat alone cannot determine the event times (counterexample 6).

## 4. Validity limits
Quasi-steady ratings (τ_c ≫ travel time), well-mixed lake, two outlets, temperate ice, rating tables at initial geometry (thermal lowering enters only through the reduced model), 36 h horizon, hourly resolved forcing. No claim for real glaciers.
