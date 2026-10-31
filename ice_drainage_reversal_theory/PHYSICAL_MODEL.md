# Physical model and scope

## State, fluxes, and units

H is lake-surface elevation [m] above a fixed basin datum, b is the persistent outlet's control-sill elevation [m], and w = H - b is its wet head [m]. V(H) is lake volume [m^3], with A(H) = V'(H) > 0 [m^2]. Q = F(w) [m^3/s] is water drawn from the lake into the persistent outlet. The rating is zero for a dry outlet and positive/increasing for w > 0. The bypass discharge u >= 0 has the same units. Proposition 1 has no external lake inflow; Proposition 2 permits the same finite-duration inflow I(t) in both compared trajectories, and Proposition 3 permits a shared prescribed inflow.

The reduced equations without inflow on a wet interval are

\[
A(H)\dot H=-F(H-b)-u,\qquad
\dot b=-E(F(H-b)).                                      \tag{P1}
\]

E is effective lowering of the *hydraulically controlling sill* [m/s], not the average melt rate of an arbitrary reach. The distinction matters: regional melting that does not move that control cannot simply be substituted for E.

The constant-area nonlinear theorem assumes

\[
0<\beta(w)=\frac{A E(F(w))}{F(w)}\le\beta_{\max}<1.       \tag{P2}
\]

This is stable drainage in the ordinary sense that lake drawdown exceeds incision when u = 0. That concept is established prior work. Proportional incision E(Q) = kappa Q has kappa [m^-2] and beta = A kappa. A linear rating F(w) = k w has k [m^2/s]; it is an exactly soluble class or local rating approximation, not a universal glacier law. The nonlinear theorems do not require that linear rating.

## A thermodynamically admissible origin of the incision law

For a well-mixed lake at constant temperature T_L > 0 above the melting point, a persistent control region of fixed effective ice area S_c [m^2], and captured sensible-heat fraction eta in (0,1), let

\[
P_{ice}=\eta\rho_w c_p T_L Q,
\qquad E(Q)=\frac{P_{ice}}{\rho_i L_f S_c}
           =\kappa Q,
\quad \kappa=\frac{\eta\rho_w c_p T_L}{\rho_i L_f S_c}.    \tag{P3}
\]

The dimensionless maximum allowed by the incoming sensible heat is

\[
\beta\le\beta_{energy}
 =\frac{A\rho_w c_p T_L}{\rho_i L_f S_c}.                 \tag{P4}
\]

A proposed parameter must satisfy eta <= 1 and this bound. A desired effect cannot be obtained by choosing kappa independently of the heat supply and control geometry.

The resulting meltwater flux downstream of the control is

\[
Q_m=\frac{\rho_i S_c E}{\rho_w}.
\]

It is added to downstream channel mass, not to lake outflow in (P1). With melt entering at the melting point, the downstream advected sensible enthalpy is P_out = (1-eta) rho_w c_p T_L Q; hence P_out + rho_i L_f S_c E = P_in. Its mixed downstream temperature is P_out/[rho_w c_p (Q+Q_m)], not (1-eta) T_L when Q_m is appreciable. This closes this lumped mass/enthalpy bookkeeping; it does not resolve a melting boundary layer or downstream momentum.

Removing well-mixed water through either outlet preserves lake temperature in this idealized budget: d(rho_w c_p V T_L)/dt = -rho_w c_p T_L (Q+u), and V_dot = -(Q+u). It does not assume a source replenishes the removed heat. Shared inflow at the same T_L also preserves that temperature, with +rho_w c_p T_L I in the energy balance and +I in the volume balance. Inflow at another temperature requires a lake-heat equation and can change kappa, which is outside the constant-kappa extension. The bypass carries its sensible heat away without contributing to the persistent sill's incision. Any dissipation in its own downstream path belongs to another control volume.

Constant eta is an assumption, not a heat-transfer theorem. A more conventional fixed-conductance cooling closure can instead have

\[
E(Q)=c Q(1-e^{-d/Q}),\quad
c=\frac{\rho_w c_p T_L}{\rho_i L_f S_c},\ d>0,
\]

with

\[
E_Q=c[1-(1+d/Q)e^{-d/Q}]>0.
\]

The nonlinear retention theorem applies only on ranges satisfying P2. This closure's beta can cross one as flow changes; selecting a stable range is essential. Neither closure is declared validated for a natural open ice channel.

Equations P3 describe the sensible-enthalpy/latent-heat budget of this lumped thermal control. They do not close the hydraulic mechanical-energy budget. Frictional heating in the control must either be negligible or explicitly added to E and the outgoing energy. The ratio g Delta z/(c_p T_L) is a useful upper scale for gravitational energy relative to incoming sensible heat over a drop Delta z; kinetic-energy changes also require accounting. Neither ratio has been evaluated for a real site here. A large ratio invalidates a sensible-heat-only interpretation; it does not invalidate the abstract theorem if an appropriate total incision law still satisfies its hypotheses.

## A bypass that closes naturally

A higher rock-controlled or otherwise non-eroding outlet can be represented by

\[
u=G((H-h_s)_+),\quad b_0<h_s<H_0,
\]

where G(0)=0 and G(z)>0 for z>0. While it is active, H >= h_s and b <= b0, so w >= h_s-b0 > 0. The persistent flow remains positive. The lake therefore reaches h_s in finite time; the bypass then runs dry. Since H only decreases with no inflow, it cannot reactivate. The theorem thus has an autonomous two-outlet realization with finite rates and no artificial drainage impulse.

This does not model a second incising ice outlet. If both outlets evolve, or the bypass melts or destabilizes the controlling ice, (P1) must be extended and the present conclusion rederived.

## Applicability that cannot be assumed away

- A single persistent hydraulic control and its rating must remain valid. Downstream controls, hydraulic jumps, flow reversal, or a sill passing beneath a controlling basin floor can invalidate F(H-b).
- The lake is mixed and its basin does not melt, deform, or calve in this model. Real temperature stratification, heating, dilution, cold-ice conduction and evolving basin volume require additional equations.
- The incision law is memoryless here. A transverse thermal field, turbulence, hydraulic roughness and evolving wetted area can change E and its derivative. The previous work established that a common multiplicative conductance shape cannot be assumed.
- No hydrofracture, subglacial pressure coupling, sediment transport, ice creep, refreezing or moving free surface is resolved. A drainage-surface reversal is not a GLOF prediction.
- Residual volume concerns a drainage episode after forcing has ceased. Proposition 1 needs zero inflow; Proposition 2 permits common finite-duration inflow with constant kappa. With persistent glacier melt or rainfall there may be no final residual state; use a finite-horizon forecast and state its forcing.
- Larger residual volume is not automatically worse or better for hazard. A bypass may reduce peak flooding while retaining more water. The published successful mitigation case must not be contradicted by an unrelated idealized example.

The mathematical results are conditional predictions for (P1). CFD can later test the adequacy of F and E and the persistence of the sign under resolved physics. It cannot turn these assumptions into site validation without relevant observations.
