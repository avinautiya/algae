# A lake feeding two evolving ice channels: capture kinematics, quasi-steady theory, and the reduced model RTM

Scope: idealised controlled problems (a lake of area A, inflow I, two outlet channels on bare temperate ice). No field data, no real-glacier forecast. Everything below is derived for the physical assumptions of `docs/mathematical_contract.md` (Manning hydraulics with critical-depth control, Dittus-Boelter heat transfer, temperate ice, well-mixed lake at the inflow temperature T_L) and is checked numerically in `VERDICT_CAPTURE.md`.

## 1. Model in the notation of this note
Lake stage H (prismatic): A dH/dt = I(t) - Q1 - Q2. Outlet i has a quasi-steady rating Q_i = R_i(H, G_i(t)) with G_i its channel geometry (bed profile z, width w). Define the head over the (effective) sill u_i = H - s_i, the **rating sensitivities** g_i = dQ_i/dH, lambda_i = g_i/Q_i, and the **rating-effective incision rate**

  E_i := (dQ_i/dt at fixed H) / g_i   [m/s]

(the rate at which a uniform lowering of the sill would have to proceed to change the capacity as the real, non-uniform melting of bed and widening of walls does). With the heat budget of the model: the water entering route i at T_L loses heat along the reach, T(x) = T_L exp(-a_i x / L_i) in the heat-limited form with a_i = h_T P L / (rho_w c_w Q_i); melt speed normal to the wetted perimeter is h_T T/(rho_i L_f), bed lowering rate = that speed at the bed share w/P of the heat, walls retreat with the rest. E_i is the response of the production kernel (computed by perturbing the geometry with the heat of the current state, no assumption on its shape).

## 2. Exact kinematics (no closure for E_i)
(K1) dQ_i/dt = g_i (dH/dt + E_i). (definition of g_i, E_i)
(K2) d ln(Q2/Q1)/dt = (lambda2 - lambda1) dH/dt + (r2 - r1), with r_i = lambda_i E_i the fractional capacity growth rate. The first term is the hydraulic (lake-level) re-partition; the second is the incision-driven one: with phi2 = Q2/(Q1+Q2), phi2' = phi1 phi2 [ (lambda2-lambda1) H' + (r2-r1) ] is a replicator-type equation whose "fitness" is the capacity growth rate.
(K3) Let Delta = s2 - s1 (elevation of sill 2 above sill 1, the *sill gap*) and require that each rating depend on the geometry only through the effective sill (translation invariance, the RTM approximation, section 5). Then d(Delta)/dt = E1 - E2 exactly: the gap evolves by the difference of the two effective incision rates, **independently of the lake and of the inflow history except through the E_i** (u_i enters E_i).
The statement "an outlet gains flow iff it lowers faster" is (K3) plus the monotonicity of the ratings. It is elementary; its use is to isolate E_i as the single physical input.

## 3. Quasi-steady lake (tau_lake << tau_incision)
If the lake relaxes quickly, A H' is slaved: Q1 + Q2 = I and H' = -(g1 E1 + g2 E2)/(g1 + g2). Then phi2' = Lambda (E2 - E1)/I with Lambda = g1 g2/(g1+g2) > 0, and the sill gap obeys the **autonomous phase line**

  d(Delta)/dt = G(Delta) := E1(Q1(Delta)) - E2(Q2(Delta)),

where Q_i(Delta) is the partition at gap Delta (solve R1(u1) + R2(u1 - Delta) = I). Consequences, all computed numerically from the production kernels (`rtm.gap_dynamics`):
* kappa: gap at which Q1 = Q2 (equal flows). Outlet 2 can overtake outlet 1 only if Delta0 > kappa and G < 0 on [kappa, Delta0].
* fixed points G(Delta_c) = 0: **capture threshold** Delta_c. At Delta_c both outlets incise equally although one carries less flow. d G/d Delta > 0 at Delta_c (unstable): for Delta0 < Delta_c outlet 2 gains share and eventually captures; for Delta0 > Delta_c outlet 1 captures and outlet 2 is abandoned. Delta_c depends on the hydraulic and heat-transfer properties of the two routes, **not on the inflow temperature** to first order (every E_i is proportional to T_L; the heat-limitation factor a_i depends on Q but not on T_L), and not on lake size.
* symmetric pair (identical routes): Delta = 0 is a fixed point with growth rate sigma = dG/dDelta|0 = dE/du (u: head over the sill): the sill gap grows exponentially, e-folding 1/sigma. A perfectly symmetric pair is unstable to capture by the lower sill under any incision law whose rate increases with head.
* the replicator form (K2) adds the hydraulic term: with the lake fixed (H' = 0) the symmetric pair is stable iff d ln(lambda E)/d ln Q < 0, i.e. for power-law ratings Q ~ u^gamma iff d ln E/d ln Q < 1/gamma. Numerically (`rtm_index.json`) thermal incision gives d ln(lambda E)/d ln Q = -0.3 ... -0.5 (stabilising) and a dissipation-only incision law +0.12 ... +0.45 (destabilising) for the same hydraulics: which regime applies depends on whether the lake level is free to follow the sills (section 3) or held (the opposite limit); real lakes lie between, and the combination is what RTM integrates.

## 4. What is and is not a bifurcation, hysteresis or multistability here
* **Ordinary memory** (reversible, relaxing): lake storage; relaxation time tau_lake = A/(g1+g2). Not hysteresis.
* **Irreversible memory**: sills only lower (no refreezing, no deposition in the model): E_i >= 0. Any capture is a ratchet: if forcing later drops, outlet 1's sill stays low; reactivation of outlet 2 requires a stage above its (unchanged or lowered) sill.
* **Bistability in the sense used here**: for constant forcing and fixed parameters the phase line of Delta has two stable regimes (Delta -> -infinity, outlet 2 captures; Delta -> +infinity with outlet 2 inactive, outlet 1 captures) separated by the unstable fixed point Delta_c (a *basin boundary*). This is a statement about the autonomous quasi-steady phase line, not about equilibria of the full system (there are none while incision continues).
* **Not claimed**: a bifurcation in a parameter (Delta_c moves smoothly with the parameters; the only non-smooth events are activation/deactivation of an outlet at u_i = 0), hysteresis under a quasi-static parameter cycle (not demonstrated), or any statement under time-varying forcing beyond RTM numerics.

## 5. Reduced translation-invariant model (RTM)
State (H, d1, d2): d_i the cumulative effective lowering. Q_i = R_i(H - s_i0 + d_i), d_i' = E_i(u_i), with R_i, g_i, E_i tabulated once (70 head values) from the production kernels at the initial geometry (`rtm.tabulate`; the sill-station lowering m_sill_i(u) is tabulated alongside so that RTM also predicts the reported "sill lowering"). Approximations, with validity conditions:
 (a) geometry evolves as a translation of the rating (changes of the *shape* of the bed/width profile are neglected beyond the first-order response contained in E_i). No a-priori smallness criterion is derived: the accumulated incision over 12 h is 0.25-1.5 of the head u over the sill in the development configurations (3.3 in the steep-short symmetric stress case T14; `rtm_index.json`), so validity is established only empirically, by the comparison with the network, which evolves the full non-uniform geometry;
 (b) E_i(u) is frozen at the initial geometry: valid while relative changes of h_T and of the exchange number a_i are small;
 (c) quasi-steady hydraulics and the production closures (Manning, critical-depth control, Dittus-Boelter x factor 1, temperate ice): no momentum, no supercritical reaches, no refreezing;
 (d) the lake is the prismatic well-mixed tank of the contract; meltwater input to the lake is neglected (an addition at the downstream node in the production model).
Cost: one tabulation (about 1.5 s per route and treatment) plus a three-variable ODE (milliseconds).

## 6. Closed-form scalings (hydraulic closures)
With h_T = nu_f 0.023 k Pr^0.3 (v D_h/nu_w)^0.8 / D_h and melt speed m = h_T T(x)/(rho_i L_f):
* at an entrance-controlled (critical) section, v ~ Q^(1/3), D_h ~ 4 h_c ~ Q^(2/3) so m ~ Q^(2/15): alpha := d ln m/d ln Q = 0.13; with the exchange factor (1 - e^-a)/a and widening the numerical alpha is 0.23-0.52 (`rtm_index.json`);
* uniform Manning flow in a wide channel: v ~ Q^(2/5), D_h ~ Q^(3/5): alpha = 0.2, and the melt speed scales with S^0.18, so steeper routes incise slightly faster at equal Q;
* dissipation-driven melt (all dissipation melts the wall): m = rho_w g Q S_f/(rho_i L_f P) ~ Q: alpha ~ 1.
In the tested geometries the thermal law incises about 300 times faster than the dissipation-only law at 1 C (inlet-sill lowering 75 mm against 0.2 mm in 12 h, symmetric case D01), and depends more weakly on Q.
