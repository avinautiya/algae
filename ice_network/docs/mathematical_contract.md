# Mathematical contract (version 1, as implemented)

SI units. Signs: edge discharge Q > 0 flows from the edge's node `a` to node `b`; heat to ice q >= 0; temperatures are relative to the ice melting point T_m = 0 C (temperate ice). This describes only what is coded and tested; anything else is out of scope or unimplemented.

## 1. Domain and exclusions
Freshwater in open channels (edges) and depressions (nodes) on bare ice. Excluded: subglacial/pressurised flow, hydrofracture, calving, sediment, salt, ice dynamics, surface ablation (the terrain `z_s` is fixed), downstream inundation, snow, refreezing, neural networks, topology change, any adaptive or biological rule. Candidate connections are a fixed, hand-declared graph in this milestone.

## 2. Nodes (storage)
Lake node i: stage H_i = z_b,i + h_i with analytic volume V_i = A_ref h^(q+1) / ((q+1) h_ref^q) (q = 0: prismatic). Outlet node: fixed stage H_o (a free-outfall sink; asserted below the exit bed of every edge that reaches it, so it can never supply flow). Water balance with mass in kilograms (rho_w = 1000 kg m-3, rho_i = 917 kg m-3):

dV_i/dt = I_i + sum_in (Q + Q_melt) - sum_out Q,     Q_melt,e = (sum_k q_ice,k) / (rho_w L_f)

Q_melt is the meltwater produced in edge e. It is added **only at the downstream node** (never also on the edge), the channel's own water storage is neglected, and the hydraulic discharge along an edge is taken constant (Q_melt/Q is reported). In kilograms: rho_w Q_melt = (melted ice mass rate) = q_ice / L_f = rho_i dV_ice/dt. Four distinct quantities are tracked: external input I; ice melted in channels (latent heat); water stored in nodes; water leaving through outlets.

## 3. Edges (quasi-steady open-channel hydraulics)
Rectangular section, width w_k, bed z_k at stations k = 0..N-1; Manning n; segments subdivided into `msub` substeps (geometry interpolated linearly).
* Regime: quasi-steady gradually varied flow. Valid when channel storage << node storage and the wave transit time L/sqrt(g h) is short compared with the time scale of node-stage change (true here: seconds against hours; **not checked automatically**).
* Standard-step on total head E = z + h + Q^2/(2 g w^2 h^2), marching upstream from the exit depth h_L = max(H_b - z_(N-1), h_c): E_k = E_(k+1) + 0.5 (Sf_k + Sf_(k+1)) dx with each station's subcritical depth found by bisection (not by fixed-point iteration, which oscillates near critical flow and made the discharge law noisy). Where no subcritical depth exists, depth is critical and the energy line is raised to z + 1.5 h_c (control point: weir, hump, brink); that raise is not friction.
* Entrance: H_a = E_0 + k_e v_0^2/(2 g). Capacity Q = min(Q_outlet, Q_inlet), Q_inlet = w_0 sqrt(g) h_c^1.5 with h_c = (H_a - z_0)/(1.5 + 0.5 k_e). Dry or humped channels give Q = 0. Direction follows the water-surface gradient; the reverse direction uses the reversed station arrays; Q = 0 at H_a = H_b.
* Wetting/drying: an edge carries water only if the upstream stage exceeds the highest bed point. Hydraulic deactivation is not removal: geometry persists.
* Momentum: not solved (no jumps, waves, surges). Supercritical reaches are not represented: stations flagged `choked` have critical depth, which overestimates depth and underestimates velocity there.
* Junctions: nodes are storage tanks; mass is conserved by construction (each edge adds Q to one node and removes it from another).

## 4. Dissipation and heat (explicit terms)
Mechanical energy lost between node stages: Phi_tot = rho_w g |Q| (H_a - H_b) (W) = Phi_in + sum_k Phi_k + Phi_exit, where Phi_in = rho_w |Q| k_e v_0^2/2 (entrance), Phi_k = rho_w g |Q| hf_k (friction head loss per segment), and Phi_exit is the remainder (jump/fall/expansion at the exit) which is deposited in the downstream node's water. Phi_in and Phi_k heat the channel water; they are the only dissipative heat source and **are not also fed into an incision law** (no double counting).

Water temperature along a segment (constant coefficients within it), phi' = Phi_seg/dx:
rho_w c_w Q dT/dx = phi' - h_T P (T - T_m); h_T = Nu k_w / D_h, Nu = 0.023 Re^0.8 Pr^0.3 (Dittus-Boelter, water being cooled; smooth-pipe correlation applied with D_h = 4A/P to an open rectangular channel — a stated approximation, sensitivity factor `nu_factor`).
Exact segment solution: T_(k+1) = T_eq + (T_k - T_eq) exp(-h_T P dx/(rho_w c_w Q)), T_eq = T_m + phi'/(h_T P).
Heat to ice in the segment: q_seg = rho_w c_w Q (T_k - T_(k+1)) + Phi_seg (>= 0). Interface condition: rho_i L_f v_m = h_T (T - T_m) (equivalently q_seg = rho_i L_f dV_ice/dt over the segment).
Nodes are well mixed: d(V T)/dt = sum_in [Q T_out + Phi_exit/(rho_w c_w)] - sum_out Q T_i + I T_I (meltwater at T_m carries no sensible heat).

**Included energy terms:** transported sensible heat rho_w c_w Q T; latent heat rho_i L_f of melting; dissipative heating (Phi_in, Phi_k, Phi_exit); external inflow enthalpy. **Excluded:** conduction in ice, cold-ice warming, refreezing, solar/radiative/turbulent exchange of channel or lake water with the atmosphere, kinetic energy of the water leaving at outlets, sensible heat of ice, viscous heating inside nodes, heat of lake-stored water other than as a well-mixed tank. Potential energy of stored water is a mechanical quantity and is not part of the thermal budget.

## 5. Geometry update (operator split)
Within a macro-step `dt_geom` the geometry is frozen and the coupled ODE (V, V T, outlet totals, cumulative heat per segment, wall-weighted heat) is integrated with LSODA at tight tolerance. Then, per segment, the melted ice volume dV = Heat/(rho_i L_f) is divided by the wetted-perimeter share of the heat (uniform h_T and T over the wetted perimeter, so equal melt speed normal to bed and wetted walls: bed share w/P, wall share 2h/P — derived, not assumed). Bed lowering dz = dV_bed/(w dx); wall retreat is spread over the full bank height D' = D + dz (a rectangular-section idealisation: the undercut band is assumed to relax to a vertical wall), dw = dV_wall/(D' dx), so the section area changes by exactly w dz + D' dw. Station values are averages of adjacent segments; the residual disagreement between the volume implied by the station geometry and the energy-based volume (the `geometry_over_energy` ratio) is reported, not clipped. A sill incised below its lake bottom stops the run (the lake would drain completely; outside scope).

## 6. Budgets (defined before the tolerances in the tests were chosen)
Water: sum V(t) + sum C(t) - sum V(0) - W_in - Q_m = 0, normalised by sum of the five absolute terms. Mass: the same multiplied by rho_w, plus the check rho_w Q_m = (heat to ice)/L_f = rho_i V_ice,energy. Energy: rho_w c_w sum S(t) + E_out + H_ice - [rho_w c_w sum S(0) + E_in + Phi_cum] = 0, normalised by the sum of absolute terms; absolute residuals are stored with their scales. Because extensive quantities are integrated as extra ODE states, these identities are linear invariants of the discretisation; their residuals test **bookkeeping and the integrator**, not the physics. Test thresholds (1e-8 to 1e-9 relative) were fixed relative to the requested ODE tolerance (rtol 1e-9) before the refinement study was run.

## 7. What is borrowed and what is new
Borrowed: the structure "edge capacity adapts to flow through a coupled conservation law" (Tero et al. 2008 [not read]); here the adaptation law is the physical melt law of sections 4-5, not a cost-minimisation rule. Not used: Murray's law, shortest-path or growth/pruning rules, a decay term on conductivity. New in this code (to be judged against SaDS/SHED, not claimed as novel): advected-sensible-heat melt, wetted-perimeter melt partition, level-coupled GVF edges in a small graph, and the audit.
