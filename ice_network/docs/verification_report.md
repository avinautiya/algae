# Verification report (kernels only; no adaptive topology, no neural network, no field validation)

All numbers below are from files in `records/` produced by the scripts in `examples/` and `tests/` (30 tests pass, 175 s). Three kinds of evidence are kept apart: **solution accuracy** (agreement with independent references and convergence under refinement), **conservation** (bookkeeping identities), and **solver equivalence** (the faster discharge solver against the unchanged one).

## 1. The replacement discharge calculation
**Governing equation (one edge, one direction).** For discharge Q find R(Q) = H_req(Q) - H_a = 0 where H_req(Q) = E_0(Q) + k_e v_0(Q)^2/(2g) is the reservoir stage that supplies the entrance of the steady gradually varied profile of discharge Q computed by the standard-step method on total head E = z + h + Q^2/(2 g w^2 h^2), E_k = E_(k+1) + 0.5 (Sf_k + Sf_(k+1)) dx, Manning Sf, marching upstream from h_L = max(H_b - z_(N-1), h_c). If R(Q_in) <= 0 the entrance is critical-flow controlled and Q = Q_in = w_0 sqrt(g) h_c^1.5, h_c = (H_a - z_0)/(1.5 + 0.5 k_e).
**What changed.** (i) Station depth: bisection (80 iterations) replaced by safeguarded Newton on g(h) = z + h + Q^2/(2gw^2h^2) - E_(k+1) - 0.5(Sf(h) + Sf_(k+1))dx with analytic g' (bracket [h_c, h_hi] kept; bisection fallback), converged at |dh| <= 1e-14 h. (ii) Discharge root: bisection (70 iterations) replaced by Brent's method on R(Q), bracket [lower, Q_in] where the lower end is found by halving Q_in until R < 0, converged at bracket <= 2e-15 Q or |R| <= 1e-14 (1 + H_a - z_0) m. The equations, the friction law, the entrance loss and the control-point rule are unchanged; the old code is kept as `hydraulics/channel_ref.py` and is the reference below.
**Non-uniqueness (disclosed).** As Q -> 0 at a free outfall the brink critical depth becomes micrometric, the Manning slope there is >> 1, and R(Q) is not monotone; the root is unique only for Q above about 1e-8 m3/s. Below that the two solvers may select different roots (largest absolute difference 7e-11 m3/s). Flows < 2^-40 Q_in are reported as 0. The effect on a water budget is < 1e-6 m3 per 10^4 s.

| Measure | Result |
|---|---|
| Cases | 645 (3 geometries, 43 levels from 1e-8 to 2 m above the sill, 5 tailwater depths incl. free outfall) |
| Max relative difference fast vs reference, Q > 1e-8 m3/s (575 cases) | 2.1e-10 (median 1.5e-13) |
| Brent iterations | mean 6.4, max 22 |
| Max entrance-equation residual at the returned Q | 1.4e-14 m |
| Cost per call | 9.7x faster (2.1 ms vs 17.9 ms for the 9-station, msub=16 channel) |
| Coupled 3 h melting run, fast vs reference, rtol 1e-10 | sill lowering 5.6e-13, delivered volume 4.7e-12, max outlet Q 5.4e-12, max storage 4.5e-13 relative difference; 21.5x faster (0.6 s vs 12.9 s). **The converged physical solution is unchanged to ~1e-11, far below the discretisation errors in section 4.** |

## 2. Solution accuracy of the fixed-geometry channel (independent references)
* Uniform flow: Q equals Manning's normal discharge to 1e-9 (1e-13 observed), depth spread 1e-14.
* Independent depth-ODE reference (DOP853 on dh/dx = (S0 - Sf)/(1 - Fr^2) with the discharge found by shooting; different formulation from the head-based solver; shares Manning and the entrance-loss closure): free outfall, submerged tailwater, flat ponded near-critical channel agree to 1e-6, 1e-6, 5e-3 (observed 6.8e-11, 7.5e-9 and 2.4e-4 at msub = 64; at the default msub = 16 the near-critical ponded case is off by 0.26 %, and with msub = 1 it is off by 96 % — the substep refinement is not optional there).
* Substep refinement converges at second order (observed 2.0).
* Weir law: a steep, nearly frictionless channel is entrance-controlled and reproduces Q = w sqrt(g) (h/(1.5 + 0.5 k_e))^1.5; a prismatic basin draining through it matches the analytic h(t) = (h0^-1/2 + c t/2)^-2 to 2e-6 over 1 h.
* Backwater: the tailwater has no effect while the exit is a free outfall (asserted to 1e-9), and reduces Q monotonically once submerged, agreeing with the depth-ODE reference to 2e-3.

## 3. Behaviour across levels, near-dry, reversal, symmetry
* Near dry: Q is non-negative, non-decreasing in head from 1e-9 m upward, zero below the friction-limited onset, smooth power-law-like onset (log-log slope between 0.9 and 2.2).
* Reversal: for a flat symmetric channel Q(dH) = -Q(-dH) to 1e-9, Q = 0 at dH = 0, monotone through zero with Q ~ |dH|^(1/2) (fitted exponent 0.4-0.6). In the network, two storage nodes exchange water in the correct direction, the flow reverses when the lower node is filled past the other, levels equalise to 1e-3 m, and total volume equals initial plus input (1e-9).
* Symmetric branches: identical channels deliver identical volumes (relative 1e-9) and the solver returns bitwise identical discharge for identical inputs at ten levels.
* Branching junction: the steady stage and both branch discharges equal the independently root-found solution of Q1(H) + Q2(H) = I (stage to 1e-6 m, discharges to 1e-4 after 24 h), total delivered over input 0.80 at that horizon.

## 4. Refinement (coupled melting case, `records/convergence_study.json`, `convergence_orders.json`)
Case: lake 2000 m2, inflow 0.05 m3/s at 1 C, 100 m channel 0.5 m wide, 6 h. Errors are relative to the Richardson-extrapolated value.

| Refinement | Observed order | Default-setting error (sill lowering / delivered / final storage / outlet Q at 6 h) |
|---|---|---|
| Geometry update step dt = 3600...225 s | 0.92-0.97 (sill), 1.03-1.09 (delivered, storage, Q) | dt = 900 s: 6.4e-4 / 3.1e-3 / 8.1e-4 / 8.7e-3 |
| Stations n = 5...33 | 0.98-1.02 (sill, a point quantity), 1.2-1.6 (others) | n = 9: 2.5e-2 / 1.7e-3 / 4.5e-4 / 1.4e-3 |
| Substeps msub = 4...32 | 2.0 | msub = 16: 3.5e-7 / 8.1e-8 / 1.0e-7 / 3.9e-7 |
| ODE rtol 1e-6, 1e-8 vs 1e-10 | convergent | 5e-7, 4e-9 (sill) |
The coupled solution is therefore **first-order in the geometry-update step and in the station spacing** (operator splitting and segment-mean melt rates), second order in the hydraulic substeps. The default dt = 900 s, n = 9 gives about 1 % error in outlet discharge and 2.5 % in sill lowering for this case; production use would need either smaller steps or a higher-order splitting, not yet implemented.

## 5. Mass: ice lost and liquid water gained (case A, 6 h, `records/budget_tables.json`)
Ice lost 3466.05 kg from the heat delivered divided by L_f; 3466.047 kg from the change of channel geometry times rho_i (3.7798 m3 of ice); liquid water gained 3466.05 kg = 3.4660 m3 with rho_w = 1000 (volume ratio 0.9170 = rho_i/rho_w). Geometry-implied versus energy-implied ice volume: ratio 0.99999993.

## 6. Dissipative heating: exactly how it enters
* Mechanical energy dissipated between node stages: Phi_tot = rho_w g |Q| (H_a - H_b), split into the entrance loss Phi_in = rho_w |Q| k_e v_0^2/2, friction Phi_k = rho_w g |Q| hf_k per segment (the same head losses that set the discharge) and the exit remainder Phi_exit. Phi_in and Phi_k heat the **channel water** (source term in its energy equation); Phi_exit heats the **downstream node water**. Melting is driven only by the water-to-ice heat flux h_T (T - T_m); there is no separate dissipation term in the incision, so nothing is counted twice.
* The switch `heat_dissipation=False` sets Phi_in, Phi_k, Phi_exit to zero in the thermal problem and in the energy budget; the friction law, head losses and hence discharge and storage are untouched. Tests: with fixed geometry the storage and discharge series are equal to integrator tolerance (1e-8; not bitwise, because the extra temperature states change the step sequence); the difference of (heat leaving + heat stored) between the two runs equals the cumulative dissipation to 1e-8; steady outlet heat flux equals rho_w g Q (H_lake - H_outlet) to 2e-3; with inflow at 0 C the heat delivered to ice (1.42 MJ) never exceeds the channel share of the dissipation (8.0 MJ), and the exit share (35.5 MJ) goes to the outlet water.
* Case A vs B (dissipation counted vs excluded, 1 C inflow): latent heat to ice 1158 vs 1156 MJ, sill lowering 0.054071 vs 0.054062 m, i.e. dissipation changes the incision of this warm case by about 0.2 %.

## 7. Conservation on nonzero-throughput cases (separate from accuracy)
Relative residuals, normalised by the sum of the absolute terms: water 4.5e-16, mass 1.7e-15, energy 7.6e-16 in case A (delivered/input 0.80, Q from 6 to 55 L/s); the same order (<= 3e-15) in all five cases and in the refinement study. Absolute energy residual 2.7e-5 J against scale 3.5e10 J. **These identities are linear invariants of the discretisation and would hold for a wrong physical model; they are not evidence of correct physics.** Accuracy is established only in sections 1-4.

## 8. Not done / limitations
No SaDS or SHED run; no resolved-flow (CFD) or laboratory comparison of the heat-transfer correlation (Dittus-Boelter applied to open channels, factor 2 sensitivity not yet quantified); no supercritical reaches or hydraulic jumps; melting has not been validated against any experiment, so the melt-rate coefficients are unvalidated; first-order time splitting; hydraulics tested only on rectangular sections and Manning friction with n = 0.04 (assumed, not measured); node temperatures well mixed. The kg and energy bookkeeping show consistency, not agreement with ice-channel observations.
