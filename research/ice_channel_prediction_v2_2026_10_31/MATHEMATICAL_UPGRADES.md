# Prediction-oriented mathematical checks

These derivations complement the earlier ice-step packet. They are standard conservation/closure arguments applied to the proposed model, not certified new mathematical results. A novel result still requires solving and testing the particular coupled problem.

## 1. Why a bulk temperature does not uniquely specify melting

This constructed example rules out universal exactness of a purely instantaneous bulk closure. It does not show that such a closure is inaccurate in every practical regime.

Consider a flat water layer, normalized height 0≤y≤1, with fixed velocity shape w(y)=2y−y². It is no-slip at the ice boundary and has zero shear at the top. Let

    f(y)=y(1−y)²(1−13y/2+8y²),
    theta_±(y)=theta_s [w(y) ± f(y)/2], theta_s>0.

Both profiles have identical boundary temperatures: zero excess at the ice and theta_s at the top. Their top conductive fluxes are also identical because f'(1)=w'(1)=0. Direct polynomial integration gives

    integral_0^1 f dy=0,
    integral_0^1 w f dy=0.

Thus both depth-mean temperatures equal 2 theta_s/3, and both flow-weighted temperatures equal 4 theta_s/5. Depth, velocity, flux, bed geometry, bulk temperatures, top boundary temperature and top conductive flux all match.

The physical ice-side water gradients, however, are

    d theta_−/dz at bed = (3/2) theta_s/h,
    d theta_+/dz at bed = (5/2) theta_s/h.

The molecular wall heat fluxes differ by a factor 5/3, even though all those bulk descriptors match. For an exactly resolved no-slip thermal wall, molecular gradient determines the wall flux; turbulent motion affects the profile feeding that gradient.

The profiles are nonnegative: with P=1−13y/2+8y², |P|≤5/2 on [0,1], and |f/w|≤5/4. Therefore |f/2|≤5w/8 and theta_±≥3 theta_s w/8. These are admissible smooth instantaneous temperature fields with the stated boundaries.

For any one deterministic prediction made from the shared descriptors, the maximum absolute error over these two constructed states is at least half their flux difference. This is an information limitation, not evidence of a newly discovered physical effect. Extra spatial/thermal state or a restrictive equilibration assumption is necessary for an exact closure on this full state class.

The practical research test is whether physically generated upstream histories produce a persistent, consequential version of this ambiguity at environmental Reynolds/Péclet numbers. They may not: rapid turbulent relaxation could make the bulk closure adequate. That outcome would eliminate a proposed memory extension, not invalidate conservation laws.

## 2. Conservative exchange between CFD and a reach model

For each fixed reporting control volume, use watts rather than mixing W/m, W/m² and total power. Let E_w be water sensible energy relative to the interface melting temperature. Define inlet and outlet enthalpy fluxes using flow-weighted temperature:

    P_in/out = integral_in/out rho_w c_p (T−T_m) u·n_flow dA.

Let P_a be net atmospheric heat into water, P_d dissipation heating actually retained in the thermal model, P_w heat from water into the ice interface, P_i conduction from the interface into cold ice, and P_m latent melting power. With the relevant moving-boundary/mass terms consistently handled,

    dE_w/dt = P_in−P_out+P_a+P_d−P_w,
    P_w=P_i+P_m,
    P_m=rho_i L_f integral_Gamma m_n dA.

Hence

    P_m=P_in−P_out+P_a+P_d−dE_w/dt−P_i.

This makes total melting and downstream cooling jointly testable. A model cannot increase melting without paying for it through supplied heat, storage or reduced export. No heat-loss multiplier should be applied independently to both water cooling and ice recession.

For discrete coupling, compute every water–ice face exchange once. Insert exactly opposite copies into the water and interface/ice ledgers. Integrate resolved CFD fluxes over the actual wetted area before projecting to a coarse reach cell. Preserve the integrated power on coarsening; retain shape information separately if needed to predict its distribution. Verify the projection on nonuniform meshes, moving wetted areas and remapping. The simple identity checked in prediction_basis.py is not a test of a full discretization.

Do not equate numerical momentum diffusion with physical latent heating. A mechanical-energy-consistent treatment must account for physical dissipation, bed work, phase-change mass and potential energy. If negligible mechanical heating is omitted, bound its contribution for the case and label that approximation.

## 3. Phase-change mass and moving-wall conditions

Let n point from ice into water, ice material velocity u_i≈0, and recession m_n>0. The interface velocity is V_Gamma·n=−m_n. Mass continuity requires

    rho_w (u_w−V_Gamma)·n = rho_i (u_i−V_Gamma)·n,

so

    u_w·n = (rho_i/rho_w−1)m_n.

The fluid normal velocity is not simply the geometrical interface velocity. A moving-wall condition that sets u_w=V_Gamma indiscriminately represents an impermeable moving wall rather than melting ice. Tangential no slip is relative to the ice material, while normal mass transfer follows the equation above. If the solver uses equal densities, explicitly state that approximation and recover the correct equal-density limit u_w·n=0 for stationary ice.

On depth integration, the water-volume source is (rho_i/rho_w)m for small slopes. Ice volume loss and added water volume therefore differ. ALE/VOF remapping and this source must not count the same new volume twice.

## 4. A conditional prediction-error bound

For vertical recession b_t=−q_net/(rho_i L_f), consider a true modeled heat-flux functional F and an approximation F_tilde with identical forcing. Within a fixed admissible regime suppose

    ||F(b)−F(b_tilde)||≤K_q ||b−b_tilde||,
    ||F(b_tilde)−F_tilde(b_tilde)||≤epsilon_q.

The chosen norm, common spatial registration, interfaces and regime must be specified. This scalar bound assumes auxiliary states are prescribed or are already represented consistently by F; a full flow/thermal-memory predictor requires an augmented-state sensitivity bound. It must not inherit this scalar certificate without that derivation. K_q has units W/m³ and epsilon_q has units W/m². Gronwall's inequality gives

    e(t)≤e(0) exp[K_q t/(rho_i L_f)]
          +(epsilon_q/K_q)(exp[K_q t/(rho_i L_f)]−1).

For K_q→0 this becomes e(0)+epsilon_q t/(rho_i L_f). It is a conditional propagation bound, not a global certificate based on a few sampled cases. Moving jumps, dry/wet changes, topology transitions and non-Lipschitz closures may invalidate a uniform K_q. Resolved comparison can estimate local sensitivity but cannot establish a global bound without additional analysis.

Using rho_i=917 kg/m³ and L_f=334000 J/kg, a persistent 100 W/m² latent-flux bias produces 0.02821 m of recession error in 24 hours when feedback is absent. To keep that contribution below 0.02 m/day requires a mean bias below approximately 70.90 W/m². These are illustrative design tolerances, not observational resolution claims or promised model accuracy.

The research value of such a bound is to connect a heat-transfer benchmark to a forecast requirement. It prevents reporting small-looking temperature errors without assessing their accumulated geometry effect.

## 5. Hydraulic causality and a limit on impact claims

The local one-dimensional shallow-water characteristic speeds are u−sqrt(g h) and u+sqrt(g h). With downstream-positive flow and Fr>1, both are downstream. Infinitesimal hydraulic changes in that region do not travel upstream through it. A melt hotspot below an isolated hydraulic control does not automatically change lake release.

To test a lake-drainage impact, couple lake storage A_L(H) Hdot=I−Q_out to a compatible hydraulic boundary and model the reach between the lake and the changing control. Do not prescribe Q_out as an input and then claim its predicted change is an independent result. Finite disturbances, migration of a control and evolving geometry require explicit analysis; the local characteristic test is a necessary diagnostic, not a global no-feedback theorem.

## 6. Scale checks that change the implementation

Use explicit definitions:

    Fr=U/sqrt(g h), Re_h=U h/nu_w,
    We_h=rho_w U² h/sigma, Bo_h=rho_w g h²/sigma,
    Ste=c_p Delta T/L_f.

Some literature defines Stefan number as its reciprocal. Save both definitions in metadata and never compare their numerical values without conversion. Distinguish hydraulic supercriticality (Fr>1) from a ripple-instability threshold in Reynolds number.

For the illustrative published flume hydraulics h=0.00165 m and U≈0.403 m/s, with nominal material constants, Re_h≈369, We_h≈3.57 and Bo_h≈0.356. This does not prove capillarity dominates the long bedform wavelength. It does show that small-depth/free-surface curvature near a jump cannot be dismissed solely because the bedform itself is long. The resolved reference should retain physical surface tension or demonstrate its insignificance at the relevant local scales.

For a synthetic h=0.05 m, U=1 m/s reach, Fr≈1.43 and Re_h≈27778, putting it in a different flow regime. A laminar match at the flume scale cannot be silently transferred to that case. These are screening calculations, not field observations.

## 7. Checks completed in this revision

Run `python3 prediction_basis.py`. Exact rational arithmetic verifies the profile moments and boundary constraints. The script checks positivity numerically in addition to the analytical bound, phase-change mass conversion, a manufactured energy-exchange ledger, the zero-feedback error-bound limit and characteristic-direction examples. Outputs are `design_check_results.json` and `temperature_counterexample.csv`.

The script is not a CFD solver. It does not calculate a new stability boundary, train a predictive closure or validate an environmental model. Those remain the research work specified in PROJECT_V2.md and the next handoff.
