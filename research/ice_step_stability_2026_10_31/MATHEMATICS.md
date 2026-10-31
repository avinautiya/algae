# Mathematical foundation for a thermal ice-step stability study

Prepared 31 October 2026. These are derivations and a proposed analysis, not a solved ice-step eigenvalue problem. Standard conservation laws, diffusion operators, shock conditions and spectral methods are credited as established mathematics. No new theorem of general fluid mechanics is claimed.

## 1. The question the mathematics must answer

Given a train of water-carved steps in ice, which perturbations decay, which grow, and which cause neighboring steps to merge? Can a dynamically stable range of spacings be predicted without imposing the largest admissible spacing? Does that range change when heat supply and ice temperature are modeled dynamically?

Separate four questions: existence of a traveling profile; infinitesimal stability of that profile; finite-amplitude evolution from an initial landscape; and observational accuracy. None implies the next automatically. A stable band need not select a unique wavelength.

The nearest methodological precedents are [Balmforth & Vakil (2012)](https://doi.org/10.1017/jfm.2011.555) for erodible beds and [Johnson et al.](https://arxiv.org/abs/1803.03484) for discontinuous hydraulic roll waves. The ice-specific comparison is [Wu & Izumi (2023)](https://doi.org/10.3390/geosciences13050128). See SOURCES.md for access levels and limitations.

## 2. Variables, geometry and declared approximation

Use a unit-width, shallow open channel, horizontal coordinate x and vertical coordinate z. Water depth h>0; depth-averaged velocity u; discharge per width q=hu; bed elevation b; temperature excess theta=T_water−T_m. Positive m is vertical ice recession, b_t=−m, with small bed slopes outside jumps. Fluxes q_w and q_i are positive from water to the interface and from the interface into cold ice, respectively, in W/m² of projected bed area. Atmospheric input q_a is positive into water. K_i is ice conductivity; kappa_i is ice diffusivity; L_f is latent heat.

For a curved interface use normal velocity and actual interface area consistently; do not use a vertical recession rate with an unconverted normal-area flux. Overhangs and detached waterfalls are outside this shallow model.

A mass-consistent starting point is

    h_t + (hu)_x = r m,
    (hu)_t + (hu² + g h²/2)_x = −g h b_x − tau_b/rho_w + r m u_in,
    rho_w c_p [(h theta)_t + (hu theta)_x] = q_a − q_w + Phi,
    rho_i L_f m = q_w − q_i,
    b_t = −m,

where r=rho_i/rho_w, newly melted water has theta_in=0, and u_in is its entering horizontal velocity (approximately zero for stationary ice). Add rain, lateral inflow, evaporation and radiation only with explicit mass/enthalpy conventions. In the primitive water-temperature equation, dilution by zero-temperature meltwater appears automatically on expanding the conservative equation.

The first spectral problem should use the leading-order small-melt-volume approximation r m≈0 in the water equations, while retaining b_t. Verify the small parameter epsilon_m=m L/(h U) and cumulative meltwater contribution for every application. This is an approximation, not exact phase-change mass conservation. The finite-reach reference must either include the mass source or quantify the resulting error.

The displayed momentum equation does not by itself specify an exactly conservative total mechanical-plus-thermal-plus-phase energy scheme on an evolving bed. Derive that discrete ledger before adding frictional or shock heating. In particular, bed work, potential energy of added meltwater and dissipation cannot be omitted and subsequently described as exact total-energy conservation.

## 3. Cold ice supplies a spatially nonlocal, frequency-dependent response

An exactly soluble verification problem is a flat, stationary slab −D<z<0, with constant properties and a fixed lower temperature. Let its equilibrium temperature be T_0=T_m+G z, G>0, and let the interface perturbation be

    eta(x,t) = eta_hat exp(s t + i k x).

Linearizing the melting-temperature boundary condition at z=eta gives

    T'_i(0) = −G eta.

Write T'_i=−G eta_hat psi(z) exp(s t+i k x). The heat equation and lower boundary give

    psi'' − p² psi = 0,
    p² = k² + s/kappa_i,
    psi(−D)=0,  psi(0)=1.

Consequently,

    psi(z)=sinh[p(z+D)]/sinh(pD),
    H_D(s,k)=p coth(pD),
    delta q_i = −K_i G H_D(s,k) eta_hat.

The branch Re(p)>=0 is convenient in Re(s)>=0. H_D has units 1/m. The result is a standard diffusion Dirichlet-to-Neumann map, derived here for use as an exact code check; it is not claimed as an original mathematical discovery.

Its limits identify when a local constant-flux approximation fails:

* s=0, k→0: H_D→1/D, the derivative of the finite-slab conduction law with respect to thickness.
* |pD|≪1: H_D=1/D + D(k²+s/kappa_i)/3 + higher terms.
* Re(pD)≫1: H_D→sqrt(k²+s/kappa_i), independent of the lower boundary.
* k=0, s=i omega at high frequency: |H_D|≈sqrt(omega/kappa_i), phase≈45 degrees under the exp(i omega t) convention.

The overall flux perturbation includes the minus sign −K_i G eta; the quoted phase is for the positive admittance H_D, not for delta q_i relative to positive eta.

### A passivity check that prevents a false discovery

Multiply psi''=(k²+s/kappa_i)psi by its complex conjugate, integrate, and use the boundary values:

    H_D = integral[−D,0] |psi'|² dz
          + (k²+s/kappa_i) integral[−D,0] |psi|² dz.

Therefore Re(H_D)>0 for Re(s)>=0. In a conduction-only perturbation with fixed water heat input, Stefan balance gives

    rho_i L_f s eta_hat = −K_i G H_D eta_hat.

Taking real parts excludes a growing mode in Re(s)>=0. Cold conduction alone cannot be advertised as an instability mechanism in this problem. Any growing coupled mode must obtain feedback from the flow/heat-supply system. A numerical implementation violating this test is wrong.

This is a linear flat-slab result. It is NOT the thermal operator for a finite-amplitude stepped bed; that requires a mapped-domain diffusion solve, including the base temperature gradients and shape derivative of the normal flux. Do not insert the flat Fourier multiplier into a curved-base eigenproblem and call it exact.

## 4. A complete flat-base diagnostic dispersion system

This is a deliberately idealized verification model, not the proposed final environmental solver. Scale x and b by h_0, u by U_0, t by h_0/U_0, and theta by its positive base value theta_0. Let h=u=theta=1, and choose constant atmospheric heating to balance the base water-to-ice flux. Balance that flux with cold-ice conduction so the base bed is stationary. This balance is an explicit constructed equilibrium.

For u>0 choose q_w=alpha u_dim theta_dim, alpha constant, quadratic drag C_d, and optional specified horizontal momentum/thermal diffusivities nu and D_T. Define

    F=U_0/sqrt(g h_0),
    B=alpha/(rho_w c_p),
    A=h_a/(rho_w c_p U_0),
    epsilon=alpha theta_0/(rho_i L_f),
    a=s+i k,
    H(s,k)=h_0 H_D(s U_0/h_0,k/h_0).

Here alpha has units J/(m³ K), because q_w=alpha u_dim theta_dim. A,B,epsilon are dimensionless. All s,k in the matrix below are dimensionless. Normal-slope driving is selected so it balances base drag. Perturbations (h_hat,u_hat,theta_hat,b_hat) satisfy

    [ a                    i k                      0                    0       ]
    [ −C_d+i k/F²          a+2 C_d+nu k²            0                    i k/F²  ]
    [ 0                    B                        a+A+B+D_T k²         0       ] v = 0.
    [ 0                    epsilon                  epsilon              s+epsilon H ]

The depth term cancels from the temperature equation only because the chosen base heating and cooling balance. A non-equilibrium base changes that row. The conduction row follows Section 3. This determinant, with H depending on s, is a nonlinear eigenvalue condition. Replacing H by H(0,k) silently removes thermal memory.

Eliminating h_hat and theta_hat, away from their zero denominators, gives

    M(s,k)[s+epsilon H(s,k)]
      − (i k epsilon/F²)[1 − B/(a+A+B+D_T k²)] = 0,

    M=a+2 C_d+nu k² + (−C_d+i k/F²)(−i k/a).

Use the full determinant near a=0 or a+A+B+D_T k²=0, because elimination can omit roots or introduce poles. No roots of this coupled dispersion equation have been calculated in this package. This is a derivation to verify an implementation before proceeding to finite-amplitude profiles, not evidence of a new instability.

## 5. Finite-amplitude traveling profiles and an essential consistency trap

For the leading-order water equations, use a moving coordinate xi=x−c t, mean slope S_0, and

    b(x,t)=−S_0 x−E t+B(xi),
    h=H(xi), q=Q(xi), theta=Theta(xi),

with periodic B,H,Q,Theta of period ell. Then

    (Q−c H)'=0,
    −c Q' + (Q²/H+g H²/2)' = g H(S_0−B') − tau_b/rho_w,
    rho_w c_p [(Q−c H)Theta]' = q_a−q_w+Phi,
    rho_i L_f (E+c B') = q_w−q_i.

E is the mean vertical recession rate. The profile is a family indexed by ell and physical controls, with unknown c and E, jump conditions, mean-flow conditions and a phase condition fixing translation. At a critical point require regularity; do not integrate through a zero denominator with a numerical floor.

Two necessary integral checks are

    integral_0^ell(q_a−q_w+Phi) dxi = 0,
    rho_i L_f E ell = integral_0^ell(q_w−q_i) dxi.

A warm periodic stream cannot keep losing heat without an energy source. Resetting its temperature each period is an unreported heater. Likewise, with exact meltwater mass addition, periodic H,Q and everywhere positive m violate integrated water continuity unless a compensating source/sink is specified. Thus either use the declared asymptotic small-melt-volume model or solve an open finite reach. These are not optional bookkeeping details.

The finite-depth cold-ice base is also not automatically steady while the surface recedes. For a frozen-time calculation show that the base changes little during the perturbation growth time. Otherwise use transient tangent dynamics rather than a stationary spectrum.

For a steadily receding cold half-space an alternative is to solve the ice base in coordinates xi=x−c t and zeta=z+S_0 x+E t. Assuming stationary bulk ice and constant diffusivity gives

    −c T_xi + E T_zeta
      = kappa_i[T_xixi+2 S_0 T_xizeta+(1+S_0²)T_zetazeta].

The interface is zeta=B(xi), with T=T_m; a sufficiently distant lower boundary tends to a specified cold temperature. Verify depth independence. For S_0=0 and a flat interface with E>0, the exact solution is

    T(zeta)=T_infinity+(T_m−T_infinity)exp(E zeta/kappa_i), zeta<0.

Its conductive flux is rho_i c_i E(T_m−T_infinity). Stefan balance yields E=q_w/[rho_i(L_f+c_i(T_m−T_infinity))]. This gives a second exact thermal check. A nontrivial stationary, zero-recession, cold half-space with two different constant temperatures at its endpoints does not have the same steady solution. Do not use the receding solution at E=0 without taking the physical limit correctly.

## 6. Hydraulic jumps must move in the perturbation problem

For the hydraulic state U=(h,q), define F(U)=(q,q²/h+g h²/2). Across a jump at X(t), on a locally continuous bed and with bounded distributed source terms,

    [F(U)−Xdot U]=0,

where [v]=v_right−v_left. At a base jump moving at c, [F(U_0)−c U_0]=0. Perturb its position by zeta exp(s t). If v is the Eulerian state perturbation, the state sampled at the moving jump changes by v+zeta U'_0. Differentiation gives

    [(F_U(U_0)−c I)(v+zeta U'_0)] − s zeta[U_0]=0.

Dropping zeta fixes the jump artificially and can change the eigenvalues. Continuous bed elevation similarly requires

    [b_hat]+zeta[B'_0]=0.

For a uniform-width jump, let j=h(u−c)>0. Conjugate depths satisfy

    h_right/h_left = (sqrt(1+8 Fr_left,relative²)−1)/2.

The loss of specific head in its rest frame is

    Delta E=(h_right−h_left)^3/(4 h_right h_left).

The mechanical power loss per width is rho_w g j Delta E. In a zero-length, adiabatic jump idealization this becomes water internal energy, giving Delta theta=g Delta E/c_p. Neglect it only after checking g Delta E/(c_p theta_scale). The ordinary passive-temperature jump condition is the leading-order small-Eckert-number limit, not exact total energy conservation.

Enhanced bed heat exchange in the finite jump/recirculation zone is a separate effect. Do not deposit all jump power into a point of ice: that creates a singular recession rule and an arbitrary melt concentration. Resolve the finite zone or use an independently supported closure with uncertainty. For a finite concentrated energy source, thermal jump linearization must include the derivative of that source as well as the hydraulic terms.

## 7. Stability calculation and what it can establish

Linearize the coupled equations and interface conditions about each converged traveling profile, including the ice temperature field when retained. Schematically,

    s M(xi) v + [A_c(xi)v]' = S_U v + S_b b_hat' + thermal-shape terms.

Use piecewise smooth intervals separated at actual hydraulic jumps. Across successive periods impose a Bloch phase beta:

    v(xi+ell)=exp(i beta)v(xi),
    zeta_(n+1)=exp(i beta)zeta_n, beta in [−pi,pi].

Period-one perturbations alone miss step pairing and longer-wavelength coarsening. The translation mode at s=0,beta=0 is expected when the model and boundary constraints preserve translation symmetry. Separate it from unstable modes; test its residual. A numerically positive real part comparable to discretization uncertainty is unresolved, not an instability.

A shock-fitting collocation eigenproblem is the most direct first implementation. Evans–Lopatinski determinant machinery is an established alternative, not a mandatory badge of complexity. Compute a residual and convergence estimate for every reported eigenvalue. If using a nonlinear eigenvalue solver for eliminated ice diffusion, check against a formulation retaining the ice degrees of freedom.

Confirm representative modes in independently discretized nonlinear initial-value calculations: perturb the profile by an eigenmode, vary its amplitude, measure early growth/decay before nonlinear saturation, and compare the moving jump positions. Then examine multi-step interactions and initial-condition dependence. Spectral stability is local; it does not prove global attraction or predict a unique final spacing.

## 8. Physical and numerical relevance criteria

Record Froude U/sqrt(g h), Reynolds U h/nu_w (state this convention), aspect h/ell, water heat-advection/relaxation ratio, ice diffusion ratio D²/(kappa_i t_forcing), small-melt-volume ratio m ell/(U h), recession vs ice-advection time, and the jump Eckert ratio. Do not extrapolate a laminar closure to turbulent channels because a nondimensional Froude number happens to match.

The thermal-memory effect must change a physically meaningful output by more than reference discretization error and plausible closure spread. A shift in a spectral eigenvalue with no surviving change in melt distribution or geometry is insufficient for the proposed environmental claim.

Stop the cold-memory extension if the physically relevant limit is G≈0, response is effectively quasistatic, or conclusions are dominated by unknown water–ice heat transfer. Keep these null results. In temperate ice, the cold-gradient coupling derived here vanishes; it cannot support a universal glacier claim.

## 9. Checks actually completed

Run `python3 verify_foundations.py` to reproduce:

* Independent finite-difference solutions of the flat-slab diffusion boundary problem at five grid resolutions. Maximum final relative admittance error 1.745e−5; observed orders 1.994–1.999 for the nontrivial cases.
* 200 moving-hydraulic-jump directional derivative checks; maximum scaled error 4.49e−10.
* Jump head-loss identity checked on the same 200 cases; maximum relative error 4.67e−14.
* 1,000 passivity regression samples, in addition to the analytical proof above.

These checks verify selected algebra and subproblem discretizations. They do not establish existence, stability, novelty, field validity or superiority of the proposed model. The full dispersion system, curved-boundary operator and finite-amplitude coupled eigenproblem remain to be implemented and verified.
