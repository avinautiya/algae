# Thermal observability and regional ice-incision forecasts

Derived 31 October 2026. These statements concern an explicit fixed-flow thermal model. They are not Navier–Stokes existence results, a new heat-transfer law, or proof that a mechanism is new in the literature. The model-specific observation-placement statement below is a candidate contribution; priority remains under review.

## Physical problem

Let water with prescribed longitudinal velocity u(y) and depth H flow over a fixed temperate ice bed. Assume u is smooth, bounded, nonnegative, and positive for 0<y<=H. It may vanish at the bed, as a no-slip laminar profile does. Relative temperature theta obeys the transverse Graetz equation

    w(y) theta_x = k theta_yy, w=rho_w cp u(y), 0<y<H,
    theta(x,0)=0, theta_y(x,H)=0, theta(0,y)=f(y).

Here x is distance, not time. There is no axial conduction, buoyancy, thermal source, changing geometry or meltwater feedback. Define flow-weighted bulk temperature F(x)=integral w theta dy/integral w dy. For a region [a,b], 0<=a<b, the integrated bed-directed heat power is

    J[f] = C [F(a)-F(b)], C=W integral_0^H w(y)dy,

with channel width W. If all this heat melts temperate bed ice, holding this steady scenario for tau seconds gives mean lowering tau J/[rho_i Lf W(b-a)]. That is a latent-heat diagnostic under fixed hydraulics. It does not implement a moving interface or added-water feedback.

## Theorem 1 — exact placement criterion for bulk thermal information

Suppose only finitely many exact bulk temperatures F(x_1),...,F(x_m) are observed at distinct nonnegative distances. On the class of smooth initial profiles satisfying 0<=f<=T_cap and the bed/free-surface thermal boundary conditions, J is determined for **every** admissible f if and only if **both a and b are observation locations**.

If an endpoint is missing, there exist two smooth, nonnegative, bounded inlet profiles with identical observations at every x_j and different J. Their solutions remain nonnegative and become smooth downstream. The difference can be made arbitrarily small while preserving exact observational equivalence. This is an existence statement; it does not assert a large ambiguity for every dataset, nor on a restricted parametric profile family.

### Proof

First establish the spectral ingredients, including the no-slip endpoint. Let V={v in H^1(0,H): v(0)=0} with energy form integral k v' h', and use weighted state space L2(w dy). The weight is bounded and positive almost everywhere. Poincare's inequality makes the form coercive, and the embedding V into L2(w dy) is compact (it factors through the ordinary compact H1-to-L2 embedding). Its positive self-adjoint inverse is therefore compact. The eigenfunctions are complete, with positive eigenvalues tending to infinity, and satisfy

    -k phi_n'' = lambda_n w phi_n,
    phi_n(0)=0, phi_n'(H)=0.

Each eigenvalue is simple: two solutions with the same lambda and bed value zero are proportional by uniqueness of the second-order ODE, or by their constant Wronskian. The first eigenfunction phi_0 can be chosen strictly positive for 0<y<=H by the Rayleigh minimization/no-node property. A zero first derivative at the bed together with zero bed value would give the identically zero eigenfunction. In particular phi_n'(0) is nonzero for every n.

Integrating the eigen-equation gives the key observation fact:

    beta_n = integral w phi_n / integral w
           = k phi_n'(0)/(lambda_n integral w) != 0.   (0)

Thus **every thermal mode contributes to bulk temperature**, even with no-slip shear. For a finite expansion f=sum d_n phi_n,

    F(x)=sum d_n beta_n exp(-lambda_n x).

If J/C belonged to the span of the observation functionals, coefficients gamma_j would give, for every n,

    exp(-lambda_n a)-exp(-lambda_n b)
      -sum gamma_j exp(-lambda_n x_j)=0.                 (1)

Collect equal locations in this exponential polynomial. A nonzero real exponential sum with p distinct exponents has at most p-1 distinct real zeros. This follows by induction: divide by one nonvanishing exponential, differentiate to remove the constant term, then apply Rolle's theorem. Equation (1) has infinitely many distinct zeros lambda_n. Thus every coefficient at every distinct location must vanish. The coefficients at a and b can vanish only if both locations are observed. Sufficiency is the energy identity itself.

For an explicit finite-dimensional existence argument, retain as many modes as the number p of distinct locations among {a,b,x_j}. The corresponding exponential evaluation matrix is nonsingular by the same zero-count argument. If an endpoint is missing, there is a finite coefficient vector h annihilating all observations but not the target. Let h(y) be its sine expansion.

Positivity is not an obstacle. For a finite combination h, h/phi_0 extends continuously to the bed with limit h'(0)/phi_0'(0), and is bounded elsewhere because phi_0>0. Normalize max phi_0=1 and choose a finite B with |h|<=B phi_0. The profiles

    f_plus/minus = (T_cap/2) phi_0 +/- epsilon h,
    0<epsilon<=T_cap/(2B)

lie in [0,T_cap], are smooth, vanish at the bed and have zero derivative at the top. They give the same observations and different targets. The positive heat semigroup associated with this diffusion form preserves nonnegativity (equivalently, the maximum principle for positive interior weight). QED.

For constant plug speed, phi_n=sin((n+1/2)pi y/H), lambda_n=alpha[(n+1/2)pi/H]^2/U and beta_n=1/[H(n+1/2)pi/H]. The explicit inequality |sin((2n+1)t)|<=(2n+1)sin(t) on [0,pi/2] supplies B=sum(2n+1)|h_n|. `channel.modal_witness` implements this fully analytical positive continuum witness independently of the finite-volume operator.

### Consequence for a forward forecast

If all bulk stations are upstream of a future target reach, no finite collection of them makes this target exactly identifiable on the full admitted profile class. Additional stations can nevertheless make the remaining ambiguity very small. This theorem does not prohibit useful forecasts, identify a practical noise floor, or cover continuous-in-x observations. Exact continuous data would contain different information.

The theorem covers the prescribed no-slip/shear-free parabolic velocity used below as well as plug flow. It does not solve the momentum equation, admit reverse flow/recirculation, or cover velocity profiles vanishing on an interval. General axial conduction would also change the operator and proof.

## Theorem 2 — sharp discrete forecast interval and a checkable dual bound

A conservative transverse discretization gives

    M theta_x + K theta = 0,

where M=diag(rho_w cp u_i dy_i)>0 per unit width and K is symmetric positive definite with nonpositive off-diagonals. K includes bed conduction and an insulated top. M^(-1)K is similar to a positive symmetric matrix; S(x)=exp(-M^(-1)K x) is positive. Initial-temperature bounds and bounded source strengths define a box l<=z<=u. Measurement information is |E z-d|<=epsilon. Let the regional heat functional be c^T z and assume the feasible class is nonempty.

For every real measurement-weight vector lambda,

    c^T z <= lambda^T d + |lambda|^T epsilon
       + (c-E^T lambda)_+^T u
       + (c-E^T lambda)_-^T l.                         (2)

The minimum of (2) equals the sharp maximum over the feasible class. Apply it to -c for the sharp minimum. Both are attained because the class is compact. Their midpoint minimizes worst-case scalar absolute error; its exact conditional error radius is half the interval width.

### Proof

Write c^T z=lambda^T E z+(c-E^T lambda)^T z. The measurement band bounds the first term by lambda^T d+|lambda|^T epsilon; maximizing the second term over the box gives (2). Equality of primal and minimized dual follows from finite-dimensional linear-programming duality for a feasible bounded problem. Direct substitution of an attaining primal witness and a dual vector checks the result without trusting a reported optimizer objective. The midpoint argument follows by considering the two attainable scalar endpoints. QED.

The general LP duality and interval-minimax principles are established mathematics, not newly claimed general theorems. The code evaluates their physically derived operators, witnesses, feasibility and primal–dual gaps. Floating-point tolerances are disclosed; this is not interval-arithmetic certification.

## Theorem 3 — conservative inclusion of uncertain surface heating

For known spatial source intervals, let nonnegative surface heat fluxes q_j be bounded constants on [s_j,t_j]. The same fixed-flow equation becomes

    M theta_x+K theta=e_top q(x).

Then the state is affine/linear in the initial profile and q_j. For a region [a,b], conservation gives

    J = W [1^T M theta(a)-1^T M theta(b)
                +integral_a^b q(x)dx].                 (3)

Thus adding each q_j as a box variable preserves Theorem 2. The exact response to one unit flux is

    R_j(x)=integral_s_j^min(t_j,x) S(x-s) M^(-1)e_top ds,

with an empty integral equal to zero. Its target coefficient is

    W [1^T M(R_j(a)-R_j(b))
                 + length([a,b] intersect [s_j,t_j])].  (4)

This coefficient is nonnegative. Its value equals the integral of the positive bed-flux response over [a,b]. Surface heat supplied inside the reporting region must appear in (3); omitting it would infer the wrong melt power from downstream cooling.

### Proof

Variation of constants gives R_j. Multiplying the equation by 1^T cancels internal transverse face exchanges; 1^T K theta is the bed heat flux. Integrating in x gives (3) and (4). Positivity follows from the positive semigroup, nonnegative input and nonnegative bed-temperature output. Combining all linear responses gives the boxed-variable formulation of Theorem 2. QED.

These sources are prescribed net heating, not a radiation-transfer or weather model. The implementation restricts them to downstream of the forecast origin and nonnegative values. Cooling/freezing require a different admissible-state treatment.

## Corollary — the thermal moment that determines the target

With a known operator and no uncertain future source, write the target kernel at a current cross-section as c. For nonzero positive c, the weighted profile moment z_target=c^T theta/(sum c_i) determines this particular regional heat power exactly. An error bound |delta z_target|<=eta contributes at most (sum c_i)eta to its power error, before model/discretization error.

This is a target-specific adjoint observation, a known goal-oriented sensing principle. It generally requires spatially resolved thermal information and does not imply that one point thermometer supplies it. A different reach/horizon/operator has a different kernel. Its practical value is to specify what thermal information a regional erosion forecast needs.

## What the mathematical results establish

The complete continuum placement criterion and positive smooth counterexample are proved for the stated prescribed unidirectional velocity class, including no-slip shear. The constructive sharp forecast with uncertain heating is exact for the declared finite-dimensional class. It supports an executable prediction tool and identifies a falsifiable physical question: how much of the admitted ambiguity survives under justified inlet profiles, realistic mixing and resolved flow?

The computation does not show that turbulent glacier channels retain large hidden-profile effects, that the fixed-flow lowering is an environmental forecast, or that these are historically first results. PDE discretization, heat-source uncertainty and model discrepancy are distinct from the conditional information range and must remain separate.
