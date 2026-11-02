# Classical analytic Robin-slab reference

This is an independent implementation of the classical one-dimensional heat-equation eigenfunction solution. It is not a new mathematical theorem, proof of field accuracy or a historical novelty claim. Neither the finite-volume solver nor the one-node ODE is imported by `robin_reference.py`.

## Problem and assumptions

Let xi=x/d run from temperate ice at0 to the debris surface at1. Define total capacity C=rhoc*A*d, through-slab conductance K=k*A/d and diffusion time tau=C/K. Coefficients and forcing are constant; the initial slab temperature is zero everywhere. For a single water contact without a parallel bare patch:

    F = rho_water cp_water Q,
    H = F [1-exp(-Gd/F)],
    P = Qa + H Tin.

H=0 when F=0. The governing problem is

    T_t = (1/tau) T_xixi,
    T(0,t)=0,
    K T_xi(1,t)=P-H T(1,t),
    T(xi,0)=0.

P must be nonnegative for this nonnegative-temperature reference. The forcing step is not classically compatible with the initial derivative at the surface; the solution is understood for t>0 with the specified initial limit. There is no refreezing, moving material/ice domain, lateral diffusion, variable coefficients or water storage.

## Eigenfunction derivation

Set B=H/K and Ts=P/(K+H). The steady state is Ts*xi. Its deviation has zero temperature at0 and a homogeneous Robin condition at1. Separation of variables gives

    T(xi,t)=Ts*xi + sum cn sin(lambda_n*xi) exp(-lambda_n²*t/tau),
    lambda_n cos(lambda_n)+B sin(lambda_n)=0.

For B>0, each eigenvalue lies between (n+1/2)pi and (n+1)pi. For B=0 it equals the lower endpoint. Orthogonality of this self-adjoint boundary problem projects the zero initial temperature onto the modes:

    norm_n = 1/2 - sin(2 lambda_n)/(4 lambda_n),
    cn = -Ts [-cos(lambda_n)/lambda_n + sin(lambda_n)/lambda_n²] / norm_n.

The code solves the equivalent phase equation

    lambda_n - (n+1/2)pi - atan2(B,lambda_n)=0.

Its derivative is 1+B/(lambda_n²+B²)>0, so the bracket contains a unique root. This avoids misleading endpoint signs caused by sin(integer*pi) rounding. Sine/cosine at the root are evaluated using the Robin identity and hypot(lambda,B), including the large-B limit.

## Observable inventories

Positive heat into the ice is K*T_xi(0,t). Stored debris energy is C times the spatial mean temperature:

    ice_power = K [Ts + sum cn lambda_n exp(-lambda_n²*t/tau)],
    stored_energy = C [Ts/2 + sum cn(1-cos(lambda_n))/lambda_n exp(...)].

Let I_top be the time integral of T(1,t). Then boundary energy is P*t-H*I_top. In the physical water/air interpretation, atmospheric energy is Qa*t and water-to-debris energy is H*(Tin*t-I_top), retaining its sign. Integrating the heat equation over xi proves

    Qa*t + E_water_to_debris = E_ice + E_stored.

For a generic affine top boundary, only P*t-H*I_top is asserted; it is not silently relabelled as the entire water-stream budget.

## Stable integrated-series identities

Naively summing cn*(1-exp(...))/lambda_n leaves an algebraic truncation tail. That tail can dominate the physically tiny early ice heat even when256 and512 modes appear otherwise adequate. Integrate the transient deficit in time:

    v(xi) = integral_0^infinity [Ts*xi-T(xi,t)] dt.

It obeys v''=-tau*Ts*xi, v(0)=0 and v'(1)+B*v(1)=0. Direct integration yields

    a=(B+3)/(6(B+1)),
    v(xi)=tau*Ts*[a*xi-xi³/6].

Consequently the exact sum rules are

    sum cn/lambda_n = -Ts*a,
    sum cn*sin(lambda_n)/lambda_n² = -Ts/[3(B+1)].

The implemented integrated quantities therefore use only exponentially decaying truncated sums:

    E_ice = K [Ts*(t-tau*a) - tau*sum cn/lambda_n exp(...)],
    I_top = Ts*[t-tau/(3(B+1))]
            -tau*sum cn*sin(lambda_n)/lambda_n² exp(...).

This removes the algebraic tail without fitting or altering the physical parameters. It does not eliminate floating-point cancellation when the true early heat is almost zero. `roundoff_cancellation_scale_J` is an arithmetic conditioning diagnostic, not a certified interval-arithmetic bound or a confidence interval. Raw results are not clipped to manufacture positive heat. Do not quote a large relative error when the reference heat is at this scale; use an absolute/steady-flux normalization or higher precision.

## APIs

    solution = RobinSlab.from_contact(C, K, F, Gd, Tin, Qa, terms=512)
    result = solution.response(elapsed_seconds)
    temperature = solution.temperature(xi, elapsed_seconds)
    cell_mean = solution.cell_average(lo_xi, hi_xi, elapsed_seconds)

For parallel bare/debris exchange, the debris surface still has an affine boundary Qa+a*Tin-b*Ts, with a and b independently derived from that contact model. Use

    solution = RobinSlab.from_linear_boundary(C, K, Qa+a*Tin, b, terms=512).

Its returned water/atmospheric split is None; `boundary_heat_J`, `ice_heat_J`, `stored_energy_J` and `top_temperature_integral_K_s` remain meaningful. The caller may separately compute bare-ice heat and the whole-stream energy. All seven model physics cannot be inferred from a scalar Robin solution.

**Do not restart this zero-initial solution at every production step if the slab already has a nonzero profile.** Time-varying coefficients, geometry, reclassification and inherited profiles need their actual state or a separate spectral projection. The analytic reference is for a fixed, zero-initial benchmark, not a replacement that discards thermal memory.

## Independent checks

Ten tests cover eigenvalue brackets and equation residuals, initial/steady limits, direct Robin boundary residuals, finite-difference PDE residuals,256/512-mode stability, energy inventories, finite-volume space and time convergence to independent analytic cell averages, generic parallel-contact boundary agreement and the dry Neumann limit. Finite-volume refinement remains necessary: a small numerical energy ledger alone does not show that early-time heat has converged.

## Why a one-node time constant is not uniformly valid

The resolved slowest decay rate is lambda_0(B)²/tau. The single-node rate is (1+B)/tau when the same physical C and K are used. As B tends to0, these approach pi²/(4*tau) and1/tau respectively. As B tends to infinity, the resolved rate approaches pi²/tau while the single-node rate grows without bound. Thus a small water-side Biot ratio does not by itself make the debris uniform when the ice-side temperature is fixed. This is a classical spectral consequence; environmental consequences require a separately tested coupled model.
