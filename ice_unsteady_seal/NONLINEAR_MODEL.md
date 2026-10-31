# Nonlinear finite-interval model and local solution criterion

The executable nonlinear model retains the entire bed profile and solves the instantaneous spatial hydraulic boundary problem at every bed/lake time step. It implements a bounded extension of the linear prototype, not a completed moving-crest or outer shock-selection algorithm.

## Equations and boundaries

In the fixed reference frame of MATHEMATICS.md, write the actual bed as B0+b and actual upstream total head as H0+h. The nonlinear perturbation equations are

\[
b_T+a b_X=-[V^3-V_0^3],\qquad
\gamma h_T=i-[q-q_0].                                    \tag{1}
\]

Hydraulics solve

\[
(B_0+b+\mathcal H(V,q))_X=-R(V,q),
\]
\[
\mathcal H=(q/V)^\beta+\tfrac12Fr^2V^2,\quad
R=q^{\alpha-1}V^{3-\alpha},\quad
D=-\mathcal H_V>0.                                      \tag{2}
\]

The finite upstream bed perturbation is b(X_left,T)=0. Upstream total head is H0+h. At the right cutoff, R(V_right,q) is held at its reference value, imposing the same normal-flow terminal relation used in the linearization. The canonical base reaches true normal flow only at infinity, so this is an explicit cutoff approximation. Positive liquid velocity and discharge and a declared criticality margin are required.

There is no newly fitted relaxation time, artificial branch switch or neural correction. The time dependence resides in the evolving bed field. No finite-speed water wave is represented: the hydraulic problem remains quasisteady under the stated small-parameter ordering.

## A conditional local existence and uniqueness argument

Consider a reference solution on a fixed finite interval with strictly positive velocity and discharge, a uniform positive lower bound on D, and finite smooth hydraulic coefficients. Stay in a sufficiently small neighborhood preserving these margins. For given q and bed B, the subcritical inverse of H(V,q) expresses velocity smoothly as a function of hydraulic head minus B. Thus the hydraulic equation can be written as a first-order head ODE whose coefficients depend continuously and locally Lipschitzly on B in the sup norm. Its right boundary is smooth in q and B_right.

The derivative of the upstream head with respect to q, keeping B fixed and using the normal terminal relation, is the positive coefficient d(X_left) in MATHEMATICS.md. Its integrand is positive because

\[
R_q+(R_V/D)\mathcal H_q\ge2R/q>0,
\qquad \mathcal H_q-D\ell>0.
\]

The implicit-function theorem therefore supplies a locally unique smooth discharge q(B,H) on this connected admissible branch. The velocity profile inherits local Lipschitz dependence. This does not assert that a solution exists for arbitrary head or across dry/critical states.

Substituting q(B,H) and V(B,H) into (1) gives a transport equation with a locally Lipschitz nonlocal source coupled to a scalar ODE. For compatible continuous initial/inflow-boundary data, its characteristic integral formulation is a contraction on a sufficiently short interval. It supplies a local mild solution, continued while the hydraulic margins and bounded-state neighborhood persist. This uses established ODE/implicit-function/transport theory; it is not a new general PDE existence theorem.

This argument establishes local mathematical coherence of the declared finite problem. It neither proves global continuation nor supplies the full nonlinear infinite-pond/moving-control match.

## Numerical implementation

The nonlinear BVP uses log velocity and log discharge to represent positive variables, a cubic interpolation of the bed perturbation and an adaptive collocation mesh. A loss of hydraulic admissibility, failed collocation solve or exhausted mesh budget stops evaluation. Values are not clipped back into the admissible branch.

Bed advection uses first-order upwinding of the perturbation. Subtracting the exactly integrated base equation preserves the declared base rather than introducing a coarse-grid forcing from its steep spatial gradient. Time stepping is RK4 with an advection CFL cap; this is an implementation choice, not an unconditional nonlinear stability proof. Inflow and outflow integrals use the same RK stages.

The finite-cutoff boundary variable is total head. The reported exact bookkeeping concerns the surrogate gamma*h state. The channel surface-stage change differs by the kinetic-head correction, which is serialized separately. Hence a small surrogate residual cannot be advertised as physical lake-water closure. The infinite upstream hydraulic limit identifies where total head and stage coincide; the remaining dynamic matching and discretization errors still require evaluation.

The linear and nonlinear codes use the same physical closure and related boundary assumptions. Agreement checks linearization and implementation. It is not independent physical validation or a demonstration of superiority over an existing coupled glacier model.
