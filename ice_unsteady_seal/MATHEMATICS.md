# A causal inner-seal model and its unresolved matching problem

This formulation targets the specific unsteady inner-region issue acknowledged in Schoof et al. (2023), section 3.5, equations 3.17–3.18 and Appendix B.2: the steady discharge map is used on a time scale on which the inner bed profile also evolves. It is a candidate extension of that model, not a completed proof of the full lake-to-channel matched problem. The mathematics below uses the paper's dissipation-driven incision closure; it does not introduce an independently validated heat-transfer law.

## 1. Ordering and conservation

Use the original dimensionless equations, with P(S)=S^alpha, h(S)=S^beta, 0<=alpha<1 and beta>0:

\[
\delta S_t+(uS)_x=\varepsilon u^3P(S),
\]
\[
\nu Fr^2 S(\delta u_t+uu_x)=-u^2P(S)-Sb_x-\nu S h(S)_x,
\qquad b_t+Ub_x=w-u^3.                                      \tag{1}
\]

For an inner region, set X=(x-x_c(t))/nu, T=t/nu, S=Sigma, u=V, b=z(t)+nu B(X,T), and c=dx_c/dt. The mass equation becomes

\[
\delta(\Sigma_T-c\Sigma_X)+(V\Sigma)_X
 =\nu\varepsilon V^3P(\Sigma).                              \tag{2}
\]

The momentum equation is

\[
Fr^2\Sigma[\delta V_T+(V-\delta c)V_X]
=-V^2P(\Sigma)-\Sigma B_X-\Sigma h(\Sigma)_X.                \tag{3}
\]

Thus delta<<1 and nu*epsilon<<1 permit spatially constant q(T)=V Sigma and quasi-steady water hydraulics while retaining B_T at order one. This is a distinguished time-scale reduction, not permission to discard storage arbitrarily. It requires bounded inner gradients and a subcritical margin. Near criticality, or for forcing so fast that delta*frequency is not small, retain (2)–(3). No transcritical or hydraulic-jump theorem follows here.

For comparison, in physical rectangular-channel variables with bed retreat E, meltwater mass source rE and r=rho_i/rho_w,

\[
h_t+q_x=rE,\quad b_t=-E,\quad
\eta_t+q_x=(r-1)E,\quad\eta=h+b.                            \tag{4}
\]

Retaining meltwater injection while dropping the associated moving-bed storage can create fictitious choking. Equation (4) is a mandatory conservation audit. The present leading-order inner model drops both small terms consistently under (2), rather than claiming a new injection-induced instability.

## 2. Nonlinear hydraulic and bed state

Define hydraulic head and friction loss

\[
\mathcal H(V,q)=(q/V)^\beta+\tfrac12 Fr^2V^2,
\qquad R(V,q)=q^{\alpha-1}V^{3-\alpha}.
\]

Equation (3) reduces to

\[
(B+\mathcal H(V,q))_X=-R(V,q),\qquad
D=-\mathcal H_V=\beta q^\beta V^{-1-\beta}-Fr^2V>0.          \tag{5}
\]

The bed equation is

\[
B_T+(U-c)B_X=w-\dot z-V^3.                                \tag{6}
\]

In the upstream moving gauge used in Appendix B, the right-hand side is (U-c)p_minus-V^3+V_minus^3. Omitting B_T recovers the published steady inner profile. Retaining it makes q depend on an initial profile and its subsequent evolution, not merely on a current head correction and two slopes.

For a fixed reference frame and fixed reference datum velocity, a stationary base satisfies

\[
aB_0'=a p_--V_0^3,
\quad B_0'+\mathcal H_V(V_0,q_0)V_0'+R(V_0,q_0)=0,
\quad a=U-c_0>0.                                         \tag{7}
\]

A concrete nonuniform seal has q0=1, alpha=beta=1/2, a=1/2, p_minus=1, V_plus=1 and Fr^2<1/2. Then

\[
B_0'=1-2V_0^3,\qquad
V_0'=\frac{1-2V_0^3+V_0^{5/2}}
 {\tfrac12V_0^{-3/2}-Fr^2V_0}.                             \tag{8}
\]

Its velocity increases from the pond limit zero toward one; bed slope changes from +1 to -1. Anchor B0=0 where V0=2^(-1/3), the crest. A finite interval with V_left>0 is a pond truncation. It is not automatically the infinite-pond matched limit.

## 3. Exact finite-interval hydraulic response

In sections 3–5 only, write the finite matching interval as [0,L] by shifting its upstream cutoff to zero. The crest need not be at this origin. The canonical base in section 2 and moving-crest gauge in section 6 instead place the crest at X=0. The executable keeps its explicit grid with upstream X_left<0 and crest X=0; its upstream evaluations correspond to the lower endpoint of the formulas below.

Let b,e,v,j denote perturbations of inner bed, total hydraulic head B+H, velocity and spatially uniform discharge. Coefficients below are evaluated on the nonuniform base. Linearization gives

\[
e=b-Dv+\mathcal H_qj,\qquad e_X=-R_Vv-R_qj.                \tag{9}
\]

At a downstream normal-flow match impose v(L)=ell*j, with

\[
\ell=-R_q/R_V=(1-\alpha)V_+/[(3-\alpha)q_0].               \tag{10}
\]

The canonical heteroclinic reaches V_plus only asymptotically. At a finite cutoff its normal-flow residual must be reported. Imposing (10) there is a specified terminal perturbation impedance, approximating the far-field match; it is not an exact finite-distance match to the infinite profile. Downstream cutoff convergence is consequently required.

Set k=R_V/D>0, l=R_q+k H_q and d_L=H_q(L)-D(L)ell. Then

\[
e_X=k(e-b)-lj,\quad e(L)=b(L)+d_Lj.
\]

With Psi(x,y)=exp[-integral_x^y k(z) dz], define

\[
(\mathsf Eb)(x)=\Psi(x,L)b(L)
 +\int_x^L\Psi(x,y)k(y)b(y)\,dy,
\]
\[
d(x)=\Psi(x,L)d_L+\int_x^L\Psi(x,y)l(y)\,dy.               \tag{11}
\]

The solution is exactly e=E b+d j. With upstream total-head perturbation e(0)=h,

\[
j=D_hh+Cb,\quad D_h=1/d(0),\quad Cb=-(\mathsf Eb)(0)/d(0),
\]
\[
v=\mathsf Pb+\chi j,\quad
\mathsf P=(I-\mathsf E)/D,\quad \chi=(\mathcal H_q-d)/D.   \tag{12}
\]

There is no fitted relaxation time. For this constitutive class, d is positive: H_q/D>=V/q; therefore l>=2R/q>0, and d_L>0 because ell*q/V=(1-alpha)/(3-alpha)<1. E is a positive averaging operator: Psi(x,L)+integral_x^L Psi k=1. These provide useful algebraic verification checks.

**Shape-sensitivity proposition.** On a finite uniformly subcritical interval, let b be smooth, nonnegative, nonzero and compactly supported away from both boundaries and the crest. Keep upstream head and downstream normal-flow slope fixed. Then

\[
j=-\frac1{d(0)}\int_0^L\Psi(0,y)k(y)b(y)\,dy<0.            \tag{13}
\]

This follows immediately from (11)–(12). A sufficiently small perturbation preserves the crest height, location and curvature, outer slopes and wetness. Thus those instantaneous descriptors do not determine q within this finite-interval model. This is a conditional shape derivative, not a proof that all admitted profiles are reachable from realistic forcing histories.

## 4. Causal state and forecast operator

In the fixed reference frame, with the reference datum velocity held fixed, the prototype bed and lake perturbations are

\[
b_T+a b_X=-m(X)v,\quad m=3V_0^2,
\qquad \gamma h_T=i-j.                                   \tag{14}
\]

Here b(0,T)=0 clamps the finite upstream bed perturbation. It is a stated finite-pond boundary condition, not a derived infinite-pond match. The lake input i is a perturbation of supplied inflow; the current h is its evolving state.

Using (12),

\[
b_T=\mathcal Lb+Gh,\quad
\mathcal L=-a\partial_X-m\mathsf P-m\chi C,
\quad G=-m\chi D_h,
\]
\[
\frac{d}{dT}\binom b h
=\begin{pmatrix}\mathcal L&G\\-C/\gamma&-D_h/\gamma\end{pmatrix}
\binom b h+\binom0{i/\gamma}.                             \tag{15}
\]

On a fixed spatial discretization this is a fully specified finite-dimensional linear forecast system. Its exact state evolution is a matrix exponential/convolution. With externally prescribed h, its discharge is

\[
j(T)=D_hh(T)+C e^{\mathcal LT}b_{init}
 +\int_0^T C e^{\mathcal L(T-s)}G h(s)\,ds,                \tag{16}
\]

and its zero-initial-state admittance is D_h+C(sI-L)^(-1)G. Initial-profile memory must be retained separately. Stability, a nonzero forced-memory kernel, and a positive performance advantage are not assumed from this formula.

For a continuum finite-interval construction with bounded smooth coefficients, the transport form admits a characteristic integral equation. The hydraulic maps in (11) are bounded in the sup norm, including the downstream trace. Substitution into (14) gives a contraction on a sufficiently short time interval; continuation and Gronwall estimates yield a unique bounded mild solution on each finite horizon for compatible continuous data. This argument does not prove stability, infinite-domain convergence, full nonlinear well-posedness, or a semigroup approximation rate for a particular discretization.

## 5. Required gauge cancellation control

For a uniform normal-flow base, R,V,D,H_q,k are constant and R_q+R_V ell=0. Hence l=k*d_L, d(x)=d_L, and

\[
v=\mathsf Pb+\ell j.
\]

In the upstream moving gauge the bed perturbation is

\[
b_T+a b_X=-3V_0^2(v-\ell j)=-3V_0^2\mathsf Pb.             \tag{17}
\]

Thus head forcing does not excite bed shape from b_init=0: the correct result is j=h/d_L. A claimed history kernel in this matched uniform control can be a gauge artifact. The fixed-datum equation (14) intentionally differs; a prototype must test (17) separately and label the distinction. A nonuniform seal requires direct analysis of the projected forcing and its observable kernel rather than assuming memory survives this audit.

## 6. Smooth-crest moving gauge and remaining closure

For a smooth strict crest, use B(0,T)=B_X(0,T)=0 and B_XX(0,T)<0. Equations (5)–(6) imply

\[
\dot z=w-V_c^3,\qquad
c=U+3V_c^2V_{X,c}/B_{XX,c},
\]
\[
B_T+(U-c)B_X=V_c^3-V^3,
\qquad\gamma h_T=Q-q-\gamma w+\gamma V_c^3.                \tag{18}
\]

For an upstream pond B~p_minus X+A_minus(T), V->0, its intercept must evolve as

\[
(A_-)_T=V_c^3-(U-c)p_- .                                 \tag{19}
\]

Arbitrarily fixing that intercept can manufacture a memory effect. Equation (18) is a conditional smooth-crest PDE–DAE formulation; it fails at vanishing curvature, crest splitting or non-smooth control changes. A full matched solution must determine upstream affine matching, downstream normal-flow matching, lake storage counted exactly once, the phase condition and shock-motion selection consistently. The limit V_left->0 can be singular because D diverges. The companion [MATCHING_THEOREM.md](MATCHING_THEOREM.md) proves the canonical hydraulic normalization limit, persistent shape sensitivity, and well-posedness of a specified full-line linear state model. It does not prove the moving-crest nonlinear match or convergence of the executable discretization.

At a finite positive-velocity cutoff, the boundary variable e is total hydraulic head. The executable reports gamma*h bookkeeping as a **surrogate total-head ledger**, and separately reports the channel surface-stage correction h-Fr^2*V_left*v_left. Roundoff in the surrogate is not a proof of physical lake-water closure. Total head and surface stage coincide in the matched upstream pond limit, but finite-cutoff kinetic and missing upstream friction corrections must be assessed explicitly.

The possible research contribution is therefore a derived causal seal law and a demonstrated, consequential correction to instantaneous discharge selection. The present finite-pond prototype is a vehicle for checking that hypothesis. It is not yet a resolution of the published asymptotic gap or a new environmental discovery.
