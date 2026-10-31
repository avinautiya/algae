# Canonical hydraulic matching and a linear state-space result

These propositions concern the canonical nonuniform inner seal in MATHEMATICS.md. They resolve a linear hydraulic matching step and define a well-posed linear state model. They do not establish nonlinear branch selection, transcritical flow, physical predictive superiority or historical originality. The motivation is the unsteady inner-region issue explicitly left open in [Schoof et al. (2023), section 3.5 and Appendix B.2](https://doi.org/10.1017/jfm.2023.130).

## 1. Canonical profile

Write r=Fr^2 with 0<=r<1/2, q0=1, alpha=beta=1/2, relative advection a=1/2 and upstream slope p_minus=1. Define

\[
F(v)=1-2v^3+v^{5/2},\quad
D(v)=\tfrac12v^{-3/2}-rv,
\quad 0<v<1.                                               \tag{1}
\]

**Proposition 1.** The equation V'=F(V)/D(V) has a unique solution up to translation with V(-infinity)=0 and V(+infinity)=1. It is strictly increasing. Choosing its translation so that V(0)=2^(-1/3), define B0(0)=0 and B0'=1-2V^3. Then B0 has a unique strict maximum at zero.

**Proof.** F(v)=(1-v^3)+v^(5/2)(1-sqrt(v))>0 and D(v)>0 on (0,1). The inverse coordinate X(v)=integral D(v)/F(v) dv diverges to minus infinity at zero and plus infinity at one. Thus it defines the asserted complete orbit. The bed derivative changes sign once, and B0''=-6V^2 V'<0. QED.

There is a translation-dependent constant c0 and a positive A_plus such that

\[
V(X)=(-X+c_0)^{-2}[1+O(|X|^{-5})],\quad X\to-\infty,
\]
\[
1-V(X)=A_+e^{-\kappa X}[1+o(1)],\quad X\to+\infty,
\qquad \kappa=\frac{7/2}{1/2-r}.                            \tag{2}
\]

For the first relation, Y=V^(-1/2) obeys Y'=-1+O(V^(5/2)); integration gives Y=-X+c0+O(|X|^-4). For the second, F'(1)=-7/2 and D(1)=1/2-r. In particular, the leading upstream coefficient is one, not an adjustable amplitude.

## 2. Hydraulic coefficients and the infinite-pond limit

On this base, use

\[
R=V^{5/2},\quad R_V=\tfrac52V^{3/2},\quad R_q=-\tfrac12V^{5/2},
\quad H_q=\tfrac12 V^{-1/2},
\]
\[
k=R_V/D,\quad l=R_q+kH_q,
\quad\Psi(x,y)=\exp[-\int_x^y k(z)\,dz].                   \tag{3}
\]

Both k and l are positive. Their endpoint behavior is

\[
k(X)\sim5|X|^{-6},\quad l(X)\sim2|X|^{-5}
\quad(X\to-\infty),
\]
\[
k_+=\frac{5/2}{1/2-r},\qquad
\frac{l_+}{k_+}=d_+=\frac25+\frac r5.                       \tag{4}
\]

Define the infinite-domain hydraulic head coefficient

\[
d(x)=\int_x^{\infty}\Psi(x,y)l(y)\,dy,
\qquad
d_- =\int_{-\infty}^{\infty}\Psi(-\infty,y)l(y)\,dy.        \tag{5}
\]

Here Psi(-infinity,y)=exp[-integral_-infinity^y k] is well defined because k is integrable at the upstream end. It decays exponentially at the downstream end.

**Proposition 2.** The quantities in (5) are finite and positive. The function d is bounded and continuous, with d(-infinity)=d_minus>0 and d(+infinity)=d_plus>0. Its finite upstream cutoff error is

\[
d(-L)=d_--\tfrac12L^{-4}+O(L^{-5}),\quad L\to\infty.         \tag{6}
\]

**Proof.** Upstream integrability follows from l~2|X|^-5 and bounded Psi. Downstream convergence follows from k->k_plus>0 and bounded l. Positivity is strict. The downstream limit is the constant-coefficient integral l_plus/k_plus. Let

\[
K_L=\int_{-\infty}^{-L}k= L^{-5}+O(L^{-6}),\quad
A_L=\int_{-\infty}^{-L}\Psi(-\infty,y)l(y)\,dy
 =\tfrac12L^{-4}+O(L^{-5}).
\]

Then d(-L)=exp(K_L)(d_minus-A_L), proving (6). Dominated convergence on the two tails and a compact middle interval establishes continuity and the upstream limit. QED.

At a finite downstream cutoff R, impose the normal-flow perturbation impedance

\[
\ell_R=V(R)/5,\quad
d_R^N=H_q(R)-D(R)\ell_R.
\]

The finite-interval coefficient is

\[
d_{L,R}=\Psi(-L,R)d_R^N+\int_{-L}^{R}\Psi(-L,y)l(y)\,dy.    \tag{7}
\]

The exact base reaches normal flow only asymptotically; (7) is a truncation, not an exact finite-distance normal-flow match. Since d_R^N-d(R)=O(exp(-kappa R)), its propagated error at a fixed upstream point is O(exp[-(kappa+k_plus)R]). A conservative O(exp(-k_plus R)) bound holds for arbitrary uniformly bounded terminal data. Constants can be chosen bounded as the left cutoff moves upstream. These estimates assume the exact base coefficients; errors from numerical base construction are additional.

## 3. A nonvanishing shape derivative after matching

For bounded continuous bed perturbations with finite endpoint limits, set

\[
(\mathsf E b)(x)=\int_x^{\infty}\Psi(x,y)k(y)b(y)\,dy.
\]

The weights integrate to one. Its upstream limit is

\[
E_-b=\int_{-\infty}^{\infty}\Psi(-\infty,y)k(y)b(y)\,dy.
\]

With b(-infinity)=0 and prescribed upstream total-head perturbation h, the matched linear flux is

\[
j=D_h h+Cb,\qquad D_h=1/d_-,\quad Cb=-E_-b/d_-.            \tag{8}
\]

Consequently a smooth nonnegative, nonzero compact bed perturbation, supported away from the crest, gives Cb<0. Scaling its amplitude sufficiently small preserves the crest height, position and curvature as well as both outer slopes. The shape sensitivity therefore persists in this canonical hydraulic infinite-pond limit.

For a fixed compact support, the unnormalized numerator's left-cutoff normalization error is O(L^-5). The normalized flux functional C has O(L^-4) left-cutoff error because d(-L) has the larger O(L^-4) error. The distinction matters when assessing numerical convergence. Downstream errors obey the corresponding terminal/cutoff estimates when the cutoff lies beyond the support.

This is a linear conditional identifiability statement. It does not establish that every admitted profile is reachable through environmental forcing, nor that a current-head rating is inaccurate for a particular dataset.

## 4. A full-line linear evolution

Let the real Banach space be

\[
\mathcal X=\{b\in BUC(\mathbb R): b(-\infty)=0,
\ b(+\infty)\text{ exists and is finite}\},
\]

with the sup norm. It is closed. Define

\[
\mathsf P b=(b-\mathsf E b)/D,\quad
\chi=(H_q-d)/D,\quad m=3V^2,
\]
\[
\mathcal L=-a\partial_X-m\mathsf P-m\chi C,
\qquad G=-m\chi D_h.                                     \tag{9}
\]

The auxiliary E need not map X into itself: its upstream limit generally is E_minus b rather than zero. However P does map X into X, since D diverges upstream, while E b-b tends to zero downstream. Also

\[
\chi\sim V\ (X\to-\infty),\qquad \chi\to1/5\ (X\to+\infty).
\]

All these operators are bounded where they should be: ||E||<=1, ||C||<=1/d_minus, D>=1/2-r, P is bounded, chi is bounded, and G belongs to X. The multiplier perturbation F=-mP-m chi C has upstream limit zero and downstream limit -3(Cb)/5, so F maps X into X. Uniform continuity follows from the smooth bounded coefficients and E'=k(E-b).

**Proposition 3.** The operator L, with the standard translation-generator domain, generates a strongly continuous semigroup on X. With gamma>0, the lake-coupled block operator

\[
\mathcal A=\begin{pmatrix}\mathcal L&G\\-C/\gamma&-D_h/\gamma\end{pmatrix}
\]

generates a strongly continuous semigroup on X times R. Hence the linear initial-value problem with continuous forcing i,

\[
b_T=\mathcal Lb+Gh,\qquad
\gamma h_T=i-j,\qquad j=D_hh+Cb,                          \tag{10}
\]

has a unique mild solution on every finite horizon.

**Proof.** The group T(T)b(X)=b(X-aT) preserves X and is strongly continuous by uniform continuity. Its generator is -a partial_X on its standard derivative domain. The remaining terms in L are bounded linear maps X->X as shown above. Apply the bounded-perturbation theorem. The lake coupling is another bounded block perturbation of the generator diag(-a partial_X,0). QED.

This is a linear well-posedness theorem, not a stability or discretization-convergence theorem. It does not supply a nonlinear wet/dry continuation or solve the moving-phase matching problem.

For externally specified h and suitable initial b, the resulting exact causal output is

\[
j(T)=D_hh(T)+C e^{\mathcal LT}b_{init}
 +\int_0^T C e^{\mathcal L(T-s)}G h(s)\,ds.                \tag{11}
\]

## 5. Genuine shape forcing and the remaining gauge question

The velocity gain chi is positive. To see this, let

\[
a_h=(R_V-D_VV')/D>0,\quad
b_h=(R_q-D_qV')/D<0.
\]

The inequalities follow from V'>0, D_V<0 and D_q>0 on this canonical base. The bounded normal-matched solution of chi'=a_h chi+b_h is

\[
\chi(x)=-\int_x^{\infty}
\exp[-\int_x^y a_h]b_h(y)\,dy>0.
\]

Moreover,

\[
\frac{G}{V^3}\to-3D_h\quad(X\to-\infty),\qquad
\frac{G}{V^3}\to-\frac35D_h\quad(X\to+\infty).              \tag{12}
\]

Thus G is not in span{1,B0'}=span{1,V^3}: it vanishes upstream, which would remove the constant term, but its ratio to V^3 is not constant. In the fixed upstream-datum model, forcing therefore excites geometry beyond a rigid translation and datum shift.

There is a further strictly model-scoped observation: C is a negative positive-weight functional and G is strictly negative, so CG>0. Strong continuity then makes the unprojected kernel C exp(LT)G nonzero for sufficiently small nonnegative T. This establishes memory in (10), not automatically a gauge-invariant physical memory law for the fully matched nonlinear seal. The uniform moving-gauge cancellation control in MATHEMATICS.md remains mandatory; projecting out moving-control conventions must be handled consistently rather than inferred from CG alone.

The remaining research questions are nonlinear branch/control selection, phase motion, reachability of perturbations, conservation and matching for the full hydraulic/thermal system, stability and growth rates, and whether the correction improves useful forecasts over the strongest existing coupled model. No experimental or historical-priority claim follows solely from these propositions.
