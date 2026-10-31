# An algebraic lower bound on canonical seal memory

31 October 2026. This is a conditional result for the canonical infinite-pond, strictly subcritical linear model developed in the matching and stability notes. It is not an observed glacier-memory law, an exact large-time asymptotic, or an accuracy improvement on a finite forecast horizon. Historical priority is not certified.

## 1. Model and asymptotic ingredients

Use alpha=beta=1/2, q0=1, relative advection speed a0=1/2, and f=Fr^2 in [0,1/2). Write the shape generator as

\[
L_0=-a_0\partial_X-\lambda_j(X)I+
\lambda_j(X)\mathcal E,
\qquad \lambda_j=3V^2/D,
\]

where E is the positive spatial averaging operator and D=(1/2)V^(-3/2)-fV. Its positive semigroup is denoted T0(T). The projected shape input is

\[
g_s=-3V^2(\chi-V/5)<0.
\]

The nonnegative renewal memory kernel is

\[
K_0(T)=CT_0(T)g_s
=\int_{\mathbb R}w(X)
[T_0(T)(-g_s)](X)\,dX,
\]

\[
w(X)=D_Hm_0(X)k(X)>0,
\qquad \int_{\mathbb R}w(X)\,dX=D_H<\infty.
\]

Here C is a negative hydraulic functional. T denotes the fast inner time, not an unconverted number of physical seconds.

The base equation gives V'=2V^(3/2)[1+o(1)] upstream. Consequently

\[
V(X)\sim|X|^{-2},\quad
\lambda_j(X)\sim6|X|^{-7},\quad
-g_s(X)\sim\tfrac{12}{5}|X|^{-6}
\qquad(X\to-\infty).
\tag{1}
\]

The final relation uses chi/V->1. In particular, the jump rate is integrable over every upstream half-line.

## 2. Positive killed-transport contribution

Split the generator into killed transport and a positive jump term. The first term of its positive Dyson-Phillips expansion gives, for every nonnegative function F,

\[
[T_0(T)F](X)\ge
\exp\left[-\int_0^T\lambda_j(X-a_0s)\,ds\right]
F(X-a_0T).
\tag{2}
\]

No cancellation can remove this contribution because every term of the expansion is positive.

Fix any finite interval I=[X1,X2] of positive length. By (1),

\[
p_I=\exp\left[-a_0^{-1}
\int_{-\infty}^{X_2}\lambda_j(Y)\,dY\right]>0.
\]

For every X in I and every T>=0, the survival factor in (2) is at least p_I. Also W_I=integral_I w(X)dX is strictly positive.

There are c_g>0 and R>0 such that -g_s(Y)>=c_g|Y|^(-6) for Y<=-R. Choose T_I sufficiently large that X-a0 T<=-R throughout I and

\[
\tfrac12a_0T\le a_0T-X\le\tfrac32a_0T
\quad(X\in I,\ T\ge T_I).
\]

Applying (2) to -g_s and restricting the hydraulic observation integral to I therefore proves

\[
\boxed{K_0(T)\ge c_I T^{-6}\quad(T\ge T_I),}
\tag{3}
\]

where the explicitly positive constant can be chosen as

\[
c_I=W_Ip_Ic_g(3a_0/2)^{-6}>0.
\]

This is a lower bound. Contributions involving jumps may produce a larger or slower-decaying response; no matching upper asymptotic has been proved.

## 3. Consequences and exact limits of the conclusion

Equation (3) implies that no positive exponential moment exists:

\[
\int_0^\infty e^{\sigma T}K_0(T)\,dT=\infty
\quad\hbox{for every }\sigma>0.
\tag{4}
\]

In particular, K0 cannot be bounded for all sufficiently large T by A exp(-sigma T) for any finite A and sigma>0. This is compatible with the already proved finite kernel mass and strong decay in STABILITY_AND_SLOW_LIMIT.md. Strong decay does not imply exponential decay.

The lower bound does **not** prove that the first time moment diverges. Indeed, the first moment of T^(-6) is finite. It also does not prove finiteness of that moment, since a lower bound supplies no upper-tail control. Claims about an average memory time or a first-order slow-forcing error still require additional analysis.

For an exponentially stable finite-dimensional LTI approximation, its ordinary impulse response away from T=0 has a bound

\[
|K_f(T)|\le A(1+T^p)e^{-\sigma T}
\]

for some finite A,p and sigma>0. Equation (3) gives K_f(T)/K0(T)->0. Therefore

\[
\frac{|K_f(T)-K_0(T)|}{K_0(T)}\longrightarrow1.
\tag{5}
\]

Such a closure cannot retain uniformly small relative error over the entire infinite-time tail. This statement does not exclude good finite-horizon or absolute-error approximations. It does not apply to every finite-dimensional nonlinear model: nonlinear decay can be algebraic. It also does not establish that a neural network is necessary or superior.

The resummed head-to-discharge kernel contains D_H K0 as a nonnegative first term, so the corresponding algebraic lower-bound conclusion also applies to that resolved memory contribution. Its direct instantaneous term must be handled separately.

## 4. Physical and asymptotic interpretation

The mechanism in this reduced model is precise: some shape response is transported from increasingly far upstream without a jump. The upstream jump rate decays rapidly enough that this no-jump contribution has a survival probability bounded away from zero at any fixed observation interval. The input profile sampled there has an algebraic upstream tail, which generates (3).

This argument uses an infinite upstream pond and the inner model's instantaneous hydraulic matching. As T grows, it samples X approximately -a0 T. At fixed nonzero matching parameter nu, this eventually leaves the region in which the inner expansion represents a finite glacier lake. Finite pond length, outer lake structure, retained hydraulic storage or other omitted physics may cut off or alter the tail. The result therefore concerns the ordered idealized inner/infinite-pond limit, not an assertion that a real glacier retains an unbounded temporal memory law.

The model originates from the fast-scale inner-seal issue described in [Schoof et al. 2023, section 3.5 and Appendix B](https://doi.org/10.1017/jfm.2023.130). Positive-kernel and transport arguments are established techniques. The conditional prediction derived here is the non-exponential tail obstruction for this specified seal closure; originality and finite-environment significance remain separate questions.
