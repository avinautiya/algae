# Canonical matched linear stability and the slow-input limit

31 October 2026. This note concerns the canonical infinite hydraulic match and its gauge-projected linear shape dynamics. It does not establish nonlinear branch selection, physical forecast superiority, or exponential stability of the full infinite-domain bed state. The critical storage value recovered below is an existing drainage criterion in this canonical setting, not a new instability claim.

## 1. Canonical base and the correct reference geometry

Use alpha=beta=1/2, q0=1, relative advection speed a0=1/2, upstream slope +1 and downstream slope -1. The physical control speed is c=U-a; it must not be confused with the relative speed a. Write f=Fr^2, with 0<=f<1/2, and

\[
F(v)=1-2v^3+v^{5/2},\qquad
D(v)=\tfrac12v^{-3/2}-fv.
\]

The base velocity satisfies V'=F(V)/D(V), with V(-infinity)=0 and V(+infinity)=1. The bed satisfies B'=1-2V^3 and B''=-6V^2V'<0. The parameter f is the squared nondimensional Froude parameter in the equations; at the downstream base state the local squared Froude number is 2f. The strict subcritical restriction must not be replaced by f<1.

The reference used here puts the intersection of the two *outer affine bed asymptotes* at horizontal coordinate and elevation zero. Thus

\[
\bar B(\xi)-\xi\longrightarrow0\quad(\xi\to-\infty),
\qquad
\bar B(\xi)+\xi\longrightarrow0\quad(\xi\to+\infty).
\]

Bars are suppressed in most subsequent formulas, but this normalization is essential. It is not the gauge that sets the actual smooth inner crest elevation to zero.

If an implementation uses B(X) with asymptotes X+A_minus and -X+A_plus, define

\[
X_o=(A_{plus}-A_{minus})/2,\quad
Z_o=(A_{plus}+A_{minus})/2,
\]

\[
\xi=X-X_o,\quad
\bar B(\xi)=B(\xi+X_o)-Z_o,\quad
\bar h=H_{up}-Z_o.
\]

A crest-zero base therefore needs this conversion before using the formulas below. In the original coordinates the shape derivative is

\[
S(X)=\tfrac25\{B(X)-Z_o-(X-X_o)B'(X)\}.
\]

Using (2/5)(B-XB') without the intercept correction silently adds datum and translation modes.

## 2. Positive upstream head and the steady homothety

For the normalized base, let hbar be its upstream water-head correction. Changing variables from xi to v gives

\[
\bar h=
\lim_{v\downarrow0}\left[
v^{-1/2}-\int_v^1\frac{(1-s^3)D(s)}{F(s)}\,ds\right].
\]

The difference between the integrand and (1/2)s^(-3/2) is explicitly negative, so

\[
\bar h=1+\int_0^1
\frac{\tfrac12s(1-\sqrt{s})+fs(1-s^3)}{F(s)}\,ds>1.
\tag{1}
\]

The integrand is finite at both endpoints: F has a simple zero at one and the numerator vanishes there. This is a positive head above the outer-line intersection, not a head measured above the smooth inner crest.

There is an exact family of normalized steady profiles:

\[
V_q(\xi)=q^{1/5}V(\xi/q^{2/5}),\qquad
B_q(\xi)=q^{2/5}B(\xi/q^{2/5}),
\]

\[
a(q)=\tfrac12q^{3/5},\qquad h(q)=q^{2/5}\bar h.
\tag{2}
\]

Substitution verifies the hydraulic equation and travelling bed equation. Both outer slopes and both normalized affine intercepts remain fixed. This family does not hold the relative speed a, or the physical speed c=U-a, fixed while changing q.

Its bed-shape derivative at q=1 is

\[
S=\tfrac25(B-\xi B'),\qquad S(\pm\infty)=0.
\tag{3}
\]

It is strictly negative at every finite xi. Indeed, S'=-(2/5)xi B'', which is negative to the left of zero and positive to the right, while both endpoint limits are zero. Also S' tends to zero, so S belongs to the usual transport-generator domain on C0(R).

## 3. Hydraulic response, translation and shape projection

Let

\[
R=V^{5/2},\quad k=R_V/D,\quad
m_0(\xi)=\exp[-\int_{-\infty}^{\xi}k(y)\,dy].
\]

The matched hydraulic operator has the form

\[
j=D_Hh+Cb,\qquad
Cb=-D_H\int_{\mathbb R}m_0k b\,d\xi,
\quad D_H>0.
\tag{4}
\]

Here C is a bounded negative functional and integral m0 k=1. The velocity perturbation is v=Pb+chi j, with

\[
Pb=\frac{b-\mathcal Kb}{D},\qquad
\mathcal Kb(\xi)=\int_\xi^\infty
e^{-\int_\xi^y k}k(y)b(y)\,dy.
\tag{5}
\]

The kernel is positive and K applied to the constant one is one. The susceptibility obeys chi>0, chi/V->1 upstream and chi->1/5 downstream. More precisely w=chi/V-1/5 satisfies

\[
w'=a_rw-\tfrac25V'/V,\quad
a_r=\{R_V-(D_V+D/V)V'\}/D>0.
\]

Its bounded normal solution is a strictly positive backward integral. Consequently

\[
g_s=-3V^2(\chi-V/5)<0,
\qquad g_s\in C_0(\mathbb R).
\tag{6}
\]

For completeness, the fixed-frame linear bed equation is

\[
b_t=L_0b-3V^2\chi j,\qquad
L_0=-\tfrac12\partial_\xi-3V^2P.
\]

Its upstream intercept remains zero and its downstream intercept r satisfies r_dot=-(3/5)j. Put W=V^3. A translated base B(xi-delta x)+delta x has b=2 delta x W, h=delta x and v=-delta x V', so

\[
CW=-D_H/2,\quad PW=-V'/2,\quad L_0W=0.
\]

Defining

\[
s=b-rW,\quad h_{rel}=h-r/2
\]

therefore yields

\[
j=D_Hh_{rel}+Cs,\qquad
s_t=L_0s+g_sj.
\tag{7}
\]

Both endpoints of s vanish. This projection incorporates the linear Rankine-Hugoniot shift delta x=r/2 and delta x_dot=-(3/10)j. It does not hold the shock speed artificially fixed.

The actual smooth crest occurs where V^3=1/2. Its elevation perturbation is r/2+s(xi_crest). Thus the head above the actual crest is h_rel-s(xi_crest), not h_rel. Comparisons conditioned on the same actual crest head must account for that extra quantity.

## 4. Positive shape semigroup and integrable memory

On Y=C0(R), rewrite

\[
L_0=-a_0\partial_\xi-\lambda_jI+\lambda_j\mathcal K,
\qquad a_0=\tfrac12,\quad\lambda_j=3V^2/D\ge0.
\tag{8}
\]

The rate lambda_j is bounded and tends to zero upstream. K alone need not map C0 into C0 upstream, but lambda_j K does. The reaction/jump part is bounded, positive off the diagonal and satisfies the maximum principle because K averages future values with total weight one. Combining it with translation gives a positive contraction semigroup T(t) on Y. This construction uses the strictly subcritical base; no critical-point continuation is asserted.

Let

\[
K_m(t)=CT(t)g_s.
\]

Since C is negative and g_s<0, K_m(t)>=0. At zero, Cg_s>0, so the kernel is nonzero. The equations select its shape; no empirical relaxation time is introduced.

Differentiating the exact steady family (2), including a'(1)=3/10, gives

\[
L_0S+g_s=0,\qquad
1=\tfrac25D_H\bar h+CS.
\tag{9}
\]

As C and S are strictly negative, CS>0; equation (1) and D_H>0 imply

\[
0<CS=1-\tfrac25D_H\bar h<1.
\]

The semigroup identity now yields

\[
\int_0^tK_m(u)\,du=CS-CT(t)S\le CS.
\tag{10}
\]

Therefore K_m belongs to L1(0,infinity), and initially its mass M obeys the bound

\[
0<M\le CS<1.
\tag{11}
\]

This proves integrability without claiming a spectral gap. The next section closes the equality question. Neither argument proves a finite first time moment or exponential tail.

## 5. Decay of the steady shape derivative

Since L0 S=-g_s>0, positivity gives T(t)S increasing with t, and S<=T(t)S<=0. Let Z be its pointwise limit. It is bounded, measurable and dominated in absolute value by -S, which vanishes at both spatial ends.

To justify passing this limit through the semigroup, split the generator into killed transport -a0 partial_x-lambda_j and the positive jump operator lambda_j K. K is a proper probability kernel at every finite point. The Dyson-Phillips expansion is a countable sum of positive finite kernels; the jump rate is bounded, so the terms have factorially summable finite-time bounds. This defines the bounded-measurable extension and preserves pointwise bounded convergence. No particle or characteristic enters from spatial infinity in finite time in this construction; each jump kernel has no mass at infinity and there is no finite-time explosion of jump number.

It follows that T(tau)Z=Z for every tau>0. This step does not assume that Z is already in C0.

For any bounded measurable Z, KZ is continuous: its integrating-factor representation is a smooth factor times a locally absolutely continuous tail integral. The fixed-point killed-transport identity is

\[
Z(x)=e^{-\int_0^\tau\lambda_j(x-a_0u)du}Z(x-a_0\tau)
+\int_0^\tau e^{-\int_0^u\lambda_j(x-a_0v)dv}
\lambda_j(x-a_0u)(\mathcal KZ)(x-a_0u)\,du.
\]

Fixing a lower characteristic endpoint makes this a continuous, locally differentiable expression for Z. Hence Z is locally C1 and satisfies L0 Z=0. Its tail bound gives Z in C0. Set v_Z=PZ. Stationarity and the hydraulic identity give

\[
Z'=-6V^2v_Z,\qquad
(Dv_Z)'-R_Vv_Z=Z'.
\]

Consequently

\[
v_Z'=\frac{R_V-6V^2-D'}{D}v_Z
=\frac{V''}{V'}v_Z,
\]

where the last equality follows by differentiating D V'=F(V). Since V'>0 at every finite point, v_Z=c1 V' and Z=c1 B'+c2. Its two zero endpoint limits require c1=c2=0. Thus Z=0.

Dini convergence on every compact interval, together with the uniform tail domination |T(t)S|<=|S|, now gives

\[
\|T(t)S\|_\infty\longrightarrow0.
\tag{12}
\]

In fact T(t) is strongly stable on all C0. A compactly supported continuous function is dominated in absolute value by a constant multiple of -S, because -S is strictly positive at every finite point. Positivity proves its norm decay. Density of compactly supported functions and contraction extend this to every initial state in C0. This establishes strong stability of L0, not an exponential rate or a spectral gap, and does not establish stability of a different coupled generator.

## 6. Exact steady-family differential gain

Taking the limit in (10), using (12) and boundedness of C, proves the equality

\[
M=CS=1-\tfrac25D_H\bar h\quad\hbox{with}\quad0<M<1.
\tag{13}
\]

For zero initial shape, eliminating s gives

\[
j=D_Hh_{rel}+K_m*j.
\]

Since M<1, the corresponding head-to-discharge memory kernel can be resummed in L1:

\[
K_{head}=D_H\sum_{n=1}^{\infty}K_m^{*n},\qquad
\|K_{head}\|_1=\frac{D_HM}{1-M}<\infty.
\]

Every term is nonnegative. Thus j=D_H h_rel+K_head*h_rel has a direct response and a derived positive history contribution, without assuming exponential decay of the full shape semigroup. This is a bounded input-output statement for the head-forced closure, not for the complete lake/shape initial-value problem.

The zero-frequency head-to-discharge gain is therefore exactly the differential gain of the steady homothetic family:

\[
G_{DC}=\frac{D_H}{1-M}
=\frac{D_H}{1-CS}=\frac{5}{2\bar h}=G_{steady}.
\tag{14}
\]

The resolvent identity

\[
\widehat K_m(z)=CS-zC(zI-L_0)^{-1}S
\]

has its final term tend to zero as z decreases to zero along the positive real axis, by the proved bounded norm decay and the Abelian averaging formula. This limit is justified here; it must not be assumed for an arbitrary infinite-domain model.

For a bounded uniformly continuous head history H with consistent past, the resummed convolution also has a controlled slow-input limit. If the input is H(epsilon t), its difference from G_DC H(epsilon t) is bounded by

\[
\int_0^\infty K_{head}(u)
\omega_H(\epsilon u)\,du\longrightarrow0,
\]

where omega_H is a bounded modulus of continuity. Dominated convergence uses only K_head in L1. No O(epsilon) error rate follows without an appropriate first time moment. A calculation started with zero history has a separate startup transient.

These are continuum linear response statements. They do not prove that finite-domain truncation, mesh refinement and infinite time can be interchanged, or that a finite-domain fitted memory time converges uniformly. Nor do they recover the entire nonlinear steady flux relation from a linearization.

## 7. Observable lake-loop characteristic roots

With lake storage gamma and inflow perturbation i, the projected lake equation is

\[
\dot h_{rel}=i/\gamma-\mu j,\qquad
\mu=1/\gamma-3/10.
\]

The scalar characteristic function for the zero-initial-shape observable loop is

\[
\Delta(z)=z[1-\widehat K_m(z)]+\mu D_H.
\tag{15}
\]

For Re(z)>=0, equation (13) gives Re(1-Khat_m(z))>=1-M>0. If mu>0, a root would require z=-mu D_H/(1-Khat_m(z)), whose real part is strictly negative. Thus (15) has no unstable or imaginary-axis roots.

If mu<0, Delta(0)<0 and Delta(r) tends to positive infinity on the real axis. Moreover r[1-Khat_m(r)] is strictly increasing, so there is exactly one positive real root. At mu=0, a zero root occurs. This recovers the canonical critical value gamma=10/3; it does not introduce a new oscillatory instability.

These statements concern the observable Volterra characteristic equation. An infinite-domain transport generator can have essential spectrum or arbitrarily slow shape tails even when it is strongly stable. Absence of roots of (15) is not a proof of exponential decay of the entire lake-coupled shape state, nor of nonlinear stability or branch/basin selection. Those require separate analysis; no such extension to the coupled generator is asserted here.

The source model, inner matching caveat and existing drainage mechanism are in [Schoof et al. 2023, section 3.5 and Appendix B](https://doi.org/10.1017/jfm.2023.130). Semigroup, integrating-factor and positive-kernel arguments are established mathematics. The contribution investigated here is their derived closure for the specified unsteady inner seal, with the limitations above retained.
