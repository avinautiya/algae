# Removing the outer translation mode from the linear seal law

This document derives an exact coordinate reduction of the canonical full-line linear model in MATCHING_THEOREM.md. It separates changing outer-control position from genuine inner-profile deformation. The result is a causal shape/discharge law consistent with the linearized outer shock speed. It does not establish nonlinear branch selection, environmental skill or historical originality.

Use q0=1, alpha=beta=1/2, relative advection a=1/2, p_minus=1, p_plus=-1 and ell=1/5. The base satisfies B0'=1-2V^3, V(-infinity)=0 and V(+infinity)=1. Its hydraulic maps P, C, chi and D_h are those already derived. Define

\[
m=3V^2,\quad W=V^3,\quad
\mathcal L_0=-a\partial_X-m\mathsf P.
\]

The fixed-reference equations are

\[
b_T=\mathcal L_0b-m\chi j,\quad
j=D_hh+Cb,\quad \gamma h_T=i-j.                            \tag{1}
\]

## 1. Exact translation covariance

Translate the base profile by a small displacement epsilon and add epsilon to its bed datum. The perturbed bed is B0(X-epsilon)+epsilon. Its upstream affine asymptote X+A_minus is unchanged, while the downstream affine intercept changes by 2 epsilon. Its first-order fields are

\[
b=\varepsilon(1-B_0')=2\varepsilon W,\quad
h=\varepsilon,\quad v=-\varepsilon V',\quad j=0.            \tag{2}
\]

The hydraulic relation and velocity map therefore imply

\[
CW=-D_h/2,\qquad \mathsf PW=-V'/2.                         \tag{3}
\]

Since W'=3V^2V', a=1/2 and m=3V^2,

\[
\mathcal L_0W=-\tfrac12W'-m\mathsf PW=0.                   \tag{4}
\]

These are exact full-line identities. A cropped numerical operator will satisfy them only to its matching/discretization accuracy; their defects must be measured, not silently set to zero.

## 2. Outer-control motion and the quotient state

Let r(T)=b(+infinity,T). The upstream bed perturbation remains zero. From (1), using P b->0, m->3 and chi->ell downstream,

\[
\dot r=-3\ell j=-\tfrac35j.                               \tag{5}
\]

Write

\[
b=s+rW,\qquad s(-\infty)=s(+\infty)=0,
\quad h_{rel}=h-r/2.                                      \tag{6}
\]

The intersecting outer affine beds have slopes +1 and -1. With their upstream intercept fixed, a change r in the downstream intercept moves their intersection horizontally and vertically by r/2. Thus h_rel is lake head relative to the **outer affine intersection**, not automatically relative to the actual smooth inner crest.

Equations (3)–(6) give the closed subsystem

\[
j=D_hh_{rel}+Cs,
\]
\[
s_T=\mathcal L_0s+g_sj,
\qquad g_s=-3V^2(\chi-V/5),
\]
\[
\dot h_{rel}=i/\gamma-(1/\gamma-3/10)j.                    \tag{7}
\]

Once j is known, recover r by (5), and then b and h by (6). No prescribed relaxation time, fitted delay or arbitrary imposed shock velocity appears.

The outer shock speed is recovered correctly. At fixed outer slopes, normal flow gives V_plus(q)=q^(1/5), and the Rankine–Hugoniot relative speed is a(q)=V_plus(q)^3/(p_minus-p_plus)=q^(3/5)/2. Therefore the physical control speed c(q)=U-a(q) has perturbation

\[
\delta c=-\tfrac3{10}j=\dot r/2.                           \tag{8}
\]

Equation (5) is consequently the linear outer shock-motion law, expressed through the downstream intercept state.

At the base crest X=0, W(0)=1/2. The actual inner crest-height perturbation is s(0)+r/2, so head above that crest is h_rel-s(0). When derivatives exist, the actual crest-position perturbation is r/2-s'(0)/B0''(0). The state s carries the distinction between the smooth inner maximum and the intersection of the outer asymptotes.

## 3. The shape-driving term is strictly negative

Let chi be the normal-matched velocity gain for changing q with the bed profile fixed. From the hydraulic linearization,

\[
\chi'=a_h\chi+b_h,\quad
a_h=(R_V-D_VV')/D,\quad b_h=(R_q-D_qV')/D.
\]

Set z=chi/V-1/5. Direct differentiation gives

\[
z'=a_rz-\tfrac25 V'/V,
\qquad
a_r=\frac{R_V-(D_V+D/V)V'}D>0.                            \tag{9}
\]

For this canonical base D_V+D/V=-V^(-5/2)/4-2Fr^2<0, while R_V,V',D are positive. The normal-flow condition gives z(+infinity)=0. Its bounded backward solution is

\[
z(x)=\tfrac25\int_x^\infty
\exp[-\int_x^y a_r]V'(y)/V(y)\,dy>0.                   \tag{10}
\]

Consequently chi>V/5 at every finite X and

\[
g_s<0,\qquad g_s(-\infty)=g_s(+\infty)=0.                 \tag{11}
\]

The source therefore drives genuine decaying shape deformation after removing the outer translation/datum mode. The uniform-flow cancellation benchmark remains a separate required control; its moving-gauge source is zero, as documented in MATHEMATICS.md.

## 4. Well-posed causal shape law and two kernel conventions

Let Y be the closed subspace of BUC(R) with both endpoint limits zero. The operator L0 has the form

\[
\mathcal L_0=-a\partial_X-\mu I+\mu\mathsf E,
\qquad\mu=m/D\ge0.
\]

Although E alone does not preserve the zero upstream limit, mu E does, because mu vanishes upstream. All the perturbation terms are bounded maps Y->Y. Transport with nonnegative killing generates a positive semigroup; adding the bounded positive operator mu E preserves positivity. Thus L0 generates a positive strongly continuous semigroup S0(T) on Y.

There are two distinct, compatible kernels.

First, retaining discharge j as the shape input, (7) gives the exact renewal relation

\[
j(T)=D_hh_{rel}(T)+CS_0(T)s_{init}
 +\int_0^T K_0(T-u)j(u)\,du,
\quad K_0(T)=CS_0(T)g_s.                                 \tag{12}
\]

Second, resolve the instantaneous discharge feedback into the shape generator:

\[
\mathcal A_s=\mathcal L_0+g_sC,\quad G_s=g_sD_h,
\qquad s_T=\mathcal A_ss+G_sh_{rel}.
\]

Its head-to-discharge representation is

\[
j(T)=D_hh_{rel}(T)+Ce^{\mathcal A_sT}s_{init}
 +\int_0^T K_1(T-u)h_{rel}(u)\,du,
\quad K_1(T)=Ce^{\mathcal A_sT}G_s.                        \tag{13}
\]

K0 is the renewal kernel for discharge driving shape; K1 is the already-resolved kernel for relative head driving discharge. Their gains, units and integrals must not be interchanged.

C is a negative positive-weight functional and g_s is strictly negative. Hence the rank-one map g_s C is positive. The bounded-perturbation theorem gives a positive semigroup generated by A_s. It follows that K0 and K1 are nonnegative. At zero they are strictly positive because Cg_s>0 and D_h>0; strong continuity makes them nonzero on an interval. In fact, the transport-with-killing lower bound preserves strict sign for these strictly signed source functions at any finite time.

This proves a nontrivial memory term for the specified translation-reduced linear seal model. The lake-coupled block operator on Y times R is

\[
\begin{pmatrix}
\mathcal A_s&G_s\\
-(1/\gamma-3/10)C&-(1/\gamma-3/10)D_h
\end{pmatrix},                                           \tag{14}
\]

a bounded perturbation of diag(-a partial_X,0). It therefore generates a strongly continuous semigroup and defines a unique linear forecast for compatible initial state and continuous forcing. Positivity of the hydraulic shape semigroups does not by itself prove stability of this lake-coupled block. The coefficient 1/gamma-3/10 can change sign.

## 5. Recovery of the steady hydraulic scaling

The exact steady family is

\[
V_q(X)=q^{1/5}V_0(X/q^{2/5}),\quad
B_q(X)=q^{2/5}B_0(X/q^{2/5}),\quad a(q)=q^{3/5}/2.          \tag{15}
\]

Both depth and kinetic hydraulic head scale as q^(2/5), and friction R is unchanged under this spatial scaling. Let A_minus=lim(B0-X), A_plus=lim(B0+X), and H_minus=lim(B0+H) upstream. Define the steady head relative to the outer intersection by

\[
h_{s0}=H_--(A_-+A_+)/2.
\]

It is positive without numerical quadrature:

\[
h_{s0}=H_++\int_{-\infty}^{\infty}(R-V^3)\,dX>0,
\qquad H_+=1+Fr^2/2,                                    \tag{16}
\]

since R=V^(5/2)>V^3. The integral converges by the endpoint asymptotics. Formula (16) follows by integrating the hydraulic head gradient -R and the two bed asymptotes.

Differentiate (15), remove the upstream datum shift and then its downstream offset mode. The resulting explicit branch derivative is

\[
s_{br}=\tfrac25[B_0-XB_0'-A_-]
 -\tfrac25(A_+-A_-)W,                                   \tag{17}
\]

which belongs to Y. Differentiating the steady equations gives

\[
\mathcal L_0s_{br}+g_s=0,
\qquad 1=D_h(2h_{s0}/5)+Cs_{br}.                         \tag{18}
\]

Thus the static branch slope is dq/dh_s=5/(2h_s0), recovering the q proportional to h_s^(5/2) law. This algebraic check does not require asserting a bounded resolvent (-L0)^(-1) at zero frequency; such an assertion can fail on an unbounded domain. The companion STABILITY_AND_SLOW_LIMIT.md supplies the separate shape-decay argument, integrable-kernel result and exact zero-frequency gain. It does not assert exponential decay of the full lake-coupled state.

The model now has a nontrivial, derived memory law, an outer-motion consistency check and the correct differentiated steady family. The remaining scientific questions include the nonlinear reduction, agreement with full unsteady hydraulics and heat transfer, control changes, environmental parameter ranges and actual forecast gains over the published model. None is settled by the linear closure alone.
