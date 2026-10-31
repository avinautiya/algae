# Finite-bypass retention and delayed response

These are self-contained conditional propositions for PHYSICAL_MODEL.md. A proof establishes its conclusion under its assumptions; it does not establish historical originality. All equalities below are analytical. No simulations were used to select them.

## 1. Nonlinear retention identity

Let A > 0 be constant. Let F be continuous, locally Lipschitz and strictly increasing on positive head, F(0)=0. Let E be continuous and locally Lipschitz for positive discharge, with E(0)=0 and, for 0<w<=w0,

\[
0<\beta(w)=\frac{A E(F(w))}{F(w)}\le\beta_{\max}<1.
\]

Assume positive-head solutions are unique. Define

\[
R(w)=\frac{1}{1-\beta(w)},\quad
\Phi(w)=\int_0^w R(z)\,dz,\quad
J=H-\Phi(w).                                           \tag{1}
\]

The bounded integral is well defined even if beta does not have a limit at zero. Consider two trajectories starting from the same H0,b0, with J0=H0-Phi(w0)>0. The reference has u=0. The other has nonnegative integrable bypass u supported on [0,t_b], total C>0, and stays wet through that interval. Both use (P1).

**Proposition 1.** Both retain positive lake water. After the bypass ends, their wet heads tend to zero (or reach zero and stop). Their final levels satisfy the exact identity

\[
H_{\infty,u}-H_{\infty,0}
 =\frac1A\int_0^{t_b}\frac{\beta(w_u(t))}{1-\beta(w_u(t))}
 u(t)\,dt>0.                                           \tag{2}
\]

If beta_min <= beta(w_u(t)) <= beta_max < 1 during the pulse, then the extra retained volume obeys

\[
\frac{\beta_{min}}{1-\beta_{min}} C
 \le V_{\infty,u}-V_{\infty,0}
 \le\frac{\beta_{max}}{1-\beta_{max}} C.                 \tag{3}
\]

The lower bound is strict if beta_min>0. More generally (2) is strictly positive whenever incision is positive on a set with positive bypass volume.

**Proof.** Write Q=F(w) and D(w)=Q/A-E(Q)=Q(1-beta)/A>0. Then

\[
\dot w=-D(w)-u/A,
\quad \dot H=-Q/A-u/A,
\]

and differentiate (1):

\[
\dot J=-Q/A-u/A+R(w)(D(w)+u/A)
       =\frac{u}{A}\frac{\beta(w)}{1-\beta(w)}.          \tag{4}
\]

Thus J=J0 for the reference and J>=J0>0 for the bypass trajectory. Since Phi(w)>=w, both H=J+Phi(w)>0 and b=J+Phi(w)-w>=J0>0 while wet. After t_b, J is constant and w decreases. A positive limit of w would have D(w)>0 and hence cannot be a limit. Therefore w tends to zero, Phi(w) tends to zero, and H tends to J. Integrating (4) yields (2), and monotonicity of beta/(1-beta) yields (3). QED.

If u begins with positive flow, the bypass trajectory initially has the additional term -u/A in H_dot and hence initially lies below the reference. Equation (2) then implies at least one later crossing. For pulses that start continuously at zero, the same initial sign follows from any ordinary positive onset for which the direct removal term dominates its higher-order incision response. The final-level result itself does not require an onset regularity claim.

**Natural second-outlet corollary.** The higher fixed sill u=G((H-h_s)_+) described in PHYSICAL_MODEL.md turns off in finite time while w>0. Proposition 1 applies to its actual trajectory. The reversal therefore does not require imposing an instantaneous or arbitrarily large drainage pulse.

### Proportional incision

When E(Q)=kappa Q and beta=A kappa in (0,1), Phi(w)=w/(1-beta), and

\[
H_{\infty,0}=H_0-\frac{w_0}{1-\beta},\qquad
\Delta V_\infty=\frac{\beta}{1-\beta} C.                \tag{5}
\]

This holds for every admitted increasing rating F, not just a linear one. The ratio beta/(1-beta) is a retention amplification factor for this ideal model. It becomes singular as beta approaches one because the stable asymptotic state becomes sensitive there, not because real glaciers can produce an unlimited effect.

For a prescribed pulse, a convenient sufficient wetness condition is

\[
C+(1-\beta)F(w_0)t_b<A w_0.                             \tag{6}
\]

Indeed, during the pulse w<=w0 and therefore w_dot>=-[(1-beta)F(w0)+u]/A. Integrating gives w>0. Together with H0>w0/(1-beta), (6) gives a nonempty model class without an assumed future observation.

**Exact illustration.** Take beta=1/2, H0=3w0, b0=2w0 and C=A w0/2. Choose any finite pulse satisfying (6), for instance t_b <= A w0/[4F(w0)]. Then H_infty,0=w0 and H_infty,u=3w0/2. All sills stay above the fixed basin datum. These are algebraic levels of an ideal model, not a simulation, a field measurement, or an accuracy score.

## 2. Arbitrary hypsometry with proportional incision

Let V(H) be continuously differentiable with A(H)>0, and E(Q)=kappa Q, kappa>0. Here allow a common prescribed inflow I(t)>=0 of finite duration and total volume M_in, so A(H) H_dot=I-Q-u. Let both trajectories remain on a common wet, nonempty drainage branch and later end at w=0 after inflow and bypass have ceased. Define

\[
\Psi(H)=H-\kappa V(H).
\]

Suppose 1-kappa A(H)>0 on the interval between the two residual levels. Assume both final roots below exist on this branch; a root outside the physical basin is not admitted.

**Proposition 2.** For total finite bypass volume C,

\[
\Psi(H_{\infty,0})=b_0-\kappa V(H_0)-\kappa M_{in},\quad
\Psi(H_{\infty,u})=b_0-\kappa V(H_0)-\kappa M_{in}+\kappa C. \tag{7}
\]

and H_infty,u>H_infty,0. The residual level, considered as a function of cumulative bypass volume c, satisfies

\[
\frac{dH_\infty}{dc}=\frac{\kappa}{1-\kappa A(H_\infty)},
\qquad
\frac{dV_\infty}{dc}=\frac{\kappa A(H_\infty)}
 {1-\kappa A(H_\infty)}>0.                             \tag{8}
\]

If A_min <= A(H_infty(c)) <= A_max < 1/kappa, then

\[
\frac{\kappa A_{min}}{1-\kappa A_{min}}C
 \le\Delta V_\infty\le
\frac{\kappa A_{max}}{1-\kappa A_{max}}C.                \tag{9}
\]

**Proof.** Conservation and incision give

\[
\frac{d}{dt}[b-\kappa V(H)]=\kappa(u-I).
\]

At the final wet/dry boundary H=b, giving (7). Psi is strictly increasing on the stated branch, so its larger argument value requires a higher level. Implicit differentiation yields (8); integrating its volume derivative gives (9). QED.

Thus the effect of C does not require zero environmental water input: the same finite inflow enters both root equations and cancels from their difference. This statement requires a common constant kappa. An inflow temperature that changes kappa is not covered merely because its water-volume history is known.

At a branch fold 1-kappa A=0 the monotone-root assumption fails. Equation (8) is not valid through such a fold. Basin evolution, multiple roots or full lake emptying require a new analysis.

## 3. Exact finite-time crossing forecast

For constant A, E(Q)=kappa Q and F(w)=k w, define

\[
a=k/A,\quad\beta=A\kappa>0,\quad\lambda=a(\beta-1).
\]

Let u be any nonnegative finite pulse ending at t_b, with C=integral u. Both compared runs use the same initial H0,b0 and the same prescribed inflow I(t), which cancels from the perturbation equations; the reference has no bypass. Work only while the rating remains valid and both lakes remain wet and nonempty. If lambda is nonzero set

\[
M_\lambda=\int_0^{t_b}e^{-\lambda s}u(s)\,ds,\quad
\mu_\lambda=-\frac1\lambda\log(M_\lambda/C).
\]

Since u>=0, 0<=mu_lambda<=t_b. With I>=0, a sufficient condition for a positive persistent head throughout and after the pulse is M_lambda<A w0.

**Proposition 3.** For t>=t_b,

\[
\Delta H(t)=\frac{a}{A\lambda}
 [e^{\lambda t}M_\lambda-\beta C].                      \tag{10}
\]

If the crossing occurs after t_b and within model scope, its time is exactly

\[
t_\times=\frac{\log\beta}{a(\beta-1)}+\mu_\lambda.       \tag{11}
\]

For beta=1 these limits reduce to

\[
\Delta H(t)=\frac{C}{A}[a(t-\mu)-1],\quad
t_\times=1/a+\mu,\quad
\mu=\frac1C\int_0^{t_b}s u(s)\,ds.                     \tag{12}
\]

The crossing term log(beta)/[a(beta-1)] is positive for every beta>0 and tends to 1/a at beta=1. In the stable class 0<beta<1 with positive residual levels, the final difference from (10) is beta C/[A(1-beta)], consistently with Proposition 1. This is an open parameter class.

**Proof.** Head satisfies w_dot=lambda w-u/A. Variation of constants gives

\[
w_u-w_{ref}(t)=-\frac1A\int_0^t e^{\lambda(t-s)}u(s)\,ds.
\]

Substitute into Delta H_dot=-a Delta w-u/A and integrate once. Fubini's theorem gives (10), using a+lambda=a beta. Its only post-pulse zero is (11). Expanding M_lambda=C[1-lambda mu+o(lambda)] proves (12). QED.

An auxiliary outlet whose flow depends on H still generates an actual pulse u(t) to which these identities apply. Forecasting that pulse *before it occurs* requires its rating and initial state; inserting its future measured volume or centroid would be hindcasting. For prescribed intervention schedules, C and mu_lambda can be known at the forecast origin.

## 4. Local input-output zero

Use a differentiable general E(Q). At a constant-head travelling reference with E(Q*)=Q*/A, let k=F'(w*)>0, e=E_Q(Q*)>0, a=k/A, and lambda=k(e-1/A). H and b decline together; head and linearization coefficients are constant. Perturbations obey

\[
\begin{pmatrix}\delta\dot H\\\delta\dot b\end{pmatrix}
=\begin{pmatrix}-a&a\\-ek&ek\end{pmatrix}
\begin{pmatrix}\delta H\\\delta b\end{pmatrix}
-\begin{pmatrix}1/A\\0\end{pmatrix}\delta u.
\]

With zero initial perturbations, the transfer function is

\[
\frac{\delta H(s)}{\delta u(s)}
=-\frac{s-ek}{A s(s-\lambda)}.                          \tag{13}
\]

There is a right-half-plane zero at s=ek. For a varying reference this is only a frozen-coefficient statement; a global Laplace transfer function cannot be assigned to arbitrary time-varying coefficients. At fixed geometry e=0 the apparent cancellation gives -1/[A(s+a)], with no inverse response.

For variable A(H), the frozen Jacobian has first diagonal entry r-a, with r=(Q+u-I)A_H/A^2 for prescribed common inflow I. Its numerator is still -(s-ek)/A, while the characteristic polynomial is s^2-(r+lambda)s+r ek. This is local algebra, not a global stability theorem or proof for an arbitrary environmental hydrograph.

The same sign has a nonlinear tracking interpretation: if a bypass control keeps H exactly on a wet reference trajectory, an internal sill perturbation satisfies delta b_dot=E_Q F' delta b to first order. Maintaining that output can hide a growing difference in the control geometry. Such tracking uses signed perturbations about an admissible operating input and must satisfy physical control limits.

## What the proofs do not establish

The stable drainage criterion, phase-plane reduction, variation of constants, transfer-function zeros and comparison principles are established mathematics. The candidate research statement is the particular finite-bypass retention/crossing phenomenon in this ice-outlet setting and its prediction consequences. Whether that statement is historically new remains open. No new turbulent, free-surface or hydrofracture law has been proved.
