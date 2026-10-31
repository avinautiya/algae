# Forecast-focused mathematical foundation

These are applications of established linear-systems, projection and adjoint arguments. No new mathematical theorem or glacier discovery is asserted. Their purpose is to make a specific physical finding testable and to connect it to useful forecasts.

## 1. Why extend the previous response calculation

The preceding package evaluates C exp(A t) B for an initial perturbation of an autonomous base state. Environmental forcing changes over time, and a noncommuting sequence of operators cannot generally be replaced by an averaged A. Nor can a transient forcing be replaced by an initial perturbation without deriving that equivalence.

For consistently normalized perturbation states x and input amplitudes f, consider

```
x_dot = A_j x + B_j f_j,  t_j <= t < t_{j+1};  x(0)=0
y(T) = C x(T)
```

The base state may vary between intervals, but the same state basis and normalization must be used throughout. This is an approximation to a physical tangent problem; refine the time intervals and verify its response against the underlying equations.

Exact zero-order-hold propagation within an interval is

```
x_{j+1} = Phi_j x_j + Gamma_j f_j
Phi_j = exp(A_j dt_j)
Gamma_j = integral_0^{dt_j} exp(A_j s) B_j ds
```

The upper-right block of exp(dt_j [[A_j,B_j],[0,0]]) is Gamma_j. This handles singular A_j and neutral modes without A_j inverse. Chronological multiplication produces the full endpoint response H_raw from the stacked forcing amplitudes. Reordering intervals changes the result in general.

## 2. A forcing norm that survives temporal refinement

Use the physically declared time-average quadratic norm

```
||f||_T^2 = (1/T) sum_j dt_j ||f_j||_2^2
z_j = sqrt(dt_j/T) f_j,  f_j = sqrt(T/dt_j) z_j
```

Then ||z||_2=||f||_T. A Euclidean norm of raw sample amplitudes would change merely because the same history was divided into more intervals. The implemented response matrix H includes the sqrt(T/dt_j) column factors.

This norm is not automatically a thermodynamic energy. Each input channel needs a physical scale or uncertainty normalization. Temperature, discharge and geometric perturbations cannot be mixed by raw numerical magnitude. Different admissible norms can produce different optimal disturbances; disclose that sensitivity.

## 3. Matched-budget counterfactuals

Let R_raw f=0 express linear budget constraints. After normalization R=R_raw S, where S contains the sqrt(T/dt_j) factors. With orthogonal projector P onto null(R), admissible normalized inputs satisfy z=Pz.

For ||z||_2<=epsilon,

```
max ||y(T)||_2^2 = epsilon^2 sigma_max(H P)^2
```

The code computes a dominant right singular direction through the symmetric Gram matrix P^T H^T H P, then replays the corresponding forcing through the interval state recurrence. Budget residuals, replay disagreement and rank are returned, rather than hidden behind a scalar gain.

For water/heat inputs the constraint rows must be derived. For example, a linearized inlet sensible-power perturbation is

```
delta P_in = rho cp [Q_0 delta theta + theta_0 delta Q]
```

and its integral supplies a heat-budget row. A zero sum of temperature samples is not a matched-heat constraint when Q varies. Finite-amplitude pairs also contain delta Q delta theta and possibly property changes; check their exact water/enthalpy integrals in the nonlinear runs.

An ellipsoidal norm and equalities do not enforce temperature positivity, discharge bounds or realistic forcing smoothness. The returned signed direction must be scaled and checked against those conditions. The program does not solve a box-constrained or general nonlinear optimal-control problem.

## 4. Closed-form tests that rule out implementation artifacts

For x_dot=f over a horizon T, output x(T), the unrestricted time-average-L2 gain is T^2. A zero-integral forcing constraint makes the endpoint response zero. This result must remain invariant when time bins are subdivided.

For x_dot=-a x+f with a>0 and zero-integral forcing, the continuous-time squared gain is

```
G = T (1-exp(-2 a T))/(2 a) - [(1-exp(-a T))/a]^2
```

This follows by projecting the kernel exp[-a(T-s)] off the constant function under the time-average inner product. Piecewise-constant forcing approaches this value from below as the input subspace is refined. The tests check that behavior. These are known mathematical examples, not ice-melt findings.

## 5. Connect model errors to incision rather than only temperature norms

For a physical quantity J=c^T x(T), a small forcing/discrepancy r in

```
delta x_dot = A(t) delta x + B_q(t) r(t)
```

has linearized effect

```
delta J = lambda(0)^T delta x(0) + integral_0^T lambda(t)^T B_q(t) r(t) dt
-lambda_dot = A(t)^T lambda,  lambda(T)=c
```

The weight lambda^T B_q identifies which locations and times of a heat-flux error matter to the declared incision or control-lowering prediction. It is an adjoint sensitivity, not a newly discovered physical law. Its absolute-value integral gives a conditional error bound when bounds on r are supplied.

B_q must represent a *single conservative interface-exchange error*: extra heat into ice is removed from water using the same power and corresponding phase-change mass. Perturbing the melt equation alone would inject energy and produce an unphysical sensitivity. Preserve geometric area and water/ice density conversions.

The current executable returns interval endpoint response matrices and optimal inputs; it does not yet derive or spatially discretize the physical B_q/adjoint. That requires the verified coupled model and a declared observable.

## 6. Scope of a nonlinear statement

If the physical perturbation equation has a remainder R with a justified quadratic bound in a neighborhood, variation of constants gives an O(epsilon^2) correction to the linear prediction over a finite horizon. The constants depend on the propagator, derivatives, allowed inputs and neighborhood.

They are not known just because several finite differences agree. Wet/dry changes, jump displacement, remeshing or control switching may destroy differentiability in the chosen coordinates. Use an appropriate interface coordinate, event sensitivity or a direct finite-amplitude comparison. Do not present the computed linear gain as a certified bound on arbitrary glacier behavior.

## 7. The physical claim to earn

Demonstrate that a specifically identified free-surface thermal–hydraulic mechanism makes a forecast-relevant response differ from the strongest suitable equilibrium/local model, beyond numerical and closure uncertainty. Verify the predicted response in independent nonlinear runs and on held-out configurations. A large gain in a manufactured matrix is only a software test.

Earlier event-averaging and reduced-network studies already produced limited or negative novelty results. This work does not erase them or recycle the statement that forcing order can matter. A new claim must be tied to the newly resolved physics and a measurable predictive consequence.
