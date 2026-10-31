# Independent mathematical and physical rejection tests

31 October 2026. This review distinguishes a useful linear prototype from a completed glacier-drainage model. The current prototype concerns a finite, positive-velocity, strictly subcritical interval, with a declared upstream bed perturbation and downstream normal-flow matching. It is not the full infinite-domain pond/seal matching problem, a transcritical free-surface solver, or an established rule for nonlinear branch selection.

The identities and control-theory arguments below are established mathematics applied as quality gates. They are not claimed as historically new theorems.

## 1. The specific open problem and the prototype's boundary

Schoof et al. (2023), section 3.5, equations 3.17–3.18, use a steady inner flux relation while discussing selection on a faster time scale. They explicitly explain that the inner bed equation then requires a time derivative and that discharge need not depend only on instantaneous head and outer slopes. Appendix B supplies the spatial equations. This is a concrete modeling gap, but identifying it does not solve it. [Primary full text](https://www.cambridge.org/core/journals/journal-of-fluid-mechanics/article/drainage-of-glacier-and-ice-sheet-surface-lakes/6A699229F6F8C928BA52E83ABFEC16B7)

A finite linearized operator can show how retained shape state affects flux and whether a steady closure misses that response. It cannot by itself establish the correct far-field conditions, moving-seal speed, nonlinear basin boundaries, hydraulic-jump admissibility or environmental transfer. Calling it a completed resolution of the paper's open problem would be premature.

A static map is already causal in the systems sense. Its possible deficiency is omitted state/history. An instantaneous spatial hydraulic problem also does not represent finite-speed hydraulic waves merely because its bed equation is time dependent.

The subsequently derived MATCHING_THEOREM.md addresses one of these gates: the canonical hydraulic normalization converges into the infinite pond, with a persistent shape sensitivity, and a specified full-line **linear fixed-datum** evolution is well posed. This narrows the remaining gap; it does not complete nonlinear moving-control matching or establish physical branch selection.

## 2. Positive thermal lag cannot cause a first Hopf instability

Consider

\[
\dot z(t)=-dz(t)+g\int_0^\infty z(t-s)K(ds),
\qquad d>0,\quad g\ge0,
\]

where K is a finite nonnegative measure of total mass M>0. Assume finite delay, or an appropriate exponential moment when making exponential-stability claims. The characteristic function is

\[
f(s)=s+d-g\widehat K(s).
\]

If d>gM, then for Re(s)>=0,

\[
|s+d|\ge\operatorname{Re}(s)+d>gM\ge g|\widehat K(s)|.
\]

No characteristic root lies in that half-plane. An exponential moment permits a small lambda>0 satisfying d-lambda>g integral exp(lambda s)K(ds), giving an exponential comparison bound. Mere integrability of a heavy-tailed kernel must not be promoted to exponential decay.

If gM>d, f(0)<0 and f(r) is strictly increasing to infinity for real r>=0. There is a unique positive real root r*. Every complex characteristic root with real part r>=0 satisfies

\[
r+d\le g\widehat K(r),
\]

so r<=r*, strictly when its imaginary part is nonzero. At gM=d, the crossing is s=0; if the first moment is finite,

\[
f'(0)=1+g\int_0^\infty sK(ds)>0.
\]

Thus increasing positive loop gain first loses stability through a real root, not a Hopf pair. Passive scalar positive thermal lag is insufficient to establish a new oscillatory onset. Related dominant-real-pole results are established in positive-system theory. [Ebihara et al. 2014](https://doi.org/10.1049/iet-cta.2014.0375)

This gate applies when the relevant feedback kernel is actually nonnegative, for example a positive inlet-temperature perturbation driving heat transfer in a suitable passive heat problem. Geometry and discharge derivatives can have sign-changing kernels. Their signs must be derived; positivity of temperature does not prove positivity of every linearized coupling. Conversely, a complex eigenpair in a multi-state model does not prove a first oscillatory instability.

## 3. Moving-bed mass and momentum benchmark

Let r=rho_i/rho_w, let E>0 be downward melting speed, and write b_t=-E and h=eta-b. For a fixed-width reach, exact depth-integrated liquid mass gives

\[
h_t+q_x=rE,\qquad
\eta_t+q_x=(r-1)E.
\]

Retaining q_x=rE while discarding h_t is inconsistent at the order at which bed melting changes water storage. In a slow-time expansion,

\[
q_{1X}=rE-h_{0T}=(r-1)E-\eta_{0T}.
\]

For a flat uniform reach, prescribed constant E, no spatial gradients or friction, and zero horizontal momentum carried by added meltwater, the conservative balances are

\[
h_t+q_x=rE,\qquad
q_t+\partial_x(q^2/h+gh^2/2)=-ghb_x.
\]

They admit the exact solution

\[
h=h_0+rEt,\quad b=b_0-Et,\quad q=q_0,
\quad\eta=\eta_0+(r-1)Et,
\]

\[
u=\frac{q_0}{h_0+rEt},\qquad
Fr^2=\frac{q_0^2}{g(h_0+rEt)^3}.
\]

At equal densities the free surface remains fixed and Fr decreases. A claimed injection-driven spatial choking effect that fails this case is a storage error. This benchmark assumes a specified injection-momentum convention; changing that convention requires changing the exact solution consistently.

Equal density does not eliminate moving-geometry or entrainment-momentum effects in a nonuniform reach. Near Fr=1, every term promoted by a distinguished asymptotic limit must be retained consistently, including slow storage, momentum and the moving boundary. A small melt source divided by a small critical denominator is not sufficient evidence of a new mechanism. This benchmark is intended for a physical melt-enabled solver, not as a direct test of a linear inner model that does not retain meltwater mass injection.

## 4. Uniform normal-flow cancellation

An exact control is available for the paper's hydraulic equations. Write

\[
\partial_X\left[B+\frac{Fr^2V^2}{2}+(q/V)^\beta\right]
=-q^{\alpha-1}V^{3-\alpha}.
\]

Linearize at uniform V0,q0, with R0=q0^(alpha-1)V0^(3-alpha)=-p0 and

\[
d=\frac{\beta(q_0/V_0)^\beta}{V_0}-Fr^2V_0>0,
\quad r_v=\frac{(3-\alpha)R_0}{V_0},\quad a=r_v/d,
\]

\[
\ell=\frac{(1-\alpha)V_0}{(3-\alpha)q_0}.
\]

For perturbations v,b,j and downstream normal matching v(L)=ell j, set z=v-ell j. Direct linearization gives

\[
-dz_X+b_X+r_vz=0,\quad z(L)=0,
\quad z(X)=-\frac1d\int_X^L e^{-a(s-X)}b_X(s)\,ds.
\]

With b(0)=0, upstream relative-head perturbation H satisfies

\[
j=\frac{H-\int_0^L e^{-as}b_X(s)\,ds}{D_q},
\quad
D_q=\frac{2\beta(q_0/V_0)^\beta+(1-\alpha)Fr^2V_0^2}
{(3-\alpha)q_0}>0.
\]

Now linearize the correctly referenced bed equation

\[
B_T+cB_X=cp_--V^3+V_-^3.
\]

For fixed c and p_minus, uniform normal matching gives v_minus=ell j, so

\[
b_T+cb_X=-3V_0^2(v-\ell j)=-3V_0^2z.
\]

The uniform discharge-dependent incision cancels against reference-bed translation. Zero initial shape perturbation therefore remains zero under head forcing: j=H/Dq. Initial shape can affect discharge; a forced shape-memory mode does not arise at first order in this matched uniform control.

This cancellation uses a moving reference and matching conditions. A fixed absolute upstream bed, a fixed absolute downstream head, or an omitted reference-incision term defines a different BVP. A fixed-frame prototype need not return zero absolute bed motion. It must instead demonstrate that the transformed shape and head obey the cancellation when that matched control is represented. Do not impose the benchmark on incompatible boundary conditions, or use incompatible conditions to advertise a physical memory mechanism.

## 5. Requirements before a predictive claim

- **Gauge covariance:** adding a common vertical datum shift changes neither discharge nor physical head difference. Translating the bed reference requires transforming lake head and its evolution equation. Report absolute level, relative head, seal height and shape separately.
- **Matching:** upstream bed perturbation zero is a finite truncation choice. Increase domain extent and compare admissible upstream/downstream matching prescriptions. A pole caused by the crop boundary is not a glacier mode.
- **Admissibility:** retain positive area and velocity, a stated margin from criticality, and a nonzero hydraulic Schur complement. Linear results outside that region require a different model.
- **Refinement:** converge spatial operators, time integration, nonlocal quadrature, eigenvalues and input-output kernels independently. Compare the linearized operator with small-amplitude nonlinear responses under identical boundaries.
- **Selection:** a stable eigenvalue calculation on one branch does not establish which branch a nonlinear initial condition selects. That requires evolving the nonuniform seal and resolving its basin/matching problem.
- **Budgets:** identify hydraulic dissipation, sensible heat, latent heat, ice conduction and injection momentum. A dissipative hydraulic identity is not a complete thermal/moving-ice energy budget.

Thermal meander instability, ripple instability and shape-evolving conduit models already exist: [Karlstrom et al. 2013](https://doi.org/10.1002/jgrf.20135), [Camporeale et al. 2017](https://doi.org/10.1103/PhysRevFluids.2.053904), [Dallaston & Hewitt 2014](https://doi.org/10.1063/1.4892389). These were screened as close alternatives; no exhaustive priority certification is asserted. A stronger contribution here would be a derived, converged state-dependent seal closure with a consequential forecast that established steady closures miss. The present prototype is a foundation for investigating that claim, not evidence that it has already been achieved.
