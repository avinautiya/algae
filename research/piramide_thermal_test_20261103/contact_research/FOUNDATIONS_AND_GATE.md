# Contact-resistance foundation and identification gate

This is an isolated, additive research component. The existing four field
comparisons and their failed all-baseline superiority verdict are unchanged.
No Tower2 camera outcome is opened, no new production experiment is authorized,
and no intrinsic contact resistance or missing-model flaw has been established.
Contact/porous-interface resistance is existing physics; historical originality
is a separate prior-art question. A successful unit test is not a field discovery.

## Definition and conservative discretization

Let `T_d` be the debris-side contact-face temperature, `T_i` the ice-side phase
interface temperature, and`q_d` the downward heat flux. Introduce a zero-capacity
contact resistance`R>=0`, units m²K/W:

```
T_d − T_i = R q_d.
```

This is separate from bulk debris thickness/conductivity. For a final debris
cell center a half cell above the face, eliminate the face node:

```
G_d = 1 / [Δz_d/(2 k_d) + R]
q_d = G_d (T_d,last − T_i)
G_i = 2 k_i/Δz_i.
```

The upper debris boundary retains its own`2k_d/Δz_d` conductance; addingR does
**not** alter that top half-cell resistance. First try the melting active set
`T_i=0C`, independently solve the two implicit material blocks and evaluate
`q_net=q_d−q_i`. If`q_net>=0`, mass is`m=dt*q_net/L`, nonnegative. Otherwise
resolve the full coupled implicit system using

```
G_interface = [Δz_d/(2k_d) + R + Δz_i/(2k_i)]^−1,
T_i = (G_d*T_d,last + G_i*T_ice,first)/(G_d+G_i),
m=0.
```

The last equation reconstructs the **ice-side** interface; the debris face is
`T_d=T_d,last−q_d*Δz_d/(2k_d)`. The contact has no heat capacity, so it cannot
create/remove energy. With signed deep-boundary heat, the exact discrete ledger
remains

```
ΔE_debris + ΔE_ice + Lm = Q_surface − Q_deep_out.
T_i<=0C, m>=0, T_i*m=0.
```

`R=0` recovers the frozen model exactly, tested on warming and cooling trajectories.
`R=+infinity` is implemented explicitly as zero interface conductance, without
undefined`infinity*0`. It is an insulated limit, not an identified field value.
For warm stationary ice and constant positive surface temperature,

```
q = Ts / (d/k + R),  T_d = R q,  m=q*dt/L.
```

For a cold stationary composite ending at0C deep ice,
`q=Ts/[d/k+R+d_i/k_i]`; analytic profiles and both heat ledgers are independently
tested. Wet-debris latent heat, moving geometry and ice advection remain absent.

## Entropy identity and inequality

Use absolute Kelvin temperatures`theta_d,theta_i>0`. The contact entropy
production is

```
sigma_R = q_d(1/theta_i − 1/theta_d)
        = R*q_d²/(theta_d*theta_i) >= 0.
```

It stays nonnegative for cooling/upward heat as well as warming. Every conductive
edge between positive temperatures has production
`G(theta_a−theta_b)²/(theta_a*theta_b)>=0`. Include both boundary half cells,
material interior edges, the debris half cell, contact and ice half cell. For
constant cell heat capacityA, define sensible entropy`A log(theta/theta_ref)`.
Backward Euler contributes the additional **numerical**, nonnegative defect

```
D_BE = A [log(r) − 1 + 1/r] >=0,  r=theta_new/theta_old>0.
```

The inequality follows because`f(r)=log(r)−1+1/r` has derivative`(r−1)/r²` and
its global minimum is zero at`r=1`. Summing the implicit energy equations divided
by each new cell temperature gives the discrete entropy identity

```
ΔS_sensible + Lm/Tm − S_boundary_in
  = dt*sum(sigma_conduction_edges) + sum(D_BE).
```

`S_boundary_in=Q_surface/T_surface−Q_deep/T_deep`. At melting,`Tm=273.15K`;
on the cold active set`m=0`. Internal conductive/contact exchanges cancel
algebraically. The implementation records the identity residual separately
from the nonnegative physical production and BE defect. It never describes
the BE defect as extra physical entropy. This is an established passive-diffusion
inequality, not a priority claim to a new mathematical theorem.

## Structural identifiability

For an ice interface held at0C, debris temperatures depend on

```
alpha = k/C,  ell=kR,  Bi=d/ell.
```

For any`s>0`, the transformation`(k,C,R)->(s*k,s*C,R/s)` leavesalpha andell
unchanged, hence the entire debris temperature history unchanged, while heat
flux and melt scale bys. An executable transient test constructs this exact
ambiguity, including known ice properties when the ice remains temperate.
Temperature observations alone therefore do **not** identify absolutek,C,R or
absolute melt in that regime. A known cold-ice domain can break the exact scaling
during informative cooling transients because the relative material effusivity
changes. That possibility is not proof of practical identification from noisy
thermistors, assumed ice initial state and uncertain geometry.

The steady profile is`T(z)/Ts=1−z/(d+ell)`. Thusd andell also have an exact
steady degeneracy: only their sum is visible at fixed observation positions.
For example`d=.50,R=.10,k=1` and`d=.52,R=.08,k=1` give identical interior
steady temperatures. Unknown basal gap, unresolved low-conductivity material,
sensor depth/contact error and genuinely intrinsic contact resistance can all
look like an effective boundary resistance. A fittedR is not a causal diagnosis.

## TRAIN-only falsification before fitting a new candidate

Under the explicitly restricted **steady, one-dimensional, temperate-interface,
positive-flux** assumptions, a sensor at or above the debris contact satisfies

```
T_sensor = q[R + delta/k] >= qR.
R <= (T_sensor + declared_temperature_error)/q.
```

The previously inspected Tower1 TRAIN diagnostic supplied by the parent review
has meanDT5 approximately.374C, declared sensor error.4C and heat flux about
23.5W/m². These give necessary`R<=.03294 m²K/W`; anR around.371 needed to explain
the full amplitude mismatch would require roughly8.72K contact jump. Therefore
**contact resistance alone cannot explain that full mismatch under these stated
assumptions**. The bound does not apply unchanged to transient/cold-interface
periods and is not a new intrinsicR measurement. It is a falsifier, not a reason
to force another fit. Do not transfer a fictional3m stake separation: the final
paper places the stake and thermistors in the same pit; PVC stake length is not
their separation.

## Source-gated inference plan — currently not a production authorization

1. Preserve the already inspected four periods as failed model-superiority
   tests. Pin the prior-art record before claiming any method is new.
2. Define a genuinely untouched target site/event and audit native camera
   operator, clocks, resets, geometry and phase assumptions without reading its
   melt values. Existing Tower2 feature availability does not establish this
   model's validity.
3. Pin a TRAIN-only thermistor block and initialization; fit only identifiable
   temperature parametersalpha andell over a finite predeclared grid, for each
   independently source-bounded geometry and bulk-capacity scenario. All imposed
   boundary sensors are excluded from loss. Do not fitR andC independently and
   pretend a unique optimum is identified. Define exact ties and edge/flat-profile
   failure rules in advance.
4. Establish practical identifiability before production: joint profile losses,
   sensitivity/Jacobian rank, held-out **thermal** predictions, geometry/capacity
   ambiguity and the constantR necessary bounds. A narrow melt-conditionedR fit
   cannot substitute for these checks. Record a no-go if the TRAIN profiles reject
   the resistance needed by the hypothesis; do not alter the hypothesis by
   secretly adding a target-calibrated scale.
5. Conditional conversion is`k=alpha*C`, `R=ell/k`; report every declaredC/geometry
   scenario. Absolute physical melt requires independent bulk-capacity and ice
   state support or an explicitly separate empirical TRAIN-melt calibration.
   A latter gain is not a material-property measurement or proof of contact physics.
6. Freeze all future predictions and comparison rules before opening new camera
   targets. Include a credible heterogeneous-layer/geometry alternative and the
   strongest feasible existing baselines, rather than only comparing withR=0.

No Tower2 production run or camera unseal is justified by this foundation alone.

## Execute verification

From`contact_research/`:

```
python3 -m unittest discover -v
```

The solver and standard-library tridiagonal implementation are local copies;
no original model, frozen source, result or criterion is edited. Tests import
the frozen model read-only solely for theR=0 equivalence oracle.
