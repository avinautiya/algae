# Signed measured-gradient heat retrieval with explicit ice cold content

This module is a physical/numerical foundation before any new Tower2 camera
release. It fits no field parameter and scores no new target. Measured thermal
gradients and moisture are conditioning observations; the result is a retrieval
or conditional hindcast, not a forecast of future weather. The heat equation,
Stefan complementarity and P1 storage formulas are established mathematics.

## Ice component and phase condition

Ice has cell-centered temperaturesTi, constantki/Ci, depth2m by default, explicit
initial profile and a deep0C Dirichlet boundary. The upper interface has
temperatureGamma<=0C, nonnegative melt massM and`Gamma*M=0`. Its signed external
fluxq is positive into the interface. There is no zero clipping ofq or an
invented liquid-water reservoir for refreezing.

Without a soil interface state, at each implicit step tryGamma=0C and solve
the ice heat equation. Its incoming sensible demand is
`q_ice=Gi*(0−T_first_new)`, `Gi=2ki/dz`. Ifq>=q_ice, the excess produces
`M=(q−q_ice)*dt/L`. Otherwise impose the **signed Neumann flux** on the full ice
block and reconstruct`Gamma=T_first_new+q/Gi`. The interface cools below0C and
M=0. On reheating, this cold content must be restored before melt resumes.

For optional known *soil* interface areal capacityA (J/m²/K), the balance is

```
A*(Gamma_new−Gamma_old)/dt + q_ice + L*M/dt = q.
```

The melting candidate must also supply`A*(0−Gamma_old)/dt`. If it cannot, solve
Gamma and every ice temperature together implicitly. The first matrix row is
`(A/dt+Gi)*Gamma_new−Gi*T_first_new=q+A*Gamma_old/dt`. Ice rows include their
ordinary full finite-volume capacities. A is not an additional ice half cell;
it represents a separately derived soil basis contribution. Counting both an
ice half-cell lump and the existing first ice control volume would double count.

The exact discrete ledger is

```
Delta E_ice + A*Delta Gamma + L*M = Q_input − Q_deep_out.
```

Deep heat is signed:0C deeper ice can supply heat to a cooled domain. The ledger
does not establish that this particular deep boundary/initial profile is correct
in the field. Units: mass kg/m² equals mm water equivalent, with no extra1000.
Geometry, ice advection and pore-water latent heat are not represented.

## Storage below the measured-gradient plane

Let measured probes be`(z4,T4)`, `(z5,T5)`, withiceinterface at`d>=z5`.
A finite-difference gradient naturally represents the plane
`p=(z4+z5)/2`. Set`a=(z5−z4)/2`, `g=d−z5`. For a piecewise-linear profile through
T4,T5 andGamma, exact sensible storage below that plane is

```
E_below = Csoil * integral_p^d T dz
        = F(T4,T5) + A*Gamma,
F = Csoil*[a*T4/4 + (3a/4+g/2)*T5],
A = Csoil*g/2.
```

The first term integrates the lower half of the T4–T5 segment; the second
integrates the entire T5–Gamma gap. `soil_storage_components()` returns F andA
separately. If the chosen physical flux plane were insteadz5, the formula would
be`F=Csoil*g*T5/2`, `A=Csoil*g/2`; the two plane conventions must not be mixed.

For declared interval fluxq_plane and observed F endpoints, define
`q_eff=q_plane−(F_next−F_previous)/dt_interval`. Feedq_eff to the ice solver with
the separate areal capacityA. Summing its exact ledger gives

```
Delta(E_ice + A*Gamma + F) + L*M = Q_plane − Q_deep_out.
```

There is no second subtraction ofA*DeltaGamma and no duplication of the measured
soil storage. All F differences, including cooling/storage release, retain
their signs. Temporal interpolation/reconstruction of source samples must be
declared by the field driver. Missing values or unsupported source endpoints
must be rejected upstream. This module silently fills nothing.

Csoil must be an explicitly source-supported/bounded bulk capacity. Solid-rock
rho*cp alone is inadequate for porous/moist debris. If a driver changesCsoil
with time, it must difference the full energyF, not merely useC*Tdot. The
interface term also becomes`Delta(A*Gamma)` ifA changes, requiring a different
solver. This implementation assumes constantA. Neither a varying moisture-based
conductivity nor sensibleC accounts for pore freezing, moisture advection or
evaporative heat transport automatically.

## Intended moisture-conditioned flux law and inference limits

The parent may separately test a candidate

```
q_plane = kappa_ref*(theta_deep/theta_ref)^p*(T4−T5)/(z5−z4).
```

The coefficient is a learned effective heat-flux closure unless independent
conductivity measurements establish otherwise. Moisture units, calibration,
depth association, source uncertainty and handling of nonphysical/nonpositive
theta must be specified before a new test. Negative gradients remain signed.
Do not clip sensor water content or thermal gradients to manufacture a positive
model. A moisture-weighted PDD comparator must receive the same extra exponent
and learning budget.

Learningp from the four already released Tower1 periods is **development**;
it cannot convert those periods into fresh validation. A conductivity/gain fit
may use only its declared TRAIN-melt interval. Recompute the ice phase dynamics
for every trialkappa: a post-solve mass multiplier is not equivalent to changing
the imposed heat flux in the presence of cold content and soil storage.

For signed gradients, increasingkappa does not necessarily give pointwise larger
heat input, so global monotonicity/unique calibration is not automatic. A
predeclared bracketed bisection can locate a continuous root, but multiple-root
or sampled monotonicity checks must accompany any uniqueness claim. Never use
new evaluation melt to choose the bracket, p, moisture proxy or storage plane.

## Verification and use

From`moisture/`: `python3 -m unittest discover -v`.

Nine tests cover exact positive-flux latent units, signed cooling, independent
one-cell Neumann and two-node soil/ice analytic inverses, cooling/reheating delay,
exact P1 soil quadrature, complete plane-energy accounting, spatial refinement
against the semi-infinite Neumann cooling solution, and invalid inputs. No field
accuracy or historical novelty follows from these verification tests.

API:

```
from flux_ice import soil_storage_components, retrieve_from_plane
# Derive F_native and constantA from actual source geometry and bulk capacity.
result = retrieve_from_plane(
    plane_flux_W_m2=q_native, observed_storage_F_J_m2=F_native,
    initial_ice_temperature_C=ice_profile,
    initial_interface_temperature_C=0.,
    soil_interface_areal_capacity_J_m2K=A,
    interval_seconds=1800., dt_seconds=300.)
```

F has one more value thanq. Every returned physical heat/mass ledger is for the
actual simulated component; a separate empirical output correction cannot
inherit its conservation claim. All existing failed model-comparison results
remain unchanged. Tower2 outcomes must stay sealed until the parent publishes
the new source, protocol, parameter choices and predictions.
