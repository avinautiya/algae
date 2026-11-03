# Source-gated transient debris conduction

This executable component adds **thermal storage and delayed basal heat delivery**
to a fixed-thickness sub-debris melt calculation. It implements established heat
conduction; neither the PDE nor this discretization is new mathematics. It has
analytic and numerical verification, **no field-accuracy result yet**. No Arolla
target, calibration coefficient, freeze, historical output or held-out score is
changed or loaded by this package.

## Physics and discrete accounting

For downward coordinate `z` on `[0,d]`, constant conductivity `k` and volumetric
heat capacity `C`, solve

```
C ∂T/∂t = k ∂²T/∂z²
T(0,t) = measured/independently specified DEBRIS SURFACE temperature Ts(t)
T(d,t) = 0 °C
T(z,0) = explicitly supplied cell-center initial profile
```

Use `N` uniform control volumes (`Δz=d/N`). A cell stores `C Δz T` joules per
square metre. Interior face conductance is `k/Δz`; each Dirichlet boundary is
**a half cell away**, so its conductance is `2k/Δz`. A backward-Euler step is

```
(C Δz / Δt)(Tj_new − Tj_old)
  = Gj−1/2(Tj−1_new − Tj_new) + Gj+1/2(Tj+1_new − Tj_new).
```

The first and last adjacent temperatures are the prescribed boundaries. The
tridiagonal system is solved with Thomas elimination using Python's standard
library. No external numerical package is required. Internal face contributions
cancel algebraically, giving the **same-discretization ledger**

```
ΔE_storage = Δt (q_surface − q_basal)
q_surface = (2k/Δz)(Ts − T0_new)
q_basal = (2k/Δz)Tlast_new.
```

Surface heat is signed: warm debris may lose energy through a cooled surface.
Initial sensible energy is counted. It is therefore valid for basal heat to
exceed contemporaneous surface supply during cooling. The ledger is roundoff
closure of the discrete equations, not proof of physical accuracy or zero
truncation error. Both global and per-step residuals are returned.

All basal heat is assigned to latent melting of ice already at 0 °C:
`melt_mm_water_equivalent = basal_heat_J_m2 / latent_heat_J_kg`.
One kg/m² equals one mm water equivalent; **no extra factor of 1000** is needed.
The declared latent heat defaults to `333700 J/kg`. This is a model assumption,
not a measured site parameter. Negative initial/boundary temperatures are
refused because refreezing and cold-ice energy demand are not implemented.

## Inputs and source gates

`Column` requires positive finite `conductivity_W_mK`,
`volumetric_heat_capacity_J_m3K`, `depth_m`, a nonempty explicit cell-center
`initial_temperature_C` list, and five source records keyed `conductivity`,
`heat_capacity`, `depth`, `initial_profile`, `basal_ice_state`. Resolution is the profile length.
`SurfaceBoundary` requires contiguous hourly `hourly_surface_temperature_C` and
a source record. Value `i` applies exactly on relative time
`[3600*i, 3600*(i+1))`; the caller must establish the absolute start, source clock,
interval convention and any visit-time scenarios before supplying this array.
The solver shortens steps at every hourly transition and at the final horizon.
There is **no interpolation, gap filling, source-clock inference or spin-up**.

Every record needs nonempty `source`, `location`, `method`, status
`independently_measured` or `independent_literature`, and explicit
`uses_evaluation_melt: false`. Target-fitted and unknown records are refused.
The surface record additionally needs
`quantity: "debris_surface_temperature"`; an air-temperature field is not an
adequate boundary. For analytic illustrations only, status `declared_synthetic`
is accepted if and only if `synthetic: true` is supplied. Declarations are
auditable requirements, **not machine verification of source truth**. Parameters
derived from the same evaluation melt must not be relabelled as independent.

`horizon_s`, `max_horizon_s` and `dt_s` are explicit positive finite options.
The horizon may not exceed the declared maximum or available surface forcing.
The maximum is a caller-supplied computational scope limit, not an automatically
established physical validity duration. All numeric booleans, NaN/Inf, invalid
coefficients, missing profile/boundary data and unrepresentable arithmetic are
refused. Accuracy requires independent spatial/time refinement, even though
backward Euler is stable for large timesteps.

## Execute

From the parent `model/` directory:

```
python3 -m unittest thermal_memory.test_solver -v
python3 -m thermal_memory thermal_memory/example_synthetic.json /tmp/thermal_example_NEW.json
```

The CLI refuses an existing output path. The example's two explicitly synthetic
initial-state scenarios demonstrate different stored heat under identical
surface forcing; their values are not tuned to field observations. Python API:

```
from thermal_memory import Column, SurfaceBoundary, simulate, ensemble
result = simulate(column, boundary, horizon_s=7200, max_horizon_s=7200,
                  dt_s=30, synthetic=False)
```

`ensemble()` uses the same boundary/options for a finite explicit list of
parameter/initial-profile scenarios. It returns each result plus min/max melt.
This is **only an envelope of the declared scenarios**, not a statistical
confidence interval or a proven bound over all intermediate parameters and
profiles. To test different boundary histories, call `simulate` separately and
retain their declared chronology and provenance.

## Verification and application limits

Tests independently verify: the exact linear steady solution (including one
cell), its constant heat flux and melt units; discrete sine-mode amplification;
second-order joint mesh/time refinement against the continuum sine mode;
first-order time refinement against exact semidiscrete decay; exact discrete
storage/flux accounting during warming and cooling; hour transitions and partial
final steps; source gates, missing forcing, unsupported refreezing, bad inputs;
initial-state memory, delayed basal flux, finite-scenario bounds; and executable
CLI plus overwrite refusal.

This is a fixed-geometry **conductive sub-debris column**, not the stream,
free-surface, vertical-cliff or grain solver. No moving thickness, advective pore
water, debris heat-capacity variation, snow physics, atmospheric energy balance,
solar penetration, lateral conduction, ice cold content or refreezing is
included. The imposed surface boundary represents external energy exchange only
if independently observed or physically supplied. Unknown surface temperature
or heat capacity cannot be silently replaced with air temperature or tuned to
the inspected melt target. Integration into a new predictor requires a new
preregistered comparison and independent source-supported boundary/initial
conditions; passing these verification tests does not establish novelty,
environmental superiority or field validation.

The basal-ice source must independently support
`state: "temperate_melting_interface"`, `temperature_C: 0`; labelling an
unobserved glacier interface as temperate is not evidence. The capacity source
must specify `quantity: "bulk_debris_volumetric_heat_capacity"`. Solid-rock
`rho*cp` alone does **not** supply bulk debris heat capacity: void fraction and
pore moisture matter. Use an independently measured effective bulk capacity or
document independent constituent properties, porosity and moisture in a
separately audited mixture conversion. The solver refuses a source explicitly
labelled rock-only but cannot establish the truth of a caller's bulk declaration.
