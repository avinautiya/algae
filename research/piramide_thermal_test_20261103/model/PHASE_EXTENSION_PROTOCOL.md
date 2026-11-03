# Version2 cold-ice phase extension before camera outcome access

Version1's failed no-refreezing phase gate is preserved. Its signed basal heat
cannot be renamed melt. Version2 changes the physical model **before any camera
outcome is opened**, using a declared ice cold-content domain and an active-set
melting interface. The mathematical components are established conduction and
complementarity, not a newly discovered physical law.

## Domains, phase condition and accounting

Retain the buried upper thermistor as a conditioned temperature boundary and the
three source-bounded static ice/sensor geometries. Couple debris conduction to a
2m ice column with fixed `ki=2.1 W/m/K`, `Ci=917*2100 J/m³/K`. Initial and deep
ice are0C assumptions. Debris is initialized from the same measured initial
profile and an assumed0C ice boundary. Sensible cold content develops explicitly.

At each backward-Euler step first solve debris and ice with interfaceTb=0C.
Write `qd` for heat entering the interface from debris, `qi` for heat leaving it
into ice. If `qd-qi>=0`, melt is `(qd-qi)*dt/Lf` kg/m². Otherwise solve **both
domains again**, using the series harmonic conductance of the interface half
cells. Recover the common interface temperature by conductance weighting; it is
below0C and melt is exactly zero. Thus `Tb<=0`, `m>=0`, `Tb*m=0`. No negative
heat, temperature or observed value is clipped, and no refreezing water is
invented. Signed cooling is stored as sensible cold content in ice/debris.

The exact discrete energy balance is

```
ΔE_debris + ΔE_ice + Lf*m = Q_surface − Q_deep_out.
```

The deep boundary heat can be negative: a0C deep boundary may supply heat to
cold ice. Ledger closure verifies the numerical accounting, not the boundary's
environmental accuracy. Both global and maximum per-step residuals are reported.

Geometry remains fixed. There is no interface retreat, ice thermal advection,
debris pore-water latent heat, moisture advection or changing material thickness.
Warm air/surface forcing is observed, not forecast. The ice extension repairs
the particular interface cooling defect; it does **not** make wet debris fully
phase-complete. Field accuracy remains conditional on these approximations.

## Predeclared input and family

Reuse the already published version1 **training-temperature-only diffusivities**
for each geometry/hold/method, without refitting or consulting camera outcomes.
This keeps the numerical/physical intervention explicit. Windows remain fixed:
Dec2 03:00→Dec16 03:00 TRAIN; Dec16 03:00→Dec30 03:00 EVALUATION2022. Preserve
native naive local clocks and30-minute forcing intervals, both left/right holds,
all three geometries and both forward-fit/derivative-OLS alpha methods.

The independently decoded TRAIN-only capacity audit supplies three sensible
capacity scenarios: dry, training-probe proxy and water-saturated. Porosity0.44
is a source assumption; the probe proxy is **not** a depth-weighted capacity
measurement. Take exact values from `data/TRAINING_CAPACITY_AUDIT.json` and
set `k_debris=alpha*C_scene`. Primary capacity is the TRAIN-probe proxy; primary
geometry remains middle_gap_20mm, left hold, forward-fit alpha. All36
geometry×hold×method×capacity results are retained. No favorable scene selection.

Use32debris cells,64ice cells, maximum300s steps, split at exact30-minute
boundaries. This is a declared discretization, not a numerical-accuracy guarantee.
Fixed-parameter refinement and ice-depth sensitivity remain necessary before a
strong environmental-accuracy claim. Parameters use source-based assumptions,
not additional camera tuning.

## Raw physical component versus empirical output calibration

Compute raw positive latent melt mass (kg/m²==mmwe) and its exact energy ledger.
For compatibility with the sealed analysis, return
`normalized_heat_K_m = raw_latent_energy_J_m2 / C_scene`. This is latent energy,
not the signed total debris/interface sensible flux. Cooling flux counts remain
reported separately and do not invalidate a model that explicitly stores it.

A later routine may fit one positive empirical gain on exactly the original
TRAIN-camera target, matching the baseline calibration budget. Freeze all later
predictions before opening EVALUATION-camera outcomes. Because the known ice
conductivity/capacity breaks debris scaling invariance, that gain is **not** a
new physical heat capacity; physically changingC would require a new coupled
solve. The energy ledger applies only to unscaled raw component output. A
gain-corrected prediction is empirical and cannot inherit the raw energy ledger.

Camera operator, endpoint roles and source checks remain unchanged: exact exposed
stake distance increments in metres ice times917 give mmwe. The ambiguous provided
`ablation` field is not used. No negative camera increment is clipped. There is
no independent camera uncertainty distribution or confidence interval.

## Freeze and execution

```
# From model/, after tests and before accessing any camera outcome:
python3 phase_experiment.py freeze PHASE_EXTENSION_SPEC.json --out frozen_piramide_v2
# Publish the spec, source files and their hashes; then:
python3 phase_experiment.py predict frozen_piramide_v2 --commit PUBLISHED_SHA40
```

The version2 freeze binds all new code/protocol/spec hashes, version1 prediction
and freeze hashes, capacity audit hash, unchanged process hashes and the identical
camera evaluation contract. Its prediction module never opens camera outcomes.
The later target comparison may fail. No novelty, glacier-wide superiority,
professional-model superiority or publication outcome is presumed.
