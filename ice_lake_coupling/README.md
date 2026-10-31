# Executable lake–channel coupling foundation

Python 3.9+ standard library only. No installation, network, training, GPU or CFD software is needed to run this module from the repository root.

```sh
python3 -m unittest discover -s ice_lake_coupling/tests -v
python3 -m ice_lake_coupling --case ice_lake_coupling/cases/synthetic_lake.json --output /tmp/my_lake_comparison
python3 -m ice_lake_coupling.verify --output /tmp/my_lake_refinement.json
```

Choose unused output paths. Runs never overwrite earlier results. The first command tests the equations and API. The second produces a fixed-control/thermal-control comparison, full CSV time series, JSON summaries and a hashed input manifest. The third checks timestep refinement on the synthetic example. All commands work without the previous research ZIP.

## Implemented

* A positive, bounded stage–area table with exact piecewise-quadratic storage and an analytic inverse. This extends the repository's analytic power-law storage option; it does not establish a new physical law.
* Extensive lake volume and sensible-energy states, in m³ and J, with temperature relative to the melting point.
* An implicit, positivity-preserving example lake step for a monotone stage-discharge law. No arbitrary wet/dry depth is inserted.
* An integrated water/enthalpy exchange object for eventual coupling to a channel solver, with signed backflow support and explicit rejection of unavailable mass/heat removal.
* A runnable synthetic forward model using an ideal critical rectangular outlet and a prescribed heat-transfer coefficient. It predicts lake stage, temperature, release, outlet heat export and control lowering. Meltwater is added downstream using water density; ice loss uses ice density.
* Input validation, dry/empty behavior, explicit failure outside table/phase-change scope, event-aligned forcing, a geometry-step restriction, conservation reporting and preserved outputs.

## What the executable is not

The demo is **not resolved CFD**, a trained PINN, a field-validated glacier forecast, a hydrofracture model or a hazard warning system. It contains no turbulent heat-transfer discovery. Its smooth critical-control law assumes a free, unsubmerged outlet. Channel storage, backwater, stratification, sediment, freezing, atmospheric evaporation and ice deformation are not represented. The atmosphere input is a prescribed net power, not a computed energy-balance model.

The generic port API is ready for an adapter; integration with Claude's new CFD solver has not been implemented or tested. The existing repository already has lake storage and sensible-heat accounting. This package independently tests a coupling interface and a table-based geometry extension without rewriting that active code.

## Equations and signs

With liquid temperature excess θ above the interface melting point, lake state is `(V, E)` with `E = rho_w cp V θ`:

```
dV/dt = I - Q
dE/dt = rho_w cp I θ_in + P_external - rho_w cp Q θ_lake
```

For one backward-Euler step, fixed geometry and monotone rating `Q(H)`:

```
V_new + dt Q(H(V_new)) = V_old + dt I
E_available = E_old + dt (rho_w cp I θ_in + P_external)
E_new = E_available V_new / (V_new + dt Q)
E_port = E_available dt Q / (V_new + dt Q)
```

The volume equation is bracketed within the measured storage domain. The splitting is first-order: the sill stays fixed during a lake step and is lowered afterward. Conservative budgets do not establish time accuracy; the refinement command assesses that separately.

For the demo outlet, `Q=Cd b sqrt(g) [2(H-z_sill)/3]^(3/2)` for positive head. Over its prescribed contact area `A_c`, the plug-flow heat fraction is `1-exp[-htc A_c/(rho_w cp Q)]`. The removed sensible heat becomes latent heat once. Ice volume is `E_latent/(rho_i L)`, added water volume is `E_latent/(rho_w L)`, and bed lowering is `ice_volume/A_c`. The output temperature divides remaining enthalpy by the *total* outgoing liquid volume including meltwater.

Heat transfer changes both geometry and reservoir release. Consequently, cumulative downstream heat export can increase even while water loses heat at the ice boundary, because more lake water is released. Check the full ledger before interpreting the sign of a comparison.

## Inputs and outputs

`cases/synthetic_lake.json` documents the schema through explicit SI-unit keys. Its dimensions and thermal coefficient are synthetic, not measured Himalayan inputs. In `forcing`, each entry applies until `until_s` and is piecewise constant; the final interval must end exactly at `duration_s`. `external_power_W` may be negative only while enough sensible heat remains. Negative temperature and unmodeled freezing are rejected.

`timeseries.csv` reports accepted-step mean fluxes and end-step lake states. At the initial instant, channel-outlet temperature/heat are blank because no step has been integrated. The sampled peak-release metric mixes an initial instantaneous value with step-mean values; refine time before using it as a peak estimate. There is no claim of substep peak detection.

Use `--dt 15` to rerun at a smaller maximum step; the actual override is recorded in `manifest.json`. A failed run records `FAILURE.json` when the output manifest exists. The geometry-step restriction controls update size, not a proven integration error tolerance.

See `INTEGRATION_HANDOFF.md` for port coupling and `EXTENSIONS.md` for the decisions on PINNs, hydrofracture, debris, glacier mechanics and hazard routing.
