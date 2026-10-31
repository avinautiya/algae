# Thermal observability of local ice erosion

This independent contribution contains a continuum proof, a resolved transverse heat model, a sharp conditional forecasting tool and a reproducible synthetic study. It changes no existing CFD solver, build, job, checkpoint or frozen result.

The question is specific: **when do bulk thermal measurements determine ice-directed heat in a particular downstream reach, and what is the best prediction if they do not?**

## Mathematical contribution

`THEOREMS.md` proves a necessary-and-sufficient observation-placement criterion for a prescribed unidirectional Graetz flow, including no-slip laminar shear. It constructs two smooth positive inlet profiles matching all the specified bulk readings yet producing different regional heat when a required endpoint is unobserved. It then derives a conservative, computable sharp range with uncertain measurements and optional surface heating for a transverse finite-volume model.

Graetz thermal development, inverse thermal ambiguity, adjoints, linear-programming duality and interval minimax reasoning are established. The candidate contribution is the specific placement characterization plus constructive regional-incision recovery. Mathematical priority has not been certified; see `PRIOR_ART.md`. A successful calculation does not establish a new thermal-memory mechanism.

## Run

Tested with Python 3.9.6, NumPy 2.0.2, SciPy 1.13.1, Matplotlib 3.9.4. The forecasting core uses NumPy/SciPy; plotting additionally needs Matplotlib. Use an isolated environment if dependencies are not already installed.

```sh
python3 -m unittest discover -s ice_thermal_observability/tests -v
python3 -m ice_thermal_observability demo --output /tmp/thermal_case_01.json
python3 -m ice_thermal_observability predict --case /tmp/thermal_case_01.json --output /tmp/thermal_forecast_01
python3 -m ice_thermal_observability.study --output /tmp/thermal_study_01
python3 -m ice_thermal_observability.audit_fixed_data --output /tmp/thermal_study_01/fixed_data_refinement.json
python3 -m ice_thermal_observability.figures --study /tmp/thermal_study_01 --output /tmp/thermal_figures_01
```

Each command preserves existing output paths. The demo is explicitly synthetic. The study copies its predeclared protocol and source hashes before running. All baselines consume the same bulk observations. Comparisons report every result; midpoint minimax is not promised to win average prediction error.

The main grid study regenerates bulk data on each grid from the same uniform inlet profile. The supplemental fixed-data audit holds the observation vector constant across meshes; it is a numerical-validation addition, not a changed test campaign. Plug data are analytical; shear data use a 768-cell synthetic reference with a separately reported refinement indicator. Neither is field validation.

## Physical model and inputs

Water satisfies `rho cp u(y) theta_x = k theta_yy` over a fixed temperate bed with an insulated or heated top. Velocity is either constant plug flow (explicit analytical reference) or a prescribed no-slip, shear-free parabolic shape. Both lie in the continuum theorem's velocity class; the explicit sine-mode witness is for plug flow. No momentum equation is solved. Axial diffusion, time-dependent thermal storage, turbulence, buoyancy, cold-ice conduction, free-surface evolution and meltwater feedback are absent.

The JSON supplies initial temperature bounds; bulk-temperature observations with absolute error bounds and provenance; channel properties; named downstream reaches; and a horizon over which this steady scenario is assumed to persist. Optional nonnegative net surface-heating intervals supply additional uncertain inputs. Observations strictly downstream of the forecast origin are rejected. Future weather or future outlet data are not silently used as forecasts.

Outputs include sharp finite-volume power/lowering ranges, scalar minimax forecasts, attaining initial/source profiles, direct conservative replay, and evaluated dual bounds. Temperature/profile bounds and measurement errors are assumptions; they are not learned environmental uncertainty. The geometry-change screen is user-declared, not a theorem about model applicability.

Initial-profile boxes represent possible upstream inlet boundary fields. The positive semigroup makes their downstream states smooth and nonnegative. This does not show that an arbitrary glacier meteorological history generates every admitted profile. Restricting the class to independently justified reachable histories may narrow the range.

## Interpretation

* Continuum theorem: a proved implication for the stated prescribed unidirectional-flow PDE.
* Sharp LP range: exact for the declared discrete model/class, checked by primal and dual witnesses.
* Grid study: numerical convergence evidence, not a certified enclosure of PDE/model error.
* Baseline comparison: synthetic evidence under the declared profile distribution, not glacier validation.
* Fixed-flow heat-to-lowering conversion: a conditional diagnostic, not a coupled moving-ice simulation.

The practical use is to ask whether thermal calibration supports a local erosion prediction and what thermal information would reduce its uncertainty. See `HANDOFF.md` for later use with Claude's verified solver.
