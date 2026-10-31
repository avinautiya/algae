# Offline forecast mathematics and environmental screening

This package advances the theory prepared for the glacier-channel project while the active CFD implementation continues independently. It edits no active solver, case, build script or job. It does not launch CFD, install software or change the user's milestone.

## Runnable contribution

Python 3.9+ standard library plus the project's `ice_cfd_checks` package from commit `65819f5cf101f485e061865b86a13f9a56d8c0b9`. Both packages are present on this branch. The delivery ZIP includes both, so it can be extracted into a fresh directory and run without downloading dependencies.

```sh
python3 -m unittest discover -s ice_forecast_theory/tests -v
python3 -m ice_forecast_theory demo --output /tmp/forecast_theory_demo_01
python3 -m ice_forecast_theory forcing --input /tmp/forecast_theory_demo_01/forcing_manufactured.json --output /tmp/forcing_result_01.json
python3 -m ice_forecast_theory regime --input /tmp/forecast_theory_demo_01/regime_synthetic.json --output /tmp/regime_result_01.json
```

Use new output paths. The examples are manufactured/synthetic and no environmental accuracy or discovery is asserted.

## What is implemented

**Budget-constrained forcing response.** The previous package analyzes an initial disturbance in a constant linear operator. This package supports interval-dependent operators and piecewise-constant forcing, exact interval propagators including singular operators, explicit physical-budget rows, a time-normalized forcing norm and a replayed optimal perturbation. It reports gain, input sequence, budget residuals and replay disagreement. It requires a validated physical operator to become a glacier-prediction analysis.

**Environmental regime screening.** Calculates freshwater density/expansion using a declared approximate EOS, dimensionless groups, mean stratification, local-versus-bulk inertia/buoyancy ratios, timescales and heat-power scales. It flags scope assumptions rather than pretending to certify an entire model. It is not a buoyant-flow solver or a melt closure.

## Input contracts

For `forcing`, supply:

* `intervals`: ordered objects containing positive `duration_s`, square `A_per_s`, and `B_per_s`. All intervals use the same state basis and input-channel ordering.
* `C`: a fixed endpoint-observable matrix.
* Optional `budget_rows_on_input_amplitudes`: rows acting on stacked interval/channel forcing amplitudes. Derive these rows from actual physical budgets; do not assume equal temperature means imply equal heat.
* Nonempty `state_normalization`, `input_normalization`, `output_normalization`, `operator_provenance` strings.

The solver uses zero initial perturbation. It supports at most 32 forcing coordinates and an augmented state/input matrix dimension of 32. It assumes a time-average-L2 norm and linear equality constraints; input positivity, box bounds, nonlinear heat-budget equality and uncertain operators need separate treatment. The reported bound applies to the supplied linear model only.

For `regime`, supply `temperature_C, depth_m, velocity_m_s, length_m, width_m, ice_position`. Temperature is on ITS-90; `ice_position` is `bed` or `ceiling`. Optional inputs are `local_patch` with `depth_m, speed_m_s, length_m`, nonnegative `potential_drop_m`, signed `external_power_W`, and named positive `properties` overrides. Unsupported salinity/pressure inputs are rejected rather than silently ignored. The EOS is a one-atmosphere screening approximation.

The mean-density stratification diagnostic does not resolve nonmonotone temperature profiles. Properties are declared constants for the screen. No regime or transition is certified by a universal Richardson or Rayleigh cutoff.

## Read the research, not only the software

* `THEORY.md`: interval propagation, budget constraints, norm consistency, adjoint connection and limits of nonlinear extrapolation.
* `PHYSICS_SCOPE.md`: process-by-process environmental assumptions and exclusion decisions.
* `RESEARCH_DECISION.md`: a focused, falsifiable physical mechanism and meaningful forecast endpoints.
* `LITERATURE_AUDIT.md`: new competing work, access depth and candidate-claim limits.
* `HANDOFF.md`: later integration without interrupting the current milestone.

No tool here proves that the model beats every existing glacier model, accounts for every process or meets a journal's publication standard. Its value is to make a specific important claim testable with stronger mathematical and physical controls.
