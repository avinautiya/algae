# Executable glacier-stream recovery forecasting toolkit

Independent additive delivery, 1 November 2026. Claude's active files, jobs,
worktrees, checkpoints and frozen results are untouched.

## Delivered models and analysis

| Component | What it predicts or checks | What it does not establish |
|---|---|---|
| `dynamic/` | Reservoir water and temperature, two melting ice controls, time-dependent route delivery, prescribed clearing, downstream heat and meltwater | Spatial routing, mobile grains, cliff collapse or field accuracy |
| `predictor/` | Minimum pulse duration and total water for first contact; deterministic parameter-corner sensitivities | Success of every hydrograph with those budgets, sustained recovery, statistical confidence |
| `spatial/` | Necessary latent-energy screen using terrain/head/heat bounds and a corrected finite-width relaxation | Positive discharge or recovery when exclusion fails |
| `audit/` | Evidence checklist and paired forecast scoring with abstentions and censored timing | Independent scientific approval or publication readiness |

The dynamic model adds a useful geometry-aware point-forecast baseline to the
previous static tools. It carries finite water and thermal storage, so it can
forecast a supplied hydrograph rather than merely test peak feasibility. Its
water, sensible-plus-latent energy and ice-volume accounting are explicit.
An implicit integrator handles fast storage adjustment and slow ice evolution.

There are **45 automated tests** across the four components, including analytic
drainage and reservoir-mixing limits, event continuity, budget consistency,
strict input handling, narrow-channel counterexamples, uncertainty semantics,
and censored scoring. Independent reviews of the dynamic physics are included.
Budget closure verifies algebraic consistency; observable convergence is a
separate requirement. The recorded synthetic examples do not validate glaciers.

## Run

    python3 verify.py
    python3 -m dynamic.cli dynamic/sample_input.json --output /tmp/recovery_dynamic.json
    python3 predictor/forecast.py predictor/sample_input.json --output /tmp/recovery_static.json
    python3 spatial/spatial_forecast.py spatial/example_input.json --output /tmp/recovery_spatial.json
    python3 -m audit readiness audit/example_readiness.json
    python3 -m audit evaluate audit/example_evaluation.json

Dependencies: Python 3, NumPy, SciPy; Matplotlib for figure regeneration. Consult
the component READMEs for schemas. They deliberately use distinct output labels
for first contact, energy exclusion and sustained operational recovery. Do not
score one target as another.

## Reproducible figures and synthetic checks

    MPLCONFIGDIR=/private/tmp/ice_novelty_mpl python3 -m dynamic.demo

The dynamic demo starts from an already abandoned toe and prescribes no clearing,
day-1 clearing, or day-10 clearing. It illustrates this reduced model's substrate
memory; it does not predict sediment mobilisation, discover the mechanism or
demonstrate environmental improvement. Source inputs and complete raw outputs
are saved in `dynamic/records/`. The spatial example includes a labelled PNG/PDF.

## Scientific focus

The testable contribution remains a forecast of the gap between debris export
and restored ordinary-flow stream heat delivery after differential ice incision.
This package supplies stronger baselines, useful bounds and scoring tools. It
does not resolve the current spatial solver's failed hydraulic partition gate.
No extra neural network, fracture module or debris law is added simply to
increase complexity before that gate is passed.

Read `CLAUDE_HANDOFF.md` for the bounded integration and finalisation sequence.
