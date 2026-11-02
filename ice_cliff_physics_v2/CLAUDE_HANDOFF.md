# Integration handoff — additive cliff physics and actual field tests

## Location and isolation

Branch: `codex/cliff-physics-v2-20261101`; package prefix: `ice_cliff_physics_v2/`.
The synthetic source/spec/case pre-execution commit is `b9d6f39d09f7ead306ea87c852da8f14c53e9e89`.
The branch inherits `67f9b6eee58ae19738cd5c1ca3d78d2a543ff010` and the earlier Swiss/external forecast evidence. It does not change `claude/ice-cliff-survival-m0`, its jobs, frozen outputs or checkpoints. No missing external packet is required: executable code, cases, observation inputs, tests and reports are in this prefix.

## Start here

Read `RESULTS.md`, `METHODS.md`, `reviews/FINAL_RESULT_AUDIT.md`, then the field benchmark reports. Keep temperature, seasonal ablation and coupled cliff claims separate. This package is a reduced common-stage model, not a replacement for the existing SWE solver or validated resolved CFD.

```sh
cd ice_cliff_physics_v2
python3 verify.py
python3 verify_field.py
python3 model.py example_temperate_profile.json result.json
python3 run_study.py --out results/coupled --freeze-commit b9d6f39d09f7ead306ea87c852da8f14c53e9e89
python3 run_boundary_study.py --out results/boundary --freeze-commit b9d6f39d09f7ead306ea87c852da8f14c53e9e89
python3 analyze_study.py
```

The case runner verifies input/source hashes, retains all failures/refusals, resumes only matching cases, and writes atomically. Do not edit the existing freeze to accept changed code. Create a versioned new study if changing production equations. The plotting and field code require the analysis dependencies; the production core is standard-library based.

## Valuable integration changes

1. Use finite inlet heat once, allocating floor heat before the shared bare/debris face operator. Preserve signed debris→water export and thermal storage in the ledger.
2. Use separate wet/dry columns and conservative donor transfers. Never reset wet material to0°C on reclassification or mix all covered material globally.
3. Keep true height-dependent ice volume separate from basal routing geometry. Do not validate actual face melt with a uniform-step routing volume.
4. Insert cover/water boundaries before sampling flux, retain actual adaptive band counts, and refuse cap exhaustion rather than silently coarsening.
5. Preserve local Reynolds/Froude/domain and per-route phasewater checks. An explicit free-drainage contract is a conditional approximation, not solved meltwater feedback.
6. If porting into SWE/CFD, use its local water temperature/discharge and preserve those ledgers independently. This package does not silently fix SWE momentum, mass-source or vertical-contact physics.

## Results to retain, including negative comparisons

-131 core tests;30 coupled cases and12 analytic-boundary cases complete.
-Initial early-time FV error11.78–16.41%; separate48-run refinement reduces it to0.331–0.469%. Do not present the one-diffusion-time accuracy as an early-time certificate.
-Field homogeneous diffusion improves interior temperatures over quasi-steady at11/11 profiles but loses to the ridge baseline. Conditional ablation transfer also loses to constant-rate transfer. Do not select weaker baselines to claim universal improvement.
-All4ablation origins are one glacier's profiles/years, not4independent glaciers.2023 leading forcing is missing; date/depth and reconstruction scenarios remain conditional.
-General storage, wet resistance and classical Robin physics are not new mathematical discoveries. A new environmental improvement still needs supported observation-level evidence.

## Next bounded scientific work

The vertical-heterogeneity test has now run prospectively:30.2%lower MAE than homogeneous conduction, but ridge remains better in9/9profiles. Its frozen superiority gate fails; retain that result. Do not reuse either inspected window as an untouched test. Do not multiply free coefficients or use the inspected test to manufacture a win. Coupled cliff accuracy still requires matched water/contact forcing and face retreat over the same physical support. Do not replace that missing test with more synthetic runs or label the component evidence as full-cliff validation.
