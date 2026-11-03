# Continue without overwriting frozen experiments

This package is additive. Do not modify another worktree, background job, original cliff solver, prior frozen evaluations or the failed Swiss verdict. The original Arolla pre-execution freeze is commit `bc7916b852b9667820ef68e7442d4759e06ff342`. Read `RESULTS_REPORT.md` first.

## Implemented and independently checked

- True Arolla ice-ablation comparison, all 1,053 rows, independently recomputed from native sources.
- Exact conditional shared-reading and visit-time uncertainty audits; no robust pairwise field-superiority result.
- Executable source-gated transient debris conduction in `model/thermal_memory/`, with basal melt component, stored heat, conservation ledger, finite-scenario envelopes and 14 tests.
- Native input availability audit, environmental process comparison, Pirámide 2022–23 archive metadata and selective source extraction.
- Reproducible field-comparison and illustrative thermal-memory figures.

Run `python3 verify_ablation.py`. The current full suite has 74 checks. Passing checks certify implementation under tested assumptions, not exceptional field accuracy.

## Next executable step: source binding, then a genuinely separate test

Read `model/environment/PIRAMIDE_NEW_DATA_GATE.md` and `PIRAMIDE_PROCESS_BINDING.json`. The new source campaign is 2022–23, distinct from the earlier 2014–15 forcing archive. Native thermistor/moisture/weather variables are opaque MATLAB timetable objects. Python loading has NOT established trusted observation times or columns. Do not guess their structure or flatten object indices as time.

`model/environment/export_piramide_process.m` is a prepared export script, not MATLAB-tested. In an environment with MATLAB, execute it against the fingerprinted native process MAT file; independently verify row times, timezone, variable names/units, sensor depths and random source rows. It deliberately does not load camera-ablation targets. If MATLAB is unavailable, obtain equivalent CSV/NetCDF from the public source or an audited decoder; do not substitute invented observations.

Before any new melt scoring, freeze source hashes, training/evaluation dates, sites, physical model assumptions and metrics. Fit thermal parameters only to training interior temperatures, with physically justified bulk heat capacity; rock density times rock specific heat is not bulk debris capacity without porosity/moisture. Exclude the source's ablation-optimized conductivity and roughness from an independent test. Do not infer ice-contact geometry from the deepest sensor alone, especially Tower 2's 44 cm sensor versus 45 cm layer.

The proposed process test is stated in the new-data gate: compare transient conduction against appropriate source thermal-profile methods and quasi-steady conduction on later, source-matched observations. Camera uncertainty, snow state and basal temperate-ice assumptions need explicit evidence. The present new solver does not handle negative-temperature/refreezing or pore-water advection. A missing input is a model-scope decision, not permission to fill zero or silently clip temperatures.

Adding heat conduction/storage is established physics, not a novelty claim. A positive later comparison would support that specified test only. Do not call it validation of the original stream/cliff model or use inspected Arolla evaluation targets to choose weights, parameters, clocks or a new favorable primary metric.
