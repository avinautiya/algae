# Independent prediction extensions and benchmark audit

This additive package was developed while Claude ran a different milestone. It changes none of Claude's source, jobs, checkpoints or frozen outputs. The four modules pass26 tests in total. There is no trained new field predictor or improved melt score in this release.

## Delivered code and executed checks

* `physics/`: energy-bounded viability gate including independently bounded external ice heat and route-specific water heat transfer. Six tests pass. The original21-day repair exclusion is reproduced. Extra ice heat could change that exclusion, but the calculated necessary threshold is not measured forcing or proof of recovery. Do not choose inputs just to defeat the bound.
* `novelty/`: an energy-conservative, serial upstream-cooling calculation. Five tests pass. The executed counterexample reverses local incision ranking without changing prescribed target flows. The demonstrated difference is only0.0438mm/day; it does not justify a consequential recovery claim. Upstream cooling is established physics, not a new general mechanism. Local melt redistribution is distinct from whole-glacier heat uptake.
* `benchmark/`: public Changri Nup data, a checksum-verified archive, loader, aligned audit, source metadata, eligibility gate and training-only evaluation protocol. Eight tests pass. The loader processed10,813hourly and450daily rows. Original measurement metadata are retained. SR50 lowering is not automatically pure ice melt. Heat-flux sign, sensor movement and exact evaluation windows remain unresolved. Shared sensor-derived inputs need an issue-time dependency audit; shared provenance alone is not proof of leakage.
* `scoring/`: a matched-support daily-melt score implementation with7 tests. It retains all frozen cases and models or refuses missing predictions; it cannot silently drop difficult cases. It reports absolute/signed errors and descriptive day wins, not significance or independent replication.

## A more relevant external comparison

DCG-MIP2026 compares15 debris-covered glacier melt models across9sites and links forcing, observations and submitted outputs: https://doi.org/10.5194/tc-20-1895-2026 . These submitted outputs have not been downloaded or scored in this package. The archived Changri Nup data supplies a concrete input/target audit, not automatic eligibility for an independent pure-melt benchmark.

This is a different observation problem from stream-induced local cliff retreat. A sub-debris column experiment cannot validate water-contact prediction or a full cliff model. Any candidate must be evaluated on the exact quantity predicted by its baseline, with matching support, units, snow treatment, issue time and calibration inputs.

## How this advances the project

The new functions can quantify whether independently measured extra heat or thermal paths could materially change the current no-go result. The field loader and eligibility manifest provide a route toward a professional-model comparison rather than repeatedly fitting the inspected Swiss intervals. They do not supply missing route temperatures or justify inventing them.

Before a field campaign: resolve the observation gates, obtain the submitted comparator outputs and their metadata, then freeze the common support, training/validation/test assignments, model choices and tuning budget. Audited datasets are not to be described as unseen. Keep all baselines, including ones that outperform the physical model. Historical originality, meaningful environmental impact, field superiority and competition outcomes remain unestablished.
