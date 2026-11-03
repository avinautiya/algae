# A defensible route to a stronger field predictor

## What can and cannot be tested now

The inspected Arolla, Swiss and Changri Nup outcomes are development evidence. They cannot become fresh confirmatory tests by changing a label, metric, density, dates or a parameterization. The Arolla uncertainty audit found no robust pairwise dominance. This does not prove all models equivalent; it shows those bounds do not identify their ranking.

**A prospective computational comparison remains possible** using genuinely unused native observations after a single declared model-selection stage. Public accessibility does not itself invalidate an external test; inspecting or using its outcomes during design does. Keep a dated ledger of target files/rows, previous scoring attempts, all tried candidates, parameter dependencies and model decisions. If a new target is opened and used to improve the model, move it permanently into development and acquire another final test.

## Independent source audit

### Pirámide's richer 2022 campaign

The [2024 public dataset](https://zenodo.org/records/11581073) advertises camera-based ablation at three towers, debris temperature/moisture profiles, atmospheric profiles and lithological density/specific-heat inputs. It also supplies CRh/CRi/NYB conductivity estimators and an optimized-k/z0 calibration routine. Thus temperature-inferred material properties are already a published baseline, not an unexamined idea. The parent data agent is auditing the actual code/metadata and camera parser; no camera targets were loaded by this review.

The final [JGR2025 paper](https://doi.org/10.1029/2025JF008360), rather than the earlier preprint, supplies the authoritative method. Its §2.3.1.4 fits k to camera melt and checks deeper temperatures; that is not held-out melt validation. Camera lowering is converted with917kg/m³ (§2.2.3). Near-surface temperatures reach below0°C (§2.2.2), so the current nonnegative-temperature/no-refreezing component cannot process the whole campaign by silently clipping. The top probe also lies under a thin gravel cover; a boundary at that sensor requires the matching reduced depth. Full T&C reproduction additionally needs meteorological gap filling and topographic corrections (§2.3.4), beyond the decoded local profile files.

Eligibility is promising but conditional. Need exact sensor depths, timestamps, geometry, bulk rather than rock heat capacity, initial-state/spinup support, native ice-ablation conversion and uncertainty, and the paper's calibration target use. Reproducing an optimized model against the same camera data used to optimize k is not independent validation. Its target-fitted parameters must be re-estimated on development towers/windows only. Prefer withholding an entire tower plus a contiguous later event; overlapping frames from the same stake are not independent samples. Three towers at one glacier cannot by themselves establish cross-glacier generalization.

### Kennicott: useful data with a crucial leakage trap

The [2025 original paper §§8.2–8.3](https://rmets.onlinelibrary.wiley.com/doi/10.1002/gdj3.70032) identifies four colocated thermal-profile/stake pairs: KDT2020b/P_AWS_1 and KDT2023a/P11,2023b/P12,2023c/P13. It excludes initial installation disturbance. Its reported bulk heat capacity is pooled from melt-constrained conductivities at all four pairs, and this pooled capacity converts diffusivities into the other reported conductivities. **Using those published C/k values to validate the same 2023 melt would leak the target.**

A valid split can instead estimate diffusivity using temperature-only development data, constrain bulk heat capacity with independent lithology or the 2020 development melt only, and freeze before opening 2023 melt. The observation/geometry and site transfer uncertainties remain. Published table values are not automatically independent because they have physical units.

Version check: the paper cites [Zenodo14625691](https://doi.org/10.5281/zenodo.14625691); the newer [v4 record16389806](https://zenodo.org/records/16389806) has a190.9MB archive. Pin one version and inspect changes rather than mixing files. The model agent is auditing the smaller historical version; this review did not duplicate that download or inspect its melt targets.

### Other sources

[Satopanth2022](https://doi.org/10.1017/jog.2022.35) already studies estimating ablation from debris-temperature profiles, so a new profile-to-melt method must compare with that established formulation. Public machine-readable independent inputs/targets were not established here. [Thermal-imagery2021](https://www.frontiersin.org/journals/earth-science/articles/10.3389/feart.2021.681059/full) documents camera stake measurements and surface-temperature methods, but its availability statement offers raw data through the author; an accessible complete archive was not identified here. Neither is marked immediately executable.

## Candidate and information setting

A physically admissible thermal-memory model may improve some data, but that is a hypothesis. The existing diffusion solver is verified established physics. Develop it against observations without calling the PDE new. Distinguish:

- **Conditional melt estimation/nowcasting:** measured surface or within-debris temperature is available during the target interval. Strong same-input baselines include steady conductive flux and the published temperature-profile estimator.
- **Forecasting:** future surface temperatures are unavailable unless independently predicted. Then compare a complete meteorology-to-surface-temperature-to-melt system against full energy-balance and calibrated temperature-index baselines, with the same information cutoff.

Giving a candidate measured future surface temperature while its comparator receives air temperature alone cannot establish that its model dynamics are superior. Report such different-information comparisons separately.

## Freeze and model-selection protocol

1. Declare one deployment question, target quantity, identical temporal/spatial footprint and fixed practical improvement margin before final scoring.
2. Assign inspected targets permanently to development. Group shared camera/stake readings, storms and site-years. Use glacier-level blocks for a cross-glacier claim; finer blocks support only narrower within-site statements and need a defensible dependence model.
3. List **all** attempted hypotheses and a fixed development-only selection rule. Nested blocked validation can select conductivity parameterizations, memory/state structures and regularization inside development. No selection uses final melt targets. Temperature-based calibration at evaluation locations must be declared an assimilation input and provided to same-information competitors; it is not target-blind parameter transfer.
4. Record transitive dependencies of every parameter. A published constant pooled using any evaluation melt target is prohibited, even if the local code never reads that target.
5. Pin code/data versions and hashes, observation operator, spinup, snow handling, sensor uncertainty and clock scenarios. Predeclare failure/abstention behavior; do not drop difficult targets after seeing their errors.
6. Freeze the selected candidate, all baselines, predictions/algorithm, metric and tests. Score exactly once; keep every prespecified result. A failed primary test remains failed even if a secondary metric improves.

## Executable gate and conservative statistics

`comparison_validator.py` checks the declared protocol and refuses source/target leakage, missing hashes, mismatched support/components/information, incomplete parameter provenance, absent block definitions and overwriting a freeze. It validates declarations; it cannot authenticate a false history or make missing data real.

`conservative_block_test()` consumes one already-aggregated loss-difference interval per **independent block**, not individual image frames. Bounds must already honor joint clock and shared-reading dependence. A block is a certified win only when the upper bound of candidate-minus-baseline loss is negative. Every uncertain/tied block counts as a non-win. The one-sided exact binomial tail tests independent block-win probability against1/2. A frozen Holm family adjusts the reported baseline comparisons. This is deliberately conservative and presumes justified independence and valid external bounds; it is not a universal superiority theorem.

With three baselines, even5of5 certified independent wins fail this Holm family at0.05;6of6 pass. Two glaciers are insufficient for that test regardless of camera frame count. Do not manufacture extra independent samples by splitting continuous data. If the available independent units are too few, report effect sizes, conditional error bounds and replication needs without a confirmatory statistical superiority claim.

### Concrete pass conditions

- All source, component, geometry, time and provenance gates pass before targets are scored.
- Primary native ice-ablation loss improves by the declared practical margin on an unused target under declared measurement/clock bounds.
- Full predeclared baseline family is reported, including a well-reproduced professional baseline with verified ice quantity.
- Stronger multi-environment claims additionally pass an appropriately powered, independent blocked comparison; a single new event supports only that event.
- If these conditions fail, freeze/report failure. Neither extra computational cost nor an ex post uncertainty interpretation repairs it.

Seven executable tests pass, including independent binary enumeration of sign-test tails, target/published-parameter leakage refusals, missing component/input gates, block pseudoreplication and family adjustment. No new field target was scored here.
