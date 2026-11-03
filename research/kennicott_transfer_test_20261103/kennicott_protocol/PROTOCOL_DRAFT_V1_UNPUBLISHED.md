# Conditional cross-year Kennicott net-ablation comparison: pre-outcome protocol draft

**State: draft before any native target magnitudes; not a published prediction freeze.** Source binding, candidate implementation and finite initial-gap scenarios must be finalized and hashed before TRAIN release, followed by prediction publication before2023 release. The existing evidence remains unchanged. Published source-paper aggregate context is known, so this is a transparently frozen retrospective external-year comparison, not a prospectively blinded field campaign.

## Declared target and information split

Native target is the **negative of reported signed mass change**, `Y=-mass_change_kg_m2`, directly in mm water equivalent:1kg/m²=1mmwe. Do not apply an additional900/917 density conversion. Source metadata says stake-verified changes between consecutive visits, without adjustment to inferred seasonal extrema. The intended quantity is net ice-ablation loss under the explicit temperate-debris/no-unmodelled-refreezing interpretation. Endpoint surface type is Debris and snow depth0. These observations support that conditional interpretation but do not demonstrate absence of all intervening snowfall, evaporation or storage. Predictions are compared to native net loss, never relabelled independently observed gross latent melt.

One global TRAIN target only: P_AWS_1,2020-07-25→2020-09-06,43source days. All local scale calibration for every candidate/comparator uses this same single TRAIN amount. No2023-specific gain, drift parameter, conductivity or target-driven window selection is allowed. Source2025 heat capacity/conductivity pooled using2023 targets is excluded. Independently source-defined2020-only parameter estimates may be included if their provenance and fit inputs are explicit.

Three evaluation targets: P11,P12,P13, each2023-06-08→2023-09-02,86source days. Native magnitudes remain sealed. Three stakes in one glacier/year are not three independent glaciers or weather experiments. No nominal95%confidence interval, bootstrap across three stakes or unconditional global-superiority claim is justified. Preserve all signed/absolute errors and sensitivity scenes regardless of outcome. A zero or negative TRAIN loss causes calibration refusal rather than absolute-value conversion; negative evaluation net loss remains signed and is not clipped.

## Clock and support

Visits have date labels only, without an independently verified hour/timezone. Profile records provide UTC plus GMT−08ADT columns. Nominal convention is **12:00 localADT**, hence20:00UTC, adopted before outcome access. Discrete date-uncertainty scenes use both UTC-date and ADT-date interpretations, with each start/end at0,12or24hours on the labelled day, for18shared-clock endpoint scenes per target. These are representative clock scenes, not a theorem that all intermediate-clock extrema are attained there for arbitrary signed/phase-dependent heat. Use an exact continuous clock enclosure if claiming complete time-box bounds; otherwise say finite-scenario sensitivity.

Keep full visit intervals. Never shorten the86-day target to the83-day observed temperature period. All2023profiles start2023-06-10T20:00UTC, so missing leading duration is36–68hours under the clock union and48hours nominally. This is a measurable3%period issue, not negligible by definition. For the2023air baseline, sourceAWSendsSeptember2at21/22UTC: late-end clock scenes may need explicit endpoint support/extrapolation flags; do not silently extrapolate beyond last native observation. Initial thermal state and adaptive initialization also need declared scenarios before2023scores. A noncausal future-window fit cannot be presented as a causal forecast.

A finite, physically stated initial-gap reconstruction can support a **conditional scenario comparison**; it is not an observed heat flux or a deterministic environmental bound. If outcomes are compared only after assigning that gap, state the assignment and report its effect on every model. Missing leading heat and geometry-transfer uncertainty may be soft conditional assumptions; an undefined target operator, mismatched target dates or unsupported below-ice sensor mapping is a hard refusal. Outcome values must never select a gap or clock scene.

## Geometry and footprint binding

| Stake | Stake debris | Profile debris | Native sensor depths |
|---|---:|---:|---|
| P_AWS_1 TRAIN |210mm|220mm|source2020profile metadata pending exact binding|
| P11 |240mm|240mm|20/70/120/170/220mm|
| P12 |290mm|320mm|90/165/240/315mm|
| P13 |255mm|260mm|20/95/170/245mm|

Thermistor accuracy is reported±0.2°C; do not interpret it as Gaussianσ. Source depth uncertainty is not established. The **profile** boundary must remain at its independently surveyed pit depth. In particular, the315mm P12-profile probe cannot be embedded in a290mm-thick stake column: that would put it below the declared ice surface. Derive basal flux at the320mm profile then apply an explicitly stated nearby-stake transfer, or reconstruct a separate290mm stake column from transferred material/boundary assumptions. Proximity alone does not make pit and stake the same footprint. Report the spatial binding and thickness-transfer assumptions. The same rule applies to2020's210/220mm difference.

## Candidate and strong comparators

The candidate is the pinned causal weak integrated positive-conductance thermal graph with its exact adaptation/phase/initialization specification supplied by the model agent. Future measured temperature profiles are contemporaneous explanatory inputs: call it a state-conditioned melt estimator, not an autonomous weather forecast. A continuous goal enclosure may be reported only for the verified fixed-coefficient/polytope/phase assumptions, with source-error/regressor uncertainty and failure/abstention flags.

Required comparators, all with the identical2020TRAIN information: source-supported native-air positive-degree-day index; constant TRAIN daily rate scaled by each whole target duration; a fixed-conductivity/temperature-gradient physical index with independently sourced or2020-only coefficients; and a nonadaptive graph using the same geometry/error/phase choices. Do not call a buried thermistor an air measurement. A source-adjusted thickness/resistance baseline is appropriate only with a declared mapping, and should be reported alongside the unadjusted comparator. Published parameter choices that used2023targets can be contextual only, never an untuned external comparator.

No comparator may use more2023target information than the candidate. A globally trained2020scale may be freely different by model, but must remain fixed across all three evaluation stakes. If a strong comparator lacks sufficient native forcing, preserve its unsupported status; do not invent weather or weaken it to create a win.

## Freeze and scoring

1. Finalize and hash feature-source files, source metadata, model/phase/gap/geometry scenes, all baseline code, operator and exact target keys. Record unresolved assumptions explicitly. No native2023magnitudes in this stage.
2. Publish normalized predictions and the protocol. Release **only** the2020TRAIN target under a hash-bound role receipt; derive one global scale per model and retain all calibration diagnostics.
3. Publish all calibrated2023predictions/scenario families and source-input hashes. Only then release the exact three2023target rows. No post-release model/scenario changes in this experiment.
4. Score in native mmwe using `campaign_scorer.py` after the source-bound target loader has applied the frozen sign convention. Primary report: all three point errors, unweighted MAE, mean signed error, all-scene sensitivity and whether the candidate beats each named baseline on every stake. Do not equate a pooled MAE win with consistent improvement. All three rows share one climate year; campaign-balanced MAE is the same as unweighted MAE here.
5. Native precision/measurement uncertainty remains unknown unless authoritative metadata provides it. If a δ sensitivity is shown, label it a conditional bound, not a source-known confidence interval. Frozen point ranking is reportable without inventing precision; robustness beyond the available uncertainty remains conditional.

This protocol permits a well-defined conditional environmental test without demanding impossible certainty about every omitted process. It does not guarantee a favorable result or establish novel physics by itself.
