# Kennicott paired thermal–melt eligibility

## Actual progress

A genuine external-glacier dataset with **three explicitly paired2023 thermistor/stake sites** is present. It is stronger than the one-stake Lirung option, but is not yet a fully supported frozen experiment. Neither numerical mass-change values nor new predictions have been inspected. The current requirement is a declared initial-forcing/geometry treatment, not a guarantee that every unknown parameter must be perfectly measured.

[Zenodo16389806 v4](https://zenodo.org/records/16389806) is190894018bytes. Only central-directory and allowlisted member byte ranges were fetched; the190.9MB archive was never downloaded. Individual files have verified ZIP CRC32 and SHA256 manifests. Whole-archive MD5 was not verified because the whole archive was not fetched. The raw outcome CSV was processed transiently for explicitly allowlisted site/date/geometry/snow fields only, not persisted; its numerical mass-change column was never converted or exported.

## Source interpretation

[Petersen et al.2025](https://doi.org/10.1002/gdj3.70032), Tables3–4 and methods, identify KDT2023a/b/c with P11/P12/P13 and2020b with P_AWS_1. The authors discard several days after thermistor installation because packing, moisture and thermal states are disturbed. Stake measurements use900kg/m³ ice density and seasonal snow accounting. Published conductivity values share heat capacity estimated from2020 and2023 melt-colocated profiles, so they leak2023 outcomes if reused as independent2023 parameters. Their thermistor diffusivities are temperature-derived; heat capacity may instead be fixed independently or constrained on2020 TRAIN data only. The source provides5-minute thermistors and local weather; the two are not proof of exact thermistor/stake co-location or constant layer geometry throughout the season.

## Verified native record binding

|Profile/stake|Native target dates|Target snow start/end|Thermal observations|Profile geometry|Stake debris thickness|
|---|---|---|---|---|---|
|KDT2023a/P11|June8→Sept2 2023,86days|0/0|27073complete5min rows|d24cm;z2,7,12,17,22cm|24cm|
|KDT2023b/P12|June8→Sept2 2023,86days|0/0|sameclockandcount|d32cm;z9,16.5,24,31.5cm|29cm|
|KDT2023c/P13|June8→Sept2 2023,86days|0/0|sameclockandcount|d26cm;z2,9.5,17,24.5cm|25.5cm|

Every profile begins **June10 20:00 UTC** and endsSeptember12 20:00 UTC. The seasonal target therefore begins about2–3days before the first thermistor record, before excluding installation disturbance. Do not change the target toJune10 without an observedJune10 stake increment. No mid-season2023 measurements occur at these three sites in the released point-balance table.

The target operator is signed `mass_change_kg_m2` over consecutive surveyed dates, numerically equal to mm water equivalent. Positive ice ablation would be minus this balance only under the explicitly snow-free and no-other-mass-change conditions. Native date strings are month/day/year, contrary to the metadata's advertised YYYY-mm-dd format. Visit times/timezone are absent for stake dates; source thermistors have both UTC and explicit ADT(UTC−8) labels. Preserve a predeclared visit-time uncertainty family instead of inventing exact UTC midnight.

P12's29cm stake thickness lies above its thermistor31.5cm depth; these cannot be imposed as a single exact column. Use separately justified profile32cm geometry plus a spatial representativeness uncertainty, or exclude P12 on a feature/geometry gate before outcomes. Selecting whichever thickness yields smaller melt error would leak outcomes.

`Debris_T2` and `Debris_T3` are hourly clocks throughSeptember2 21:00/22:00 UTC, starting2020, with respectively1and10missing values over their entire records. Exact2023-window missingness and spatial transfer still require checking before baseline integration. Air stations are not at each thermistor. Amazon Creek offers longer weather but is off-glacier; spatial transfer must be explicit and shared by all models.

## Legitimate pre-outcome experiment options

**OptionA: strict native thermal retrieval.** No-go for the entire seasonal target as currently released: the initial target days have no thermistor forcing, and the installation-reset period cannot simply be deleted from the observed seasonal total. Zero-filling or retrospectively rescaling the target by coverage fraction is invalid.

**OptionB: conditional whole-season prediction, three held-out stakes.** Potentially feasible. Freeze a source-justified weather-driven or uncertainty-bounded initial2–3day segment and installation equilibration treatment; then use all subsequent observed thermistors under a conditional-retrieval label. Freeze initial state, depth/geometry scenarios, signed heat/phase rules and all numerical tolerances. Use a transfer amplitude fixed from separate development data, or2020-only calibration, not the2023 published pooled heat capacity. Declare common native-air/calendar/thermal-dose comparators, using identical forcing availability and the same allowed calibration information. Preserve all three target outcomes unseen until all seasonal predictions and source/error envelopes are committed. Score three site-level seasonal errors; do not count5-minute records as independent melt targets.

**OptionC: three-fold leave-one-stake-out.** Mathematically possible after all operators are fixed, but each fold trains on two same-season stakes and evaluates one. This is small-sample spatial cross-validation, not independent climate generalization; global hyperparameters must already be fixed. Do not describe it as three prospective independent glacier validations.

A2020TRAIN→2023TEST split has real temporal separation and new stakes. However2020b has only two thermistors according to the published summary, and its July25→September6 stake interval also begins before the July27 thermal-profile record. Native sensor headers and coverage must be checked rather than assuming it supports the same inverse. If using the source's2020-only capacity value, label its melt-trained provenance explicitly and never use the pooled2023 value.

## Required release gate

Before requesting actual mass outcomes: freeze one option; complete feature-only source checks and all initial-period/geometry assumptions; produce full-period predictions for every candidate/baseline; report unsupported intervals rather than extrapolating silently; and predeclare numerical/missing-input uncertainty. If the early-gap envelope is too wide to discriminate models, report that outcome rather than shrink it after seeing melt. This source can support a useful external test, but the current audit does not establish superior field accuracy or novelty.

## Files and reproduction

- `ZIP_METADATA_ONLY.json`: every member's path, compression, sizes, CRC and byte offset.
- `METADATA_FILES_MANIFEST.json`: five verified source metadata texts.
- `PROFILE_SUMMARY_METADATA_ONLY.json`: source pairing/location/geometry, with derived conductivity removed from inputs.
- `PAIRED_STAKE_METADATA_ONLY.json`:12paired-site date/snow/geometry rows; mass values excluded.
- `FEATURE_CLOCK_AUDIT.json`: actual three thermistor and two weather clock/missingness checks.
- `list_remote_zip.py`, `fetch_metadata_only.py`, `audit_paired_metadata.py`, `audit_feature_clocks.py`: bounded-range reproducible audits. No external local modules outside this folder are required. Running requires approved public network access.
