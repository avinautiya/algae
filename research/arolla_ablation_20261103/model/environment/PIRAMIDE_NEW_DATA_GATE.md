# Pirámide 2022–23 process-validation gate

**Decision:** promising genuinely different field campaign; **not yet executable native field prediction**, because native MATLAB timetable clocks/columns remain opaque in the available Python loader, camera uncertainty is unverified, and thermal-property independence must be enforced. No ablation member was downloaded; no target prediction or score was calculated.

## Verified sources and extraction

[Official record](https://doi.org/10.5281/zenodo.11581073), [API](https://zenodo.org/api/records/11581073), [archive preview](https://zenodo.org/records/11581073/preview/Debris_properties_estimations.zip?include_deleted=0), [associated preprint](https://doi.org/10.22541/au.172168451.12963281/v1).

The record supplies a 62,677,571-byte ZIP, published MD5 `964e16d7cba6117565384c7ab9b81c69`. Selective HTTP ranges retrieved only its 64 KB tail/index, 45 KB source-code prefix, and 2.6 MB process-observation member. Extracted-member CRC32 values match the ZIP index. **The whole-archive MD5 was not verified**, because the whole archive was not downloaded. SHA256 fingerprints and metadata are in `PIRAMIDE_PROCESS_BINDING.json`; source-member offsets/sizes in `PIRAMIDE_ARCHIVE_INDEX.json`.

The official preview lists three small camera-ablation files (31.1/28.3/13.1 KB), `DATA_Piramide.mat` (2.6 MB), lithology properties, moisture, and calibration/meteorology. `Readme.txt`, actual thermal code and numeric lithology were read locally. The largest 58.7 MB calibration meteorology member was not downloaded. Native process arrays are MATLAB `MCOS timetable` objects; `scipy.io.loadmat` exposes opaque references, not validated time/column arrays.

## Paper-confirmed measurement facts (sections 3.2–3.3, table 1)

The campaign is November 2022–March 2023, distinct from DCGMIP 2014–15. Five thermistors and three moisture probes accompany each of three towers. Debris thicknesses are 50/45/24 cm. Thermistor depths are 0/12.5/25/37.5/50; 0/11/22/33/44; 0/6/12/18/24 cm. Moisture depths are 30/40/50; 19/29/39; 12/17.5/21.5 cm. Temperature sampling is 30 minutes, accuracy ±0.4°C; moisture sampling five minutes, accuracy ±3%. Cameras observe marked stakes hourly 07–19 and at 03, with manual corrections and density 917 kg/m³. Two towers fell, interrupting records. A camera-ablation error model is not established here. One printed moisture end-year precedes its start; native clocks must resolve it. Closest comparisons already include homogeneous/two-layer profile methods and ablation-optimized conductivity. Do not label these methods new. [Preprint](https://doi.org/10.22541/au.172168451.12963281/v1).

## Independently decoded parameters and leakage exclusion

`LithologyData.mat` contains numeric rock properties:

| Tower | Density kg/m³ | Specific heat J/kg/K |
|---|---:|---:|
| 1 | 2684.933333 | 796.866667 |
| 2 | 2657.313333 | 816.126667 |
| 3 | 2657.753333 | 815.631667 |

These are source lithology-weighted estimates, not calorimetry or bulk heat capacity. The supplied calibration code assumes porosity 0.44; independent local porosity measurement is not established. It loads stake ablation explicitly, searches conductivity against ablation, and masks snow using daily albedo >0.3. Thus supplied **optimized k/z0 and calibration outputs are prohibited as independent melt-validation inputs**. The profile-method path loads debris thermistors; it offers a target-independent fitting route if recomputed on declared training times only. Using profile estimates fitted over the complete evaluation period would still violate a prospective temporal prediction design.

Do not transfer source code's units blindly: its ablation-conversion lines must be reconciled with raw camera displacement, provided mm-water-equivalent fields and the stated 917 kg/m³ density before scoring. Ten-centimetre stake markings are a length reference, not demonstrated camera precision.

## Concrete next executable model specification

This is a proposed implementation contract, not a completed model or new physics claim:

1. Execute `export_piramide_process(source_mat, output_dir)` in MATLAB. It loads **only** `debtemp`, `sm`, `trh`, `wind`; exports twelve CSVs and source-clock/column metadata. It preserves labels and does not invent UTC. Validate source timezone, ordering, cadence, missing flags and units against the native timetable. A timezone-empty export requires an authoritative clock interpretation before cross-sensor binding; relative-time thermal tests can proceed only after confirming common clocks.
2. Use Tower1 as the primary, predeclared stable tower. Freeze a calendar training block and a later disjoint evaluation block **after feature-only coverage/snow checks and before reading camera outcomes**. Do not choose blocks by melt-model fit. Require a numerical spin-up/initial-state buffer and exclude that buffer from evaluation.
3. Solve a conservative 1-D finite-volume debris conduction column with source layer thickness. Impose observed near-surface temperature as the upper boundary, and initialize from the measured vertical profile. Basal temperate ice is a declared hypothesis checked against bottom temperature, not inferred from melt agreement. At Tower2 the lowest sensor lies 1 cm above the stated ice interface; do not silently make it the interface.
4. Fit at most one constant diffusivity from **training-only interior temperatures**, using fixed geometry and bounded sensor errors. Derive k from diffusivity and independently specified/bounded bulk heat capacity. Primary comparison: one-layer conduction versus the source homogeneous profile method; stronger comparison: the source two-layer profile method, recomputed training-only. Hold interior thermistors out for time-prediction assessment; do not evaluate a sensor used as an imposed boundary.
5. Moisture is observed conditioning information, not a measured advective heat flux. A moisture-dependent extension must use a declared energy-conserving storage/transport law. Substituting a time-varying C into C*T_t without accounting for transported water enthalpy can violate conservation; no such extension is authorized by this gate alone.
6. Freeze predicted basal flux/melt and interval operators before opening camera-ablation members. Resolve whether basal flux is signed heat into temperate ice, net refreezing or gross melt; require matched ice area and units. Include cold-ice/snow gates. Compare with training-only constant melt, PDD and established source thermal models on identical supported windows.
7. Camera displacement-to-melt conversion, manual-correction flags, temporal resolution, onset/reset rules and uncertainty must be audited before a field-accuracy claim. If uncertainty is absent, present nominal errors explicitly without invented confidence intervals and test defensible independently justified error bounds.

## Useful contribution and claim ceiling

This dataset can support process fitting without using ablation as the training target, then a distinct field-ablation evaluation, **if** native exports and phase/observation gates pass. This is stronger identification support than Arolla's single training interval. It is still a weather/surface-temperature-conditioned column model, not a cliff-survival or glacier-scale forecast.

The main methodological neighbours are the provided source methods, [Laha et al.](https://doi.org/10.1017/jog.2022.35), and full debris energy-balance models. No new mathematical or environmental novelty follows merely from adding moisture, layers or temperature assimilation. Any claim would need a new independently supported result beyond those baselines. No superiority is presumed.
