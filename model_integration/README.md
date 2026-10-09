# Algae → glacier-model integration

This additive package turns the project's optics into an **offline, testable model-coupling component**. It does not establish improved melt predictions, implement a prognostic bloom model, or modify E3SM/MAR. Existing production scripts are unchanged so concurrent Phase 1 work is unaffected.

## Research basis and defensible gap

Glacier algae are not absent from all existing models. [Cook et al. (2020)](https://tc.copernicus.org/articles/14/309/2020/) already combined algae optics, field observations, remote sensing and melt attribution. [Williamson et al. (2020)](https://pmc.ncbi.nlm.nih.gov/articles/PMC7084142/) connected pigment photophysiology with darkening. [SNICAR-ADv4 (2022)](https://tc.copernicus.org/articles/16/1197/2022/) includes glacier algae, bubbly ice and layered impurity distributions; its section 3.3.2 reports discrepancies at high algal concentrations and uncertainty in algae optical properties. This supports testing an optical improvement, not assuming one.

[Whicker-Clarke et al. (2024)](https://eesm.science.energy.gov/publications/effect-physically-based-ice-radiative-processes-greenland-ice-sheet-albedo-and-surface) integrated physically based ice radiation and satellite-derived properties into E3SM. [Williamson and Tedstone (2026)](https://www.nature.com/articles/s43247-026-03758-8) modeled Greenland-wide blooms and documented major uncertainties in loss processes and spatial ground truth. This project does not solve those ecological limitations merely by improving absorption.

**Candidate contribution:** test whether chemically informed pigment spectra, cellular packaging and explicitly represented subpixel heterogeneity improve transferable albedo/absorbed-SW predictions over established algae optics, while propagating calibration uncertainty. Novelty and superiority are hypotheses until demonstrated. A molecular calculation alone is not a climate-model improvement.

## Modules and exact conventions

- `coupling.py`: paired spectral budgets; guarded host-albedo replacement/anomaly modes; area mixing; interval energy integration; explicitly conditional offline potential-melt diagnostics.
- `biosnicar_adapter.py`: evaluates existing Phase 4 optics with algae present/absent, keeping ice and dust identical. Supports `ours` and `tierA`; a separate adapter is required for each direct/diffuse configuration and process. Uses Phase 4's private `_Builder`, so real-backend regression testing is required when that interface changes.
- `export.py`: immutable NPZ spectra plus versioned JSON metadata/checksum. Records source/config/data provenance and calibration draw IDs; no pickle, silent extrapolation or universal-LUT claim.
- `benchmark.py`: four-model held-out scoring, sample/group overlap rejection, paired group-bootstrap RMSE differences and provenance hashes. This scores predictions; it does not fit models or prove a submitted split is independent.

Spectral arrays are `(..., wavelength)` and flux is **bin-integrated W m−2**, not a spectral density. If host flux is W m−2 nm−1, integrate using the host's actual bin boundaries first. Albedo and flux must have identical shape: broadcasting across unknown member/time axes is prohibited. Downwelling SW and non-SW energy terms are positive toward the surface. Melt is kg m−2, numerically mm water equivalent; latent heat is 334,000 J kg−1.

The adapter retains the existing 300–2500 nm diagnostic band and allocates supplied broadband SW to it. That is a labeled band-limited approximation, not a full-spectrum/cloud correction. For rigorous host coupling use the host's resolved spectral flux and corresponding albedos in `spectral_budget`. Direct/diffuse absorption must be evaluated separately and summed using their respective bin fluxes. At night skip RT and use zero SW; do not evaluate a daytime solver beyond its valid zenith range. The underlying bridge rounds SZA to degrees; quantify this numerical approximation before production use.

Cell abundance is **cells per mL meltwater**, not biomass, cells per mL solid ice, or cell column density. Existing bridge handles its meltwater conversion. Never convert GA_BLOOM ng dry weight/mL to cells without measured, uncertain dry mass per cell. With film-only algae, mineral dust stays uniform over the 2 cm crust (`phase2/biosnicar_bridge._layer_concs`); other dust profiles are not represented.

## Host-coupling modes (`host_modes.py`, added on `claude/sweet-noether-6pryo4`)

- **Every host must declare** `HostContract(host, algae_treatment ∈ {absent, explicit, implicit_observed}, albedo_source ∈ {prognostic, observed})`. Nothing is inferred.
- **Mode A (forward parameterisation):**
  - The host computes its own albedo.
  - The optics add the algal change for a **prescribed** abundance. There is no bloom-growth model in this repository.
  - An `absent` host allows `anomaly` or `replace`.
  - An `explicit` host allows only `replace`, and only with `own_algae_scheme_disabled=True`.
  - An observed-albedo host is refused.
- **Mode B (satellite-informed state estimation):**
  - The observed albedo is used unchanged.
  - The optics attribute part of the darkening: the counterfactual no-algae albedo is observed + Δα, signed, and flagged when outside [0, 1].
  - Every output carries `validation_status` = "not an independent validation".
- **Named host for real tests:** the reference point SEB `phase4/seb.py` (KAN_M forcing), through `reference_seb_paired`.
  - Output labels keep these quantities distinct: `modelled_increment_mwe`, `potential_melt_mwe`, observed ablation (never produced here) and runoff (not modelled).
  - No E3SM/MAR integration exists.
- **Provenance:** reconciled from `codex/glacier-model-integration` at 0fc0096. Its `.github/workflows/model-integration.yml` was deliberately NOT brought over.

## Minimal coupling example

```python
from model_integration.biosnicar_adapter import AlgaeOpticsAdapter
from model_integration.coupling import couple_host_albedo
# cfg is an existing phase4.emulator.EmulatorConfig; use validated Phase 1 inputs.
adapter = AlgaeOpticsAdapter(cfg, phase1_l2="phase1/results/level2")
out = adapter.evaluate(cells_ml_meltwater=5000, f_n=0.5, r_um=5000,
                       dust_ppb=300000, sza_deg=45, sw_down_w_m2=500)
# FULL replacement only when modeled ice/dust/state/illumination match the host.
new_albedo = couple_host_albedo(host_albedo, out["albedo_without_algae"],
                                out["albedo_with_algae"], mode="replace",
                                host_algae_treatment="implicit_observed")
```

Anomaly mode requires an explicitly algae-free host. Observed satellite albedo already contains biological darkening; adding a second algae decrement is forbidden. Replacing the full albedo scheme is different from adding a correction and must be tested against the host's original scheme. Dust, crust, snow and pond changes must not be attributed to pigments.

Area-average solved patch albedos/absorbed energies, not cell concentrations before nonlinear RT. Carry joint calibration draws through every patch/time/model evaluation, preserving common draw IDs across pixels. Independent per-pixel error bars do not represent regional uncertainty when optical calibration is shared.

## Independent experiment required for the paper

1. **Freeze hypotheses and splits before fitting.** Primary endpoint: held-out absorbed-SW RMSE, provided independent measured downwelling/upwelling SW exists. Secondary: spectral/broadband albedo and abundance; actual melt only with suitable independent measurements and a validated host SEB. Define site-year/spatial blocks and inventory every observation used in spectra fitting, hyperparameter choice, discrepancy/tau fitting and masks. Hold out entire sites/years where possible. An outer test group may not influence inner calibration, tau, prior selection or feature processing.
2. **Run four fair baselines:** (a) identical ice/dust with no algae; (b) existing empirical algae optics (`tierA`, with documented version); (c) empirical satellite-albedo approach trained on training data only; (d) pigment/cell model. Add an independently configured SNICAR-ADv4 or host-native algae baseline if claiming improvement over that model: Tier A in BioSNICAR is not automatically a reproduction of SNICAR-ADv4/E3SM. Use identical forcing, observation support and exclusion rules for all models. Compare full pipeline predictions as well as controlled optics-only experiments.
3. **Avoid circular validation.** A satellite pixel used to retrieve albedo cannot independently validate that same albedo. AWS upwelling SW used to prescribe albedo cannot also be its test target. Use ground observations not used for inversion/calibration, match times and footprints, and audit overlap with S6 training data. Date mismatch, BRDF versus hemispherical albedo and spatial heterogeneity need explicit uncertainty/sensitivity tests.
4. **Resolve real energy forcing.** Use timestamped, quality-controlled interval means and actual durations, not hard-coded 15 UTC or clear-sky RF multiplied by measured-overpass ratios. Integrate `(alpha_without(t,lambda)-alpha_with(t,lambda))*F_down(t,lambda)` through time. Test hourly SZA, direct/diffuse weighting and spectral cloud effects. Reject missing records; do not silently interpolate across gaps. Six irregular cloud-free snapshots are not a seasonal integral.
5. **Validate melt with a host SEB.** Hold non-SW meteorology identical in paired runs; solve surface temperature, cold content, penetrative absorption, turbulent and conductive fluxes, refreezing/retention as appropriate. Compare to observed ablation with measurement uncertainty. The included reservoir diagnostic is only potential melt: positive column absorption does not automatically become surface melt or runoff.
6. **Report failures and uncertainty.** Report all four scores, number of independent groups, paired RMSE differences and confidence intervals. Use posterior predictive intervals from complete joint calibration/model/forcing draws; the sampling bootstrap is not a replacement. Test pigment-only, packaging-only, Fe increment and subpixel ablations. If sophisticated optics do not outperform baselines, report it; never tune the held-out data to make an impressive result.

### Data route

Start with [GEUS/PROMICE curated AWS data](https://promice.org/download-data/) and the [GEUS THREDDS catalog](https://thredds.geus.dk/thredds/catalog/catalog.html). Freeze the dataset version/DOI, download checksum, station coordinates and QC flags. Read the [official AWS metadata](https://thredds.geus.dk/thredds/fileServer/aws/metadata/AWS_data_readme.pdf) before choosing radiation/height variables. Station availability alone does not establish co-located independent algae/pigment observations. Surface-height change requires snow/ice classification, density conversion, instrument/tilt QC and uncertainty; do not label all height change as algal ice melt. The 2026 GA_BLOOM paper links data/code at Zenodo record 20138073, useful for a documented future biomass interface rather than unvalidated unit conversion.

## Benchmark command and input contract

CSV columns: `sample_id,group_id,observed,no_algae,established_algae,empirical_satellite,pigment_cell`. One row per uniquely identified held-out observation; units and supports identical for all columns. No fabricated demo observations are supplied.

Manifest JSON must contain `quantity`, `units`, `observation_source`, `training_ids`, `training_groups`, `model_provenance` (all four model keys, with code/data/config and calibration provenance), and `calibration_frozen_before_test: true`. Train/test groups must be disjoint. Include all data-informed choices in the training inventory, not merely optimization samples.

```sh
python -m model_integration.benchmark predictions.csv manifest.json --output benchmark_report.json
python -m unittest discover -s model_integration/tests -v
```

Run separate reports for albedo, absorbed SW and independently measured melt. Scores cannot be converted from one endpoint to another by renaming the quantity. Bootstrap intervals resample entire groups; fewer than ten groups produces an explicit warning, and a single group produces no confidence interval.

## Merge and readiness gates

Starting commit: `e87a6d14124e1653072a0d7c61a1d593fa16a3fc`. This branch adds only `model_integration/` and a dedicated CI workflow; it does not alter Claude's active jobs, production tolerances or existing result files. Review/cherry-pick after reconciling concurrent changes. Existing `phase4/multi_scene.py` remains legacy and its overpass-ratio/seasonal-summary outputs are **not endorsed** by this component.

Required before paper claims: real BioSNICAR parity test against Phase 4 RF on validated inputs; wavelength/dust/illumination checks; numerical sensitivity at large SZA; independent benchmark outputs with audited source manifests; complete uncertainty propagation; host-specific SEB validation. Also resolve existing calibration, inversion and satellite QC issues identified in the project audit. The analytic tests and mock adapter tests are software checks, not scientific evidence of improved predictions. Do not say this project is perfect, award-guaranteed, first-ever algae modeling, or integrated into E3SM/MAR on this branch.

## Coupling interface v1 (`model_integration/host_modes.py`, interface version "algae-coupling/1")

| Field | Units / range | Notes |
|---|---|---|
| abundance | log10 cells mL⁻¹ meltwater, 1–6 | prescribed (Mode A) or retrieved (Mode B); no bloom-growth model |
| f_n | 0–1 | fixed at the prior mean (not identifiable) |
| r_um (bubble optical radius) | 300–20 000 µm | host or retrieval state |
| dust | 3×10⁴–1.2×10⁶ ppb, uniform in the 2 cm crust | host or retrieval state |
| SZA | frozen emulator nodes 44–56° (±1°); otherwise refused | no extrapolation |
| spectral support | 300–2500 nm, BioSNICAR sub-Arctic summer clear-sky, direct beam | broadband albedo / Δα |
| outputs | albedo (1), Δα (1, signed), modelled increment (m w.e.), potential melt (m w.e.) | observed ablation is never produced; runoff is not modelled |
| flags | `validation_status`, out-of-range counterfactual flag, host contract | Mode B always "not an independent validation" |
| provenance | emulator file sha256 (in scene/station records), calibration draw IDs (`posterior_mean` only, D1) | joint draws pending |

**Runnable example:** `python3 -m model_integration.example_reference_seb` (KAN_L 2022; outputs in `records/host_coupling/`). Runtime 4 s.

Numerical checks recorded in `numerics.json`:
- n_sub 1 vs 4 seasonal melt differs by 0.16 %;
- algae-free albedo identical across optics to 1×10⁻⁵.

**What was changed in a host:** nothing in E3SM, MAR or any external host. The tested interface is the offline reference point SEB in this repository.
