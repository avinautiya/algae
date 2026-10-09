# Molecular-contribution protocol (follow-up study; frozen 2026-10-09 before any comparison below was run)

This is a follow-up protocol written **after** negative results:

- H1, H2, H4: NOT SUPPORTED;
- the physics forward model is too bright in the visible;
- TD-DFT D ≈ measured MAC on H1.

It is not a retrospective preregistration. Results go to `records/molecular_contribution/` (new directory). Nothing historical is overwritten.

## Questions

| Q | Question | Evidence that would count |
|---|---|---|
| Q1 | Do raw quantum-chemical predictions agree with **independent** laboratory spectra of the pigment? | Raw (uncalibrated) band positions and shape vs measured spectra that the TD-DFT inputs (molecule, functional, basis, solvent, FWHM) were not tuned to. |
| Q2 | With biology, ice, dust, illumination and the observation operator held identical, does substituting the molecular absorption treatment change and improve glacier-relevant optics? | Matched pigment-level substitution (Experiment A). |
| Q3 | Does molecular information improve prediction under new biological/environmental conditions where empirical optics would need recalibration? | Untouched observations from a new campaign/site-year with independent biology. |

Success on one question does not establish the others.

## Facts established from the current code before this protocol (inspection, not experiments)

- **Calibration (`phase2/tddft_calibration.py`):**
  - shape parameters (dE, w) are fitted to unit-area HPLC diode-array spectra of the isolated uncomplexed phenolics (Williamson et al. 2020 deposit, 265–600 nm);
  - magnitude ln f and Fe fraction φ are fitted to the whole-extract MAC (Williamson 2020, 260–750 nm, per kg phenol equivalents), with the measured Fe increment D (Procházková 2025) carrying the complexed absorption.
- **Downstream only the calibrated shape reaches the glacier model.** The TD-DFT absolute oscillator strengths and the molecular mass never do: f absorbs the glucoside → phenol-equivalent conversion and the f-error.
- **`measured_mac_C`** uses the Williamson phenolic extract MAC directly, the same data that set tier D's magnitude.

## Data roles

| Data | Previous use | Role here |
|---|---|---|
| Williamson 2020 HPLC isolated-peak spectra (peaks 2–4) | calibration shape target | Q1 raw-shape comparison: **not independent** of model selection (the same data informed earlier calibration choices). Reported as "raw vs calibration data", not validation. |
| Williamson 2020 extract MAC (53 samples) | calibration magnitude target; `measured_mac_C` | Q1 raw visible-absorption comparison (shape only; concentration basis is phenol equivalents) |
| Procházková 2025 Fe increment | tier D | not used for Q1 (measured increment kept distinct from computed complexes) |
| S6 2017 albedo (H1), S6/S Greenland HCRF (`heldout_v2`), KAN stations (H2/H4), PROMBIO 2021/23/24 (H5) | all examined | **development diagnostics only** for Q2 |
| UPE_U 2018 counts + UAS | not acquired, not examined | **candidate Q3 test set**, if acquisition, pairing and an observation operator are verified BEFORE examining targets |
| Any other untouched campaign | — | Q3 requires one; otherwise Q3 = **underdetermined** |

## Fixed candidate treatments (no additions after results)

| ID | Absorption treatment | Calibration |
|---|---|---|
| T-B3 | TD-B3LYP, 6-31G*, PCM water, 30 roots (`phase1/results/level2`, provisional geometry) | identical procedure, AR(1) residuals |
| T-CAM | TD-CAM-B3LYP (TDA, 15 roots), same molecule/geometry family (`phase1/results/legacy_2026-10-09/level2_cam`) | identical procedure, AR(1) |
| T-MEAS | measured Williamson extract MAC (no TD-DFT) | none |
| T-A | BioSNICAR Tier A empirical cell optics | none: Experiment B only (different biological model) |

- Each TD treatment is used as **tier C** (uncomplexed) and **tier D** (with fitted φ × measured Fe increment).
- Equal calibration data and equal parameter freedom (dE, w, ln f, φ, two discrepancy SDs, ρ grid) for T-B3 and T-CAM.
- No new parameters.

## Experiments and endpoints

1. **Raw vs calibrated (Q1, §5 of the request):**
   - per treatment: raw band maxima (eV) vs HPLC maxima; raw visible (400–700 nm) MAC fraction;
   - calibrated parameters (posterior mean/SD, correlations);
   - calibrated MAC difference between T-B3 and T-CAM.
   - **Endpoint:** max |Δ calibrated MAC| / MAC and the solar-weighted (300–750 nm) relative difference.
2. **Experiment A (Q2), pigment-level substitution:**
   - identical CellModel (species geometries, phenol-equivalent intracellular concentration, Williamson photosynthetic pigments, packaging, Mie g), identical BioSNICAR column, dust, illumination (SZA 50°, clear sky), observation quantity (hemispherical spectral albedo; S2 SRF bands; broadband 300–2500 nm);
   - states: log B ∈ {3, 4, 4.5, 5}, r ∈ {1000, 3000, 8000} µm, dust ∈ {3×10⁴, 3×10⁵} ppb, plus zero algae;
   - **Endpoints:** Δ broadband albedo, Δ S2 bands, Δ absorbed SW (W m⁻²) between treatments, compared with the observation-uncertainty floor (albedo 0.01; HCRF band noise 0.01–0.02; absorbed SW 10 W m⁻²).
3. **Experiment B:** complete approaches (T-A vs T-B3-D) are already scored in H1/H2/H4 and are not re-labelled as molecular tests.
4. **Development diagnostics:** Experiment A treatments on the 47 S6 plots at the reference nuisance state, reusing `phase4/forward_diagnostics.py` outputs where they exist.

## Acceptance criteria and decision gates

- **A (molecular predictions independently supported):** requires independent laboratory spectra. If none exist, A = **unvalidated**.
- **B (molecular treatment changes glacier optics meaningfully):** |Δ broadband albedo| between T-B3 and T-CAM (or T-B3 and T-MEAS) ≥ 0.01, or |Δ absorbed SW| ≥ 10 W m⁻², at any physically plausible state of Experiment A. Otherwise "molecular choice not consequential after calibration".
- **C (substitution improves independent observations):** only on untouched data, contrast ≥ 0.01 albedo with a block interval excluding 0. With no untouched data, C = **underdetermined**.
- **D (transfer without recalibration):** a new campaign with no refit. Otherwise **untested**.
- **E (melt):** requires H3; **untested**.

## Uncertainty

- Calibration posterior draws are propagated jointly where used (same draw through MAC → cell → albedo).
- Differences between treatments are reported against:
  - (i) the calibration posterior spread;
  - (ii) the numerical error of direct BioSNICAR (exact; no emulator);
  - (iii) the observation floor.

## Stopping conditions

- If B fails (molecular differences vanish after calibration), no expensive chemistry sensitivities are scheduled for prediction purposes. Remaining chemistry is justified only for mechanistic questions, which are stated.
- If B passes but no untouched data exist, the study stops at "underdetermined" and specifies the measurements needed.

## Compute

All runs go through `common/run_budgeted.py`. Calibrations are cached on disk by content fingerprint (`phase4/results/cache/calibrations`). Running chemistry jobs are not interrupted.
