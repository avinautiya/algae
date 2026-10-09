# Glacier-model validation protocol (H1–H4)

**Frozen:** 2026-10-09, before any value of the new validation datasets listed below was examined. The only exception is confirming file structure (column names and row counts of `Albedo_master.csv`; first 3 rows printed while checking the format).
**Scope:** the evidence needed for claims that this project improves how glacier surface models represent algae-driven darkening.
**Changes:** any change after data are examined is logged in `docs/preregistration_deviations.md` with a date and reason. The primary endpoints and decision rules below are not changed after results are seen.

## 0. Principles

- The four hypotheses are separate. H1 success does not establish H2, H3 or H4. Better abundance retrieval (`records/heldout_v2/`) establishes none of them.
- A target is **independent** only if it was not used in any fitting, selection or calibration of the model being scored. That includes nuisance-state inference, hyper-parameter selection, discrepancy fitting and the choice of variant.
- Fitted discrepancy, simulated truth or model output is never a validation target.
- When independent data are inadequate, the hypothesis is reported as **UNTESTED**.
- Every model gets the same inputs, masks, periods and training opportunities. If one needs more fitted parameters or more information, that is reported.
- Uncertainty is computed by block resampling over the dependence unit stated per hypothesis. With fewer than 5 blocks, the result is "insufficient evidence for a generalisation interval" and the per-block results are listed instead.

## 1. Models compared (identical non-algal treatment within each hypothesis)

| ID | Model | Algae optics | Other state |
|---|---|---|---|
| M0 | no algae | none (algae concentration 0) | same ice and dust state as the others |
| M1 | established empirical algae optics | `tierA_empirical`: BioSNICAR `ice_algae_empirical_Chevrollier2023` cell optics | same |
| M1b | measured in vivo pigment optics | `measured_mac_C`: Williamson et al. (2020) MACs + packaging | same |
| M2 | empirical satellite-albedo approach (H2/H4 only) | none: Sentinel-2 BOA → broadband albedo by a published narrowband-to-broadband conversion, fixed coefficients, no local fitting | — |
| M3 | pigment/cell-informed (primary) | `tddft_D` (frozen calibration: AR(1), posterior mean, as in `records/heldout_v2`) | same |
| M3-alt | structural alternatives | `tddft_C`, `tddft_D_iid` (reported, never selected) | same |

- **M1 is not claimed to be a reproduction of any published host model** (SNICAR-ADv4, E3SM, MAR). Claims of improvement over a named host require reproducing that host's documented configuration, which is not part of this protocol.
- **M2 conversion:** the Sentinel-2 narrowband-to-broadband coefficients of Naegeli et al. (2017, Remote Sens. 9:110), as published. If those coefficients cannot be obtained verbatim from the source, M2 is a four-band irradiance-weighted mean with its formula stated, and is labelled "simple empirical", not "published".

## 2. H1 — improved independent surface-albedo prediction (plot scale)

- **Claim tested:** with abundance *given* (measured counts), M3 predicts measured spectral and broadband albedo better than M1 and M1b.
- **Data:** S6, 13–25 July 2017, hemispherical spectral albedo (ASD with cosine receptor), `Albedo_master.csv` in biosnicar-py (Cook et al. 2020 archive; revision fe74eeef). Co-located counts from `cook2020_archive_cell_counts.csv`.
- **Test population:** every plot with a finite albedo spectrum over 350–2500 nm AND a count (expected n ≈ 47, including zero counts).
  - Excluded: spectra with more than 5 % missing values, or any albedo value > 1.05 or < 0 between 400 and 1300 nm (instrument failure). Each exclusion is listed with its reason.
- **Independence caveat (declared before results):**
  - The HCRF spectra of the SAME S6 2017 plots were used to select σ and the radius prior in earlier experiments (`phase4/field_validation.py`). The albedo spectra were never used in any fit.
  - Within H1 the albedo target is new, but the surfaces are not new. H1 is therefore "new target quantity, previously seen sites", not "transfer to a new site".
  - Transfer is untestable for albedo: no albedo spectra with counts exist for S Greenland.
- **Nuisance state (ice radius, density, dust):** not measured per plot.
  - Each model gets the same population prior for these: the S6 dust prior and the ice-radius population distribution estimated by empirical Bayes on the TRAINING days' albedo spectra only.
  - Prediction for a test plot = model albedo at its measured abundance, marginalised over the nuisance prior. Nothing is fitted to the scored spectrum.
- **Folds:** leave-one-day-out over the sampling days (one block per day, about 9 blocks).
- **Primary endpoint:** mean absolute error of broadband albedo, measured vs predicted.
  - Both are integrated over 350–2500 nm with the same weights: the BioSNICAR clear-sky irradiance at the plot's SZA; the model's spectral albedo is evaluated on the measurement's wavelengths.
  - Contrast: MAE(M1) − MAE(M3), and MAE(M1b) − MAE(M3).
- **Minimum meaningful improvement:** 0.01 broadband albedo. Rationale:
  - It is comparable to the stated accuracy of field cosine-receptor albedo, about ±0.01–0.02 (levelling, receptor cosine error).
  - It corresponds to about 3 W m⁻² of daily-mean absorbed shortwave at 300 W m⁻² daily-mean SW↓, i.e. about 30 % of the retrieved daily algal forcing at S6 (10 W m⁻²).
- **H1 supported only if all of the following hold:**
  - the day-block bootstrap 95 % interval of the MAE contrast against M1 excludes 0;
  - its point estimate is ≥ 0.01;
  - the M3 predictive 90 % intervals cover 80–98 % of plots (neither collapsed nor vacuous).
- **Secondary (reported, not decisive):**
  - bias;
  - spectral RMSE over 400–700 nm and over 700–1300 nm;
  - CRPS of broadband albedo;
  - log predictive density;
  - residual spectra;
  - performance stratified by count tercile and for zero-count plots;
  - M0 as a floor.

## 3. H2 — improved absorbed shortwave under identical forcing (station scale)

- **Claim tested:** albedo predicted from Sentinel-2 + model, combined with MEASURED SW↓, predicts measured absorbed shortwave (SW↓ − SW↑) at PROMICE ablation stations better than the alternatives.
- **Data:** PROMICE KAN_L (670 m) and KAN_M (1270 m) hourly data (GEUS Dataverse doi:10.22008/FK2/IW73UU), and Sentinel-2 L2A scenes over each station, June–August 2016–2019.
- **Test population:** station-days with a Sentinel-2 acquisition meeting all of:
  - station pixel valid (no cloud, shadow or snow in SCL, SCL = 11 within 30 m);
  - PROMICE snow_height < 0.02 m;
  - tilt-corrected `dsr_cor`/`usr_cor` available. 2019 is excluded from the primary endpoint because `tilt_y` was missing (P4-SEB-2) and is reported separately.
- **Target:** measured absorbed SW over ±1 h around the overpass (hour-start convention) AND as a daily mean.
  - Predictions use the measured SW↓ of the same hours.
  - The station's own SW↑ and albedo are NEVER inputs to any prediction.
- **Primary endpoint:** MAE of absorbed SW over ±1 h of overpass (W m⁻²), contrasts M1 − M3 and M2 − M3.
- **Blocks:** station-year (up to 6 blocks for 2016–2018).
- **Minimum meaningful improvement:** 10 W m⁻² at overpass. That is about the radiometer uncertainty (about 3–5 % of a SW↓ of 600–700 W m⁻²) and about 1.5 % of SW↓.
- **H2 supported only if:** the block interval excludes 0, the point estimate is ≥ 10 W m⁻², and H1 is not contradicted.
- **Footprint caveat:** the station radiometer footprint (tens of m²) differs from a 10–20 m pixel. Disagreement is reported against the within-pixel heterogeneity measured at S6 (0.83 dex in abundance).

## 4. H3 — improved observed surface ablation within a validated SEB

- **Claim tested:** using M3 albedo instead of M1/M2/M0 albedo in the same SEB (`phase4/seb.py`) improves the predicted surface lowering between consecutive usable scenes.
- **Prerequisites (otherwise H3 = UNTESTED):**
  - the SEB with MEASURED albedo reproduces observed lowering within the ablation-sensor disagreement on the same station-years;
  - two independent ablation records (pressure transducer and stake sonic) agree within 25 % over the interval.
- **χ:** the subsurface shortwave fraction is not fitted on any H3 test interval. χ = 0 is the primary setting; χ = 0.3 is a sensitivity scenario fitted on 2016–2018 KAN_M and therefore NOT independent there.
- **Endpoint:** absolute error of lowering in m w.e. per interval, against EACH ablation record separately (no record selection).
- **Expected status:** the 2019 KAN_M records disagree by 2.5×, and the measured-albedo SEB over-predicts 2016–2018 by 1.5–1.7×. H3 is therefore expected to be UNTESTED or inconclusive at KAN_M. KAN_L is assessed against the same prerequisites before any algae model is run there.

## 5. H4 — satellite algae interpretation adds information beyond empirical satellite albedo

- **Claim tested:** for the H2 station-days, broadband albedo from the retrieval + forward model (M3) predicts the station's measured albedo better than the empirical narrowband-to-broadband satellite albedo (M2) from the same pixel.
- **Target:** measured station albedo over ±1 h of overpass (`usr_cor`/`dsr_cor`, hourly means).
- **Endpoint:** MAE of broadband albedo, contrast M2 − M3, with the same blocks and evidence rule as H2.
- **Minimum meaningful improvement:** 0.01.
- **Additional requirement:** if M3 is no better than M2, the satellite algae interpretation adds no albedo-predictive information at these stations; the abundance product is then reported only for its own (diagnostic) value.

## 6. Uncertainty (all hypotheses)

- **Primary analyses:** frozen posterior-mean calibration (as in `records/heldout_v2`).
- **Follow-up uncertainty analysis:** joint calibration draws carried through optics → radiative transfer → retrieval → albedo/SW, with the SAME draw index across all pixels and sites of a draw.
  - It is reported next to the primary analysis, not instead of it.
  - Structural scenarios (AR(1) vs iid residuals, χ, Fe term) are kept separate.
- **Proper scores** (CRPS, log density) are reported together with sharpness (mean interval width) and coverage. Coverage alone is never reported as calibration.

## 7. Data manifest

All datasets, versions, checksums, licences, uses and conflicts: `docs/data_manifest.md`. Every use of a dataset in fitting, selection or evaluation is recorded there.
