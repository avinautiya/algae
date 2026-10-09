# Master summary of quantitative findings (for the Regeneron STS write-up)

Every number below comes from a committed file in `records/` or `phase1/results/`; the source is given
per block. "Pending" marks runs still in progress (status at the end). Units: abundance in log₁₀ cells
mL⁻¹ (meltwater), forcing in W m⁻², melt in cm water equivalent per day.

## 1. Molecular level: does TD-DFT reproduce the pigment? (Phase 1 + calibration)

Source: `phase1/results/level2/`, `records/phase1_production/calibration_level2_B3LYP.json`, Fig. S0.

| Quantity | Value |
|---|---|
| Method | TD-B3LYP/6-31G*/IEF-PCM (water, non-equilibrium), 30 singlets, at a B3LYP/PCM geometry that did **not** meet the optimisation criteria (1 of 5 met; `stop_criterion.json`). Status `provisional_geometry`; spectral effect of relaxation pending (jobs L2_OPT, L2_B3LYP_TDA15_RELAXED) |
| Lowest bright states | 420 nm (f 0.074), 381 nm (0.147), 333 nm (0.198), 296 nm (0.271) |
| Fit to HPLC shape of the isolated pigment (R²) | **0.974** (placeholder spectrum: 0.91) |
| Band shift ΔE | **+0.060 ± 0.004 eV** (B3LYP slightly too red; a uniform fitted shift does not establish the absence of charge-transfer error) |
| Band FWHM | 0.62 ± 0.005 eV |
| Strength factor f | 39.5 ± 2.2 |
| f / stoichiometric glucoside→phenol (4.53) | **8.7** (not ≈ 1: mainly the 4-AAP assay's low response to this pigment, which EPA 420.1 reports as a minimum; cancels in the forcing) |
| Fe-complexed fraction φ (fit to the S6 extract) | 0.49 ± 0.04 |
| Fe share of 400–700 nm absorption (tier D) | 53 % |
| CAM-B3LYP at the same geometry (TDA, 15 roots; compared with B3LYP-TDA, 15 roots) | **pending**: band-position sensitivity check only; TDA intensities are not used downstream |
| Level 1 (purpurogallin core), both functionals | **pending** |
| Fe(III)–purpurogallin complexes (TD-DFT, spin check) | **pending** |

## 2. Molecular-to-field test: which pigment optics explain the field spectra? (central result)

Source: `records/pigment_path_comparison/path_comparison_metrics.csv`. 59 counted samples (41 at S6 2017, 18 in southern Greenland 2021), leave-one-out, measured S6 dust prior, σ / radius prior / τ re-selected per path.

| Pigment optics | Log evidence | Bias (all) | S6 bias | S Greenland bias | RMSE | Spearman ρ | 95 % coverage: posterior / with τ |
|---|---|---|---|---|---|---|---|
| **TD-DFT tier D** (calibrated + Fe complex) | **288.4** | +0.16 | +0.25 | **−0.06** | 0.49 | 0.76 | 0.83 / 0.98 |
| Measured extract MAC | 284.0 | +0.05 | +0.16 | −0.21 | 0.48 | 0.76 | 0.83 / 0.98 |
| TD-DFT tier C (no Fe) | 189.8 | −0.31 | −0.11 | −0.78 | 0.71 | 0.64 | 0.83 / 0.97 |
| (Tier D without dust) | 233.1 | +0.42 | +0.52 | +0.20 | 0.60 | 0.78 | 0.75 / 0.97 |

- **Tier D vs measured MAC:** Bayes factor e^4.4 ≈ 80 in favour of tier D.
- **Tier D vs tier C:** e^98 in favour of tier D. Iron complexation is required to explain the field spectra.
- **Caveat (BioSNICAR's empirical measured-cell optics, Tier A).** Tier A fits the field spectra better still (log evidence 391.9, σ = 0.010), but recovers the counted abundance worst: bias +0.56 dex, RMSE 0.71, against +0.16 and 0.49 for tier D.
  - Spectral evidence alone therefore does not rank abundance accuracy. Tier D is the best optics for abundance among those tested.
  - Its residual misfit is consistent with its red tail. Above 620 nm tier D cells stay dark, whereas measured cells become transparent (optics audit, ratio 1.6 at 650–700 nm).

## 3. Field validation, calibration and site effects (Phase 4)

| Finding | Value | Source |
|---|---|---|
| Structural model error τ (LOO max likelihood) | 0.42 dex (tier D), 0.37 (measured MAC) | path comparison; README "Structural model error" |
| Coverage of 95 % intervals | 0.83 → 0.98 with τ | same |
| Mineral dust (measured S6 prior) | Raises log evidence by 55–70 for every path (D +55, measured +57, C +70); S6 bias +0.52 → +0.25 (tier D) | path comparison |
| Dust at the independent site (not measured there) | 4 bands do not constrain dust. The southern bias ranges from −0.28 to +0.19 dex across the dust priors, i.e. ±0.2 dex structural | `records/dust_sensitivity_williamson2020/` |
| Algal film (top 2 mm) | Disfavoured (log evidence −20) | bias study |
| Site calibration factor | Rejected: removes S6 bias, but drives the independent site to −0.30 to −0.37 dex | bias study |
| Satellite scale, S6 2017 (measured MAC) | Pixel RMSE 0.74 dex ≈ within-pixel spread of counts 0.83 dex; site-mean agreement within 0.04–0.38 dex (count SE 0.24–0.42); coverage 1.00 | `records/satellite_validation_williamson2020/` |
| Satellite scale, tier D (production) | Pixel bias −0.03, RMSE 0.73 dex (floor 0.83); site means within 0.09–0.27 dex of the day-mean counts (count SE 0.24–0.42); coverage with τ 1.00 | `records/satellite_validation_tddft_tierD/` |

## 4. Global uncertainty and sensitivity (Phase 3, production)

Source: `records/phase3_full/`. 1000 Latin-hypercube samples; Saltelli design N = 1024 (28 672 runs); instantaneous clear-sky forcing at SZA 45°.

**Monte Carlo** (median, 95 % range):
- Forcing:
  - Tier D: 12.1 W m⁻² (0.33–211)
  - Tier C: 7.4 (0.21–131)
  - Tier A (BioSNICAR default algae): 10.8 (0.31–184)
- Efficiency, W m⁻² per 10⁴ cells mL⁻¹:
  - Tier D: 32.0 (15.7–51.4)
  - Tier C: 19.3 (9.9–30.7)
- Pigment packaging effect (C − B): −55 W m⁻² (−113 to −3). Packaging removes most of the absorption the same pigment would have dissolved.
- Fe-complexation effect (D − C): +4.9 W m⁻² (0.1–79); efficiency +66 %.
- Broadband albedo, tier D: 0.47 (0.23–0.57).

**Sobol total-effect indices S_T** (converged; Fig. S3: S_T(abundance) stable at 0.98 ± 0.07 from N = 512):

| Output | 1st | 2nd | 3rd | Molecular group S_T |
|---|---|---|---|---|
| RF tier D | cell abundance 0.98 | cell volume 0.013 | dust 0.008 | < 0.001 |
| Efficiency tier D | abundance 0.31 | cell volume 0.27 | surface density 0.16 (dust 0.14) | 0.001 |
| Packaging effect | abundance 0.88 | intracellular concentration 0.11 | cell volume 0.03 | – |
| Fe effect | abundance 0.97 | cell volume 0.02 | surface density 0.01 | – |

**Interpretation.** After calibration to measured spectra, the molecular (TD-DFT) uncertainty no longer matters for the forcing: S_T < 0.001. Abundance controls absolute forcing. The per-cell efficiency is set by cell size, surface density and dust.

## 5. Daily forcing and melt vs published estimates

Source: `records/literature_comparison/literature_comparison.csv`. Tier D; MC over the class abundance spread and all Phase 3 PDFs; measured PROMICE KAN_M hourly irradiance on the published day.

| Case | Our daily RF (W m⁻²) | Our melt (cm w.e. d⁻¹) | Published melt |
|---|---|---|---|
| Cook 2020, Hbio (2.9 × 10⁴), 21 Jul 2017 | 29 (8–74) | 0.76 (0.22–1.92) | 1.35 ± 0.01 (RF), 1.37 ± 0.48 (EB); RF 116 |
| Cook 2020, Lbio (4.7 × 10³) | 5.9 (2.1–17) | 0.15 (0.06–0.44) | 1.01 ± 0.01, 0.95 ± 0.41; RF 65 |
| Williamson 2020, high (9.0 × 10³), 26 Jul 2016 | 9.4 (3.6–23) | 0.24 (0.09–0.60) | 1.86 ± 0.99 |
| Williamson 2020, low (186) | 0.13 (0.02–0.97) | 0.003 (0.000–0.025) | 0.03 ± 0.00 |

**Reading.**
- Our forcing is the algae's own packaged absorption relative to dusty clean ice.
- The site-differencing estimates attribute all darkening between algal and clean sites to algae. That includes weathering-crust structure, dust and water.
- They also disagree with each other: Williamson's melt per cell is about 4× Cook's.
- Cook's 116 W m⁻² cannot be a 24 h mean: it would melt 3.0 cm d⁻¹, not 1.35.
- Our high-biomass melt brackets Cook's value. At lower abundance we are 2–7× lower.
- **Optics audit** (`records/optics_audit/AUDIT.md`) found no bug. Our per-cell absorption is 0.93× measured in vivo cells (Chevrollier et al. 2023), and three independent optics give the same darkening to within 5 %. Williamson's 1.86 cm d⁻¹ exceeds the black-cell geometric limit of its own stated cell size and sampling by about 6×. Cook's RF is the full algal-site vs clean-ice albedo difference.

## 6. Seasonal map series (2019, six clear scenes)

Source: `records/multi_scene_2019_tddft/` (Fig. S7).
- Area: 6 × 6 km around S6, 20 m pixels, about 90 000 bare-ice pixels per date.
- Optics: tier D with dust; τ = 0.42 dex.
- Daily mean = overpass RF × (daily mean / overpass-hour SW), using PROMICE KAN_M hourly SW of that day. 2019 has only the uncorrected `dsr`; the ratio is insensitive to a tilt calibration.

| Date | Median log₁₀ B | Median BBA | RF at overpass, ours / Tier A (W m⁻²) | Daily-mean RF, ours | Melt, ours / Tier A (cm w.e. d⁻¹) |
|---|---|---|---|---|---|
| 8 Jul | 3.53 | 0.52 | 22.6 / 35.8 | 11.9 | 0.31 / 0.45 |
| 15 Jul | 3.56 | 0.53 | 24.0 / 29.3 | 13.0 | 0.34 / 0.38 |
| 23 Jul | 3.61 | 0.51 | 26.1 / 41.1 | 13.6 | 0.35 / 0.51 |
| 2 Aug | 3.66 | 0.49 | 27.3 / 43.3 | 11.7 | 0.30 / 0.45 |
| 12 Aug | 3.62 | 0.54 | 23.7 / 32.0 | 10.7 | 0.28 / 0.35 |
| 29 Aug | 3.61 | 0.62 | 16.0 / 13.0 | 6.2 | 0.16 / 0.13 |
| **Season (6 dates)** | | | | **11.2 (6.2–13.6)** | **0.29** (Tier A daily RF 14.6) |

- **Abundance and albedo.** Abundance rises to an early-August peak (median 4.6 × 10³ cells mL⁻¹). Albedo is lowest on 2 Aug.
- **Forcing timing.** Forcing peaks in late July. By late August it falls with the lower sun (SZA 58°) and the brighter surface.
- **Fit.** 0 % of pixels fail the χ² test with either optics. The per-pixel Bayes factors mildly favour Tier A (median ln BF −0.4 to −1.4), so 4-band 20 m pixels cannot discriminate the optics.

## 7. Stated limitations (each quantified)

1. **Dust.** No dust measurement at the independent site: ±0.2 dex on its abundance.
2. **Pixel scale.** 10 m pixels average 0.83 dex of plot-to-plot variability. Satellite products are area means.
3. **Strength factor.** f / stoichiometric = 8.7 separates assay response from TD-DFT intensity error only jointly. This doesn't affect the forcing.
4. **Geometry.** The Level 2 geometry is not converged (1 of 5 criteria met). A small ground-state energy change between steps does not bound excitation-energy error; the effect is being measured by relaxing the geometry (job L2_OPT) and recomputing the spectrum (L2_B3LYP_TDA15_RELAXED vs L2_B3LYP_TDA15).
5. **Irradiance.** PROMICE KAN_M irradiance is used for S6, about 20 km away (literature comparison).
6. **Stibal et al. 2017.** Not included: the supporting information needs a browser download (`data/empirical/stibal_2017/README.md`).
7. **Chemical model.** One neutral tautomer and conformer, 6-31G(d), implicit solvent only; protonation, tautomer, basis and explicit-solvent effects are not measured.

## Status of runs (as of this file's last commit)

- Phase 1: see `python phase1/jobs.py status` and `docs/repair_ledger.md`. The repair audit of commit 814295a found defects in checkpoint, completion and provenance handling; results listed in this file that depend on them are being re-verified, and entries are corrected as each repair lands.
- Downstream: complete (satellite check tier D, 2019 seasonal series, Phase 3 full, literature comparison, optics audit).
