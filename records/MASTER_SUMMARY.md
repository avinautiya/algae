# Master summary of quantitative findings (for the Regeneron STS write-up)

Every number below comes from a committed file in `records/` or `phase1/results/`; the source is given
per block. "Pending" marks runs still in progress (status at the end). Units: abundance in log₁₀ cells
mL⁻¹ (meltwater), forcing in W m⁻², melt in cm water equivalent per day.

## 1. Molecular level: does TD-DFT reproduce the pigment? (Phase 1 + calibration)

Source: `phase1/results/level2/`, `records/phase1_production/calibration_level2_B3LYP.json`, Fig. S0.

| Quantity | Value |
|---|---|
| Method | TD-B3LYP/6-31G*/IEF-PCM (water, non-equilibrium), 30 singlets, at a B3LYP/PCM geometry (optimisation stopped at ΔE = 5 × 10⁻⁶ Eh; `stop_criterion.json`) |
| Lowest bright states | 420 nm (f 0.074), 381 nm (0.147), 333 nm (0.198), 296 nm (0.271) |
| Fit to HPLC shape of the isolated pigment (R²) | **0.974** (placeholder spectrum: 0.91) |
| Band shift ΔE | **+0.060 ± 0.004 eV** (B3LYP slightly too red; no charge-transfer error) |
| Band FWHM | 0.62 ± 0.005 eV |
| Strength factor f | 39.5 ± 2.2 |
| f / stoichiometric glucoside→phenol (4.53) | **8.7** (not ≈ 1: mainly the 4-AAP assay's low response to this pigment, which EPA 420.1 reports as a minimum; cancels in the forcing) |
| Fe-complexed fraction φ (fit to the S6 extract) | 0.49 ± 0.04 |
| Fe share of 400–700 nm absorption (tier D) | 53 % |
| CAM-B3LYP at the same geometry | **pending** |
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
| Satellite scale, tier D | **pending** (rerun after imagery-bucket fix) | |

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
- Our high-biomass melt brackets Cook's value. At lower abundance we are 2–7× lower, which is a discussion point.

## 6. Seasonal map series (2019, six clear scenes)

**Pending** (`phase4/multi_scene.py`; rerun after imagery-bucket fix).

## 7. Stated limitations (each quantified)

1. **Dust.** No dust measurement at the independent site: ±0.2 dex on its abundance.
2. **Pixel scale.** 10 m pixels average 0.83 dex of plot-to-plot variability. Satellite products are area means.
3. **Strength factor.** f / stoichiometric = 8.7 separates assay response from TD-DFT intensity error only jointly. This doesn't affect the forcing.
4. **Geometry.** The Level 2 geometry was stopped at ΔE 5 × 10⁻⁶ Eh between steps, far below TD-DFT error.
5. **Irradiance.** PROMICE KAN_M irradiance is used for S6, about 20 km away (literature comparison).
6. **Stibal et al. 2017.** Not included: the supporting information needs a browser download (`data/empirical/stibal_2017/README.md`).

## Status of runs (as of this file's last commit)

- Phase 1: CAM-B3LYP Level 2 TD-DFT running. Queue afterwards: Level 1 (B3LYP + CAM-B3LYP), Fe(III) catecholate, Fe(III) tropolonate + spin checks.
- Downstream: satellite check (tier D) and 2019 seasonal series queued (`scripts/production_downstream.sh`).
