# Empirical inputs: sources, licences and provenance

Every model input is read at run time from a file in this folder, or computed from those files by a documented operation in `phase2/empirical_data.py` or `phase2/tddft_calibration.py`. A few published values are hard-coded there (the species dimensions in Procházková et al. 2021, Table 2), with their citation. Every quote below was checked against the downloaded source text.

## Files in this folder

| File | Content | Source | Licence |
|---|---|---|---|
| `williamson2020_pigments_per_cell.csv` | Mean (SD, SE) pigment mass per glacier-algal cell (ng cell⁻¹): phenolics, chl a, chl b, 6 carotenoids. 53 surface-ice samples, S6, 2016 | Williamson et al. (2020) deposit, UK Polar Data Centre GB/NERC/BAS/PDC/01248, `average_pigment_concentrations.csv` (unchanged) | Open Government Licence v3 |
| `williamson2020_mass_absorption_coefficients.csv` | In vivo MACs (m² mg⁻¹, 250–750 nm): chl a, chl b, carotenoids (ppc/psc), phenolic extracts | Same deposit, `mass_absorption_coefficients.csv` | OGL v3 |
| `williamson2020_phenolics_mac_with_error.csv` | Phenolic extract MAC with its regression standard error (53 samples) | Same deposit, `phenolics_mac.csv` | OGL v3 |
| `williamson2020_phenolics_hplc_abs_spec.csv` | Diode-array absorption spectra of the four main HPLC peaks of the phenolic extract | Same deposit, `phenolics_hplc_abs_spec.csv` | OGL v3 |
| `williamson2020_phenolics_per_sample.csv` | Per-sample phenolic concentration, cells mL⁻¹, biovolume per cell, phenol per cell (64 samples) | Same deposit, `phenolic_pigmentation_concentrations.csv` | OGL v3 |
| `williamson2020_count_biomass_S6_2016.csv` | Cell abundance and biovolume per sample and habitat, S6, 2016 | Same deposit, `count_biomass_all_2016.csv` | OGL v3 |
| `prochazkova2025_fig4_PG_PGFe_absorbance_digitized.csv` | Absorbance of purpurogallin (PG) and Fe-purpurogallin (PG-Fe) in water, 255–800 nm, 1 nm | Digitised from Procházková et al. (2025), Fig. 4 (image file `EMI4-17-e70149-g007.jpg`) | CC BY 4.0 (article) |
| `cook2020_archive_cell_counts.csv` | 47 S6 samples (13–24 July 2017): cells counted (80 large squares), cells mL⁻¹ | Cook et al. (2020) data archive, Zenodo 10.5281/zenodo.3564501, `Cell_Counts/Updated_GrIS_Cell count results_250518.xlsx`, sheet "Cells per ml" (SB samples with a spectrum) | Open (archive README asks for citation of Cook et al. 2020) |
| `cook2020_archive_hcrf.csv` | Nadir HCRF spectra of the same samples, 350–2499 nm | Same archive, `Albedo_Reflectance_Processed/HCRF_master.csv` | as above |
| `cook2020_field_metadata.csv` | The published BioSNICAR-GO retrievals (`algae_cells_inv_model_S2`) for 31 S6 samples, used only as a literature baseline | biosnicar-py `data/additional_data/Spectra_Metadata.csv` | MIT |
| `chevrollier2023_tableS4.csv` | 26 rows of Table S4: measured ice- and snow-algae cells mL⁻¹, solar zenith of each spectrum, and the authors' model parameters | Parsed from the text of the Supplementary Material of Chevrollier et al. (2023), J. Glaciol. 69:333, doi:10.1017/jog.2022.64 | CC BY 4.0 |
| `chevrollier2023_hcrf.csv` | ASD FieldSpec4 HCRF spectra, 350–2500 nm, of the 20 southern-Greenland samples (5–6 Aug 2021) | Zenodo 10.5281/zenodo.18826013 (`HCRF_ASDFieldspec4_spectra_Chevrollier_et_al_2026.nc`, tags `<id>_HCRF`) | CC BY 4.0 |
| `ice_ssa_measurements.csv` | 19 measured specific surface areas of bubbly ice | Cooper et al. (2021); Dadic et al. (2013), micro-CT data file for Fig. 8 (UW ResearchWorks) | values from publications |
| `biosnicar_field_ARF.csv` | Anisotropic reflectance factor (HCRF/albedo) spectra, 51 field sites | biosnicar-py `data/additional_data/ARF_master.csv` | MIT |
| `esa_s2_srf_TN-15-0007_v4.0.csv` | Official Sentinel-2A/B/C spectral response functions, 1 nm | ESA COPE-GSEG-EOPG-TN-15-0007 v4.0 (SentiWiki) | Copernicus open licence |
| `promice_KAN_{L,M,U}_day_2019.csv` | Daily air temperature and altitude, 2019 | PROMICE, GEUS Dataverse doi:10.22008/FK2/IW73UU | CC BY 4.0 |
| `promice_KAN_M_hour_JJA_radiation.csv` | Hourly downwelling shortwave (`dsr`, `dsr_cor`), cloud cover `cc`, air temperature, June–August of all years (24 736 h) | Same dataverse, `KAN_M_hour.csv` (columns and months subset) | CC BY 4.0 |

## The former gaps, and what replaced them

### 1. TD-DFT error ranges → empirical calibration of the Phase 1 spectrum (`phase2/tddft_calibration.py`)

The phenolic pigment of the field algae has measured spectra in the Williamson et al. (2020) deposit:

- **Isolated (uncomplexed) chromophore.** Diode-array spectra of the HPLC peaks. Peaks 2–4 share the purpurogallin-type spectrum. Williamson et al. (2020) report it as "identical absorbance features in both the UV and across the visible spectrum (λmax = 304 nm, secondary peak at λ = 389 nm)".
- **Whole extracts.** The MAC per mass of phenolics, 250–750 nm, with its standard error.

**Units.** Phenolics were quantified with US EPA Method 420.1 (4-AAP), which the PNAS paper cites as its reference 49. So the per-cell phenolic mass and the extract MAC are both in **phenol equivalents**. Before this calibration, tiers B–D multiplied a MAC per kg of glucoside by a concentration in phenol equivalents.

**Method (cut posterior).**
- **Stage 1:** the band shift ΔE and Gaussian FWHM are fitted to the HPLC shape (emcee, with a fitted discrepancy SD).
- **Stage 2:** the scale f and the Fe-complexed fraction φ are fitted to the extract MAC in log space (maximum likelihood, Gaussian approximation, φ ∈ [0, 1]).
- **Role of f:** it absorbs the oscillator-strength error and the glucoside-to-phenol-equivalent conversion.
- **Output:** the Phase 3 molecular PDFs are the marginal posteriors.

### 2. Tier D Fe-phenolic stand-in → measured Fe-purpurogallin spectrum

Procházková et al. (2025), Environ. Microbiol. Rep. 17:e70149, doi:10.1111/1758-2229.70149, Fig. 4: "Comparison of the spectral absorbance between the purpurogallin standard (PG, green) and condensed iron‐purpurogallin (PG‐Fe, red), both dissolved in water."

**Digitisation.** Axis ticks were located in the image, curve pixels were identified by colour, and the median of each pixel column was taken. The result matches the published curves to within the line width.

**How it enters the model.** Tier D = f·[M(λ) + φ·I_M·D(λ)]:
- D is the measured Fe-induced change of PG absorbance per unit integrated PG absorbance.
- I_M is the integral of the calibrated chromophore MAC over the same window.
- φ is fitted to the S6 extract MAC: whole extracts are dark (they contain the complex), while the isolated peak is "only yellowish" (Procházková et al. 2025).

**Stated assumption.** The two solutions in Fig. 4 contain the same amount of PG. The paper gives no concentrations. Their equal maxima (1.04 and 1.05) are consistent with this, but could also mean the curves were normalised. φ absorbs any residual scale difference.

### 3. Equal intracellular concentration → measured size dependence

Williamson et al. (2020) give per-sample phenolics per cell and mean biovolume per cell (64 samples).

**Fit.** ln(pg cell⁻¹) = a + b·ln(V/V_ref), with variance = Poisson counting error (1/N cells counted) + s₀², fitted by maximum likelihood. Concentration then scales as V^γ, with γ = b − 1.

**Result:** γ = −0.73 ± 0.73 (jackknife SE).
- **Why jackknife:** the estimate is sensitive to a few low-abundance samples. Leaving one out moves it between −1.3 and −0.6. Keeping only samples with ≥ 300 cells mL⁻¹ gives +0.4.
- **Phase 4 (point estimate):** c(A. nordenskioeldii) = 19.6 kg m⁻³ and c(A. alaskanum) = 41.6 kg m⁻³.
- **Phase 3:** propagates γ ~ N(−0.73, 0.73).
- **Interpretation:** the data are consistent with equal concentration (γ = 0, within 1 SE). An independent site-level regression on Halbach et al. (2022) Tables S4/S6 points the same way (alaskanum/nordenskioeldii ≈ 1.6), also not significantly.
- **Earlier evidence**, from Halbach et al. (2022), Sci. Rep. 12:17643, doi:10.1038/s41598-022-22271-4: "The log-transformed purpurogallin content per cell showed a positive correlation with biovolume, though not significantly (R 2 = 0.49, p = 0.074)".

### 4. Ice-radius prior bounds → measured SSA, and the field spectra

**BioSNICAR's bubbly-ice scattering.** Per kg of ice it is σ_air(r)·V_air/ρ (`column_OPs.get_layer_OPs`), so the radius that matches a measured specific surface area is r = 3(1 − ρ/917)/(ρ·SSA).

**Measured SSA** (`ice_ssa_measurements.csv`):
- **Cooper et al. (2021)**, The Cryosphere 15:1931, doi:10.5194/tc-15-1931-2021. Site: "∼ 1 km from the ice sheet margin at 840 m above sea level. (67.15° N, 50.02° W)". The values are for ice 12–124 cm deep: "The optimal reff values are ∼ 9.3 and ∼ 10.6 mm with corresponding specific surface areas ∼ 0.35 and ∼ 0.31 m2 kg−1". These r_eff are radii of ice spheres in air, **not** bubble radii; only the SSA carries over.
- **Dadic et al. (2013)**, JGR Earth Surf. 118:1658, doi:10.1002/jgrf.20098. Micro-CT of Antarctic blue and white ice, 17 subsamples with ρ > 830 kg m⁻³. For comparison: "The SSA of ice at 10 m depth, which has not been subjected to seasonal temperature cycles, is 0.16 m2 kg–1".

**Measured-SSA prior.** ln SSA ~ N(−0.97, 0.35), with the two studies weighted equally. At 690 kg m⁻³ that is a median bubble radius of 2.8 mm.

**Field-spectra prior (Phase 4).** These measurements characterise ice 0.1–1 m deep, not the weathering crust at the surface. The field validation therefore also estimates the population distribution of the surface radius from the 64 field spectra, by empirical Bayes: log-normal hyperparameters jointly with σ, maximising the marginal likelihood of the other samples. The data prefer it strongly over the measured-SSA prior (log evidence 230 vs 183 for our model, 366 vs 327 for Tier A).

### 5. Phase 3 ice-radius PDF from our own retrievals → measured SSA

Phase 3 now samples SSA ~ logN(−0.97, 0.35) (measurements only) and converts it at the bottom-layer density.

### 6. Fixed choices → measurements

**Cell asymmetry parameter g.** This is now Mie theory for the equal-volume sphere (`cell_optics.mie_g`), with:
- the measured cell refractive index, 1.38. Chevrollier et al. (2023) SI: "…indicating a maximum transmission for a refractive index of 1.38, really close to theoretical estimations (Dauchet and others, 2015). For ice algae, this could not be measured and we used the same value of 1.38.";
- the host medium, ice (n = 1.31);
- k from each cell's own absorption.

The result is g ≈ 0.98–0.99. This agrees with measurements on green microalgae: "the associated asymmetry factor gλ was larger than 0.95 and did not change significantly with wavelength" (Pilon & Kandilian 2016, Adv. Chem. Eng. 48, doi:10.1016/bs.ache.2015.12.002, citing Kandilian et al. 2013). BioSNICAR's 0.96 is undocumented.

**Clear-sky transmissivity.** T = 0.919 ± 0.036. It is fitted to 2301 clear-sky hours (PROMICE cloud cover ≤ 0.1, SZA ≤ 75°) of tilt-corrected downwelling shortwave at KAN_M, June–August of all years. The former 0.75 underestimated SW↓, and therefore forcing in W m⁻², by about 35 % at SZA 50°. Phase 3 propagates T ~ N(0.919, 0.036).

**2 cm algal layer.** This is the sampling depth that defines cells mL⁻¹ in the count data used for the priors:
- Williamson et al. (2018): "the top 2 cm collected using a metal ice saw and trowel";
- Halbach et al. (2025): "surface ice was collected by scraping off the top ~2 cm".

So the model layer matches the definition of the measured abundance. Two caveats:
- **Chevrollier et al. (2023)** scraped "1–6 cm" (depth not recorded per sample), which adds unquantified error to the southern-Greenland validation counts.
- **Vertical distribution:** how the algae are distributed within the 2 cm layer is not measured; the model treats it as uniform.

### 7. Stibal et al. (2017) → still not obtained; two other datasets added

The Stibal et al. (2017) supporting-information table (`grl56634-sup-0002-2017gl075958_data_si.xlsx`) is only on Wiley/AGU, which returned HTTP 403 to every automated route. It needs a browser download or a request to the authors. The validation was instead extended from 31 samples (the biosnicar copy) to 64 samples with spectra:
- **S6 2017:** 46 samples, from the primary Cook et al. (2020) count workbook, with Poisson counting errors.
- **Southern Greenland 2021:** 18 counted samples (Chevrollier et al. 2023), an independent site, year, team and instrument.

Two samples, 22_7_SB6 and 22_7_SB7, appear in biosnicar-py's metadata with 0 cells but are **not** in the count workbook. They were never counted, so they are not used. The archived spectrum of 20_7_SB4 is empty.

## Other inputs (unchanged from the previous revision)

**Species volumes.** Halbach et al. (2022): "…in almost equal proportions of both species with average cell volumes of 2307 ± 421 and 822 ± 133 µm 3 , respectively."

**Species shapes.** Procházková et al. (2021), Microorganisms 9:1103, Table 2, mean length:width:
- A. nordenskioeldii: 2.37, giving 10.75 × 25.44 µm.
- A. alaskanum: 1.46, giving 8.95 × 13.08 µm.

**Pigments per cell and in vivo MACs.** From the Williamson et al. (2020) deposit, used as in their script (`photoinhibition_main_script.Rmd`, lines 275–295).

**Abundance prior.** log₁₀ B ~ N(3.559, 0.778), from 180 natural surface-ice samples at S6 in 2016.

**Community-fraction prior.** f_n ~ Beta(17.37, 11.42), a moment fit to three published values:
- Williamson et al. (2018): "Ancylonema typically demonstrated the greatest relative abundance (∼65%)";
- Halbach et al. (2025): "~66% were filamentous Ancylonema cf. nordenskiöldii";
- Halbach et al. (2022): "almost equal proportions" (taken as 0.50).

**Reflectance factor k.** N(0.898, 0.175), from 51 field ARF spectra.

**Density.** Cooper et al. (2018), The Cryosphere 12:955, doi:10.5194/tc-12-955-2018:
> "even lower density (0.33–0.56 g cm−3, µ = 0.45 g cm−3) unsaturated weathering crust"

> "low-density (0.43–0.91 g cm−3, µ = 0.69 g cm−3) ice to at least 1.1 m depth"

## What remains assumed (and is stated where it is used)

1. Equal PG amount in the two solutions of Procházková et al. (2025), Fig. 4 (see 2).
2. The size-scaling exponent γ is applied to chlorophylls and carotenoids as well, because only phenolics have per-sample size data.
3. The equal-volume sphere for Mie g, and the refractive index of 1.38, which was measured on snow algae and applied to ice algae (as in Chevrollier et al. 2023).
4. A uniform vertical distribution of algae within the 2 cm layer.
5. The same bubble radius in both ice layers. The field-spectra prior describes the effective (surface-dominated) radius.
