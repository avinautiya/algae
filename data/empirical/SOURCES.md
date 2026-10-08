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

## Unit and stoichiometry audit

Every link from a measured quantity to the model's absorption per cell, and from cell counts to cells per m², is checked:

| Link | Units | Status |
|---|---|---|
| Williamson extract MAC | Regression slope of absorbance on concentration, in L g⁻¹ cm⁻¹ (decadic, 1 cm path, phenol equivalents). Converted to m² mg⁻¹ as slope × 2.3 × 10⁻⁴ (`phenolic_pigmentation_mac_calculation.R`, lines 71–73): 1 L g⁻¹ cm⁻¹ = 10⁻⁴ m² mg⁻¹ decadic, × ln 10 for Napierian | Exact. Used as m² kg⁻¹ (× 10⁶) |
| Per-cell phenolics | ng phenol equivalents per cell (same 4-AAP assay, US EPA 420.1) | MAC × mass per cell: the phenol-equivalent unit cancels exactly, so the `williamson2020` path needs no conversion |
| TD-DFT MAC | m² kg⁻¹ of glucoside (M = 426.33 g mol⁻¹) | The calibrated scale f converts it to per kg phenol equivalent. If the 4-AAP assay counts one glucoside as one phenol, the stoichiometric part of f is 426.33/94.11 = 4.53, and f/4.53 is the TD-DFT intensity error alone (reported as `f_over_stoichiometric` in `tddft_calibration.json`) |
| Fe complex (tier D) | φ is dimensionless. D is per unit integrated PG absorbance, scaled by the integrated glucoside MAC | The MAC stays per kg phenol equivalent. Fe mass is not added to the pigment mass, consistent with the 4-AAP quantification. The fitted product φ·D is invariant to the unknown concentration ratio in Procházková Fig. 4 |
| Chl a, chl b, carotenoids | MAC per mg pigment (HPLC-quantified) × ng pigment per cell | Consistent |
| Intracellular concentration | ng cell⁻¹ / µm³ cell⁻¹ → kg m⁻³ (× 10⁶) | Consistent. Pooled phenolics 0.0432 ng (53 samples) / 1962 µm³ (S6 counts) = 22.0 kg m⁻³. The independent per-sample file gives Σng/Σµm³ = 23.7 kg m⁻³ |
| **Cell counts** | Field counts are cells per mL of **meltwater** (haemocytometer on melted samples; 1 mL = 1 g). BioSNICAR converts its input to cells kg⁻¹ as conc/917 × 10⁶, i.e. per mL of **solid ice** (`column_OPs.mix_in_impurities`) | **Corrected**: the bridge passes conc × 0.917 (`biosnicar_bridge.MELTWATER_TO_BIOSNICAR`), so the column number of cells is the measured count × ρ·dz. Before this, every model abundance was off by 1000/917 (0.04 dex) |
| Dust | ppb = ng g⁻¹ ice. Cook et al. (2020) give µg g⁻¹ ice, which they computed from µg mL⁻¹ "assuming 1 mL of ice to weigh 0.917 g" | Consistent (× 10³) |

## Sampling depth of the S6 counts

Cook et al. (2020) state only that "ice from within the viewing area of the spectrometer was removed using a sterile blade". Their archived discussion manuscript (Zenodo 10.5281/zenodo.3564501, `Peer_Review/Round1/Cook_et_al_Algae_Melting_GrIS_Tracked_changes.pdf`) is explicit about the same measurements: "These measurements were followed immediately by the physical removal of the upper 2 cm of the ice surface within the same patches." So the S6 counts, like Williamson et al. (2018) and Halbach et al. (2025), are cells per mL of the top 2 cm.

## Stated assumptions for the methods section (with measured sensitivity)

1. **Equal purpurogallin amounts in the two solutions of Procházková et al. (2025) Fig. 4.**
   - The paper gives no concentrations. The equal ~315 nm maxima (1.04 vs 1.05) are consistent with equal amounts, but could also mean the curves were normalised.
   - **Effect:** only the product φ·D is fitted to the S6 extract MAC. A different concentration ratio rescales D and the fitted φ inversely, leaving the tier D MAC unchanged (as long as φ stays within [0, 1]). Only the interpretation of φ as "fraction complexed" depends on the assumption. Verified numerically in `phase2/tests/test_phase2.py::test_fe_scale_invariance`.
2. **The size-scaling of intracellular concentration (c ∝ V^γ) is applied to chlorophylls and carotenoids.**
   - Only phenolics have per-sample size data.
   - **Effect:** using the pooled (γ = 0) concentrations for chl/carotenoids instead changes the forcing at 10⁴ cells mL⁻¹ (60 % A. nordenskioeldii, reference ice) by 0.03 % (35.08 vs 35.09 W m⁻²). Phenolics dominate visible absorption in the packaged cells.
3. **Spherical (equal-volume) approximation for the asymmetry parameter g**, with the refractive index 1.38 measured on snow algae and applied to ice algae (as in Chevrollier et al. 2023).
   - **Effect:** Mie g(400–700 nm) is 0.993–0.996 for n = 1.36–1.42 and an ice or water host. Reference forcing changes by < 0.15 % (43.40–43.46 W m⁻²), even with BioSNICAR's fixed 0.96. Cell scattering is negligible next to ice scattering.
4. **Uniform algae through the 2 cm sampling layer.**
   - **Test:** the same number of cells per m² placed in the top 2 mm (`IceSpec(film_dz=0.002)`), compared at the same 3-layer discretisation.
   - **Effect:** darkening is 4–8 % stronger (bubbly or granular ice, 450–650 kg m⁻³), equivalent to a few hundredths of a dex in retrieved abundance. The field-validation effect is in `phase4/bias_study.py`.
5. **One bubble radius for both ice layers.** The field-spectra prior describes the effective (surface-dominated) radius.
6. **Albedo spectra (only if Stibal et al. 2017 is ingested by hand):** k ~ N(1, `albedo_k_sd`), where the default 0.02 represents instrument calibration.

## Sample positions for the satellite-scale check

`tedstone2020_s6_2017_sample_locations.csv` holds 20 positions of counted S6 samples (15, 21, 22 and 23 July 2017). The source is `uav_sb_locations.csv` in Tedstone A. et al. (2020), *Multi-spectral unmanned aerial system imagery, S6, south-west Greenland, July 2017: Levels 2 (ground reflectance) and 3*, UK Polar Data Centre, doi:10.5285/77ca631f-a3a4-4f26-bc90-57bb17baa6fc, OGL v3.

The coordinates are **UTM zone 23N (EPSG:32623)**, not 22N. Read as 22N they fall 6° west on land, with B8 > B2. As 23N they map to 67.0776 N, 49.348 W, the S6 site, and every Sentinel-2 pixel there is classified as snow/ice. `phase4/satellite_validation.py` reprojects the positions to the tile CRS (22WEV, EPSG:32622).

The file keeps the original `allocation_comment` column. Two 15 July positions may be swapped (SB1 ↔ SB5). They are 18 m apart, i.e. in neighbouring pixels.

## Literature comparison: published forcing and melt (quotes)

**Cook et al. (2020)**, The Cryosphere 14:309, doi:10.5194/tc-14-309-2020 (CC-BY 4.0). S6, 21 July 2017.

Abundance classes:
> "Hbio = 2.9×10⁴ ± 2.01×10⁴; Lbio = 4.73×10³ ± 2.57×10³; CI = 625 ± 381; and SN = 0 ± 0 (1 SD)"

Daily forcing and melt:
> "Integrated over the entire day, this indicated a daily mean biological radiative forcing of 116 and 65 W m⁻² for Hbio and Lbio surfaces, respectively … to estimate 1.35 ± 0.01 (standard error, SE) cm w.e. of melting due to algae in Hbio areas on 21 July. For Lbio sites, biological melting on 21 July 2017 was 1.01 ± 0.01 (SE) cm w.e."

Energy-balance cross-check:
> "1.37 ± 0.48 (SE) cm w.e. attributed to Hbio and 0.95 ± 0.41 (SE) cm w.e. attributed to Lbio … 26.15 ± 3.77 % (SE) of the local melting attributed to algae in the Hbio surfaces and 21.62 ± 5.07 % (SE) for Lbio surfaces."

Regional runoff:
> "algal growth led to an additional 4.4–6.0 Gt of runoff from bare ice in the south-western sector of the GrIS in summer 2017, representing 10 %–13 % of the total"

Note: read as a 24 h mean, 116 W m⁻² would melt 116 × 86400 / 3.34×10⁶ = 3.0 cm w.e. d⁻¹, not 1.35. The published forcing and melt are not mutually consistent under that reading. Melt is therefore the like-for-like comparison quantity.

**Williamson et al. (2020)**, PNAS 117:5694, doi:10.1073/pnas.1918412117. S6 area, 26 July 2016.

> "Integration over the complete diel cycle revealed the potential for glacier algal assemblages to contribute from 0.03 ± 0.00 cm w.e.⋅d⁻¹ in low-biomass areas (mean ± SE, n = 27) up to 1.86 ± 0.99 cm w.e.⋅d⁻¹ melt production in high-biomass patches of surface ice (mean ± SE, n = 103)"

with the classes
> "low (186 ± 276 cells⋅mL⁻¹, n = 27), medium (3,711 ± 2,333 cells⋅mL⁻¹, n = 34), or high (8,989 ± 4,773 cells⋅mL⁻¹, n = 103)"

and
> "glacier algae direct only ∼1 to 2.4% of incident energy to photochemistry versus 48 to 65% to ice surface melting"

**Our comparison** (`phase3/literature_comparison.py`) uses the same classes and days. Abundance is log-normal with the class mean and SD; all other inputs come from the Phase 3 PDFs. Forcing is integrated over the measured hourly PROMICE KAN_M irradiance of that day (21 Jul 2017; 26 Jul 2016). KAN_M (1270 m) is about 20 km from S6 (≈1000 m), so the irradiance is approximate for S6 (a stated assumption).
