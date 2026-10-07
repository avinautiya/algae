# Empirical inputs: sources, licences and provenance

Every number that replaces a former placeholder is either (a) read at run time from a file in this folder, or (b) a published value hard-coded in `phase2/empirical_data.py`, with the citation and a verbatim quote below. Derived quantities (means, SDs, fits) are computed in code from these files, so they can be audited and recomputed.

## Files in this folder

| File | Content | Source | Licence |
|---|---|---|---|
| `williamson2020_pigments_per_cell.csv` | Mean (SD, SE) pigment mass per glacier-algal cell (ng cell⁻¹): phenolics, chl a, chl b, 6 carotenoids; 53 surface-ice samples, S6, 2016 | Williamson et al. (2020) data deposit, UK Polar Data Centre GB/NERC/BAS/PDC/01248, file `average_pigment_concentrations.csv` (unchanged) | Open Government Licence v3 |
| `williamson2020_mass_absorption_coefficients.csv` | In vivo mass absorption coefficients (m² mg⁻¹, 250–750 nm) for chl a, chl b, photoprotective (ppc) and photosynthetic (psc) carotenoids, phenolics | Same deposit, `mass_absorption_coefficients.csv` (unchanged) | OGL v3 |
| `williamson2020_count_biomass_S6_2016.csv` | Cell abundance (cells mL⁻¹) and biovolume (µm³ mL⁻¹) per sample and habitat, S6, 2016 | Same deposit, `count_biomass_all_2016.csv` (unchanged) | OGL v3 |
| `cook2020_field_metadata.csv` | 31 surface samples (July 2017, near S6): haemocytometer cell counts (`measured_cells`, cells mL⁻¹) and the published BioSNICAR-GO inversion (`algae_cells_inv_model_S2`) | biosnicar-py `data/additional_data/Spectra_Metadata.csv`, rows with a count | MIT (biosnicar-py) |
| `cook2020_field_hcrf_counted_samples.csv` | Nadir HCRF spectra, 350–2499 nm, of the same 31 samples | biosnicar-py `data/additional_data/HCRF_master_16171819.csv`, matching columns; wavelength column added (the file has 2150 rows = 350–2499 nm) | MIT |
| `biosnicar_field_ARF.csv` | Anisotropic reflectance factor (HCRF / albedo) spectra, 51 field sites | biosnicar-py `data/additional_data/ARF_master.csv` | MIT |
| `esa_s2_srf_TN-15-0007_v4.0.csv` | Official Sentinel-2A/B/C spectral response functions, 1 nm | ESA, COPE-GSEG-EOPG-TN-15-0007 v4.0 (SentiWiki), converted from the xlsx | Copernicus open licence |
| `promice_KAN_{L,M,U}_day_2019.csv` | Daily mean air temperature `t_u` and GPS altitude, April–Sept 2019 | PROMICE AWS data, GEUS Dataverse doi:10.22008/FK2/IW73UU (columns subset) | CC-BY 4.0 |

## Inputs and where each comes from

### Sentinel-2 band shapes (formerly BioSNICAR's top-hat approximation)
- **Now:** ESA's official SRFs (`esa_s2_srf_TN-15-0007_v4.0.csv`), averaged into BioSNICAR's 10 nm bins (`empirical_data.s2_srf_480`). The spacecraft is taken from the scene ID (S2A/S2B/S2C).
- **Band reflectance:** the forward model uses Σ α·SRF·F↓ / Σ SRF·F↓.

### Species cell sizes (formerly placeholders)
- **Volumes:** Chevrollier et al. (2022), Sci. Rep. 12:18101, doi:10.1038/s41598-022-22271-4 (Greenland).
  > "…in almost equal proportions of both species with average cell volumes of 2307 ± 421 and 822 ± 133 µm 3 , respectively."
  (A. nordenskioeldii, A. alaskanum)
- **Shape:** mean length:width of measured populations from Procházková et al. (2021), Microorganisms 9:1103, doi:10.3390/microorganisms9051103, Table 2:
  - A. alaskanum: 8.7×12.1, 7.6×12.1, 8.5×12.8, 8.4×11.4 µm
  - A. nordenskioeldii: 12.9×29.4, 13.6×37.5, 13.7×28.3 µm
- **Resulting cylinders:**
  - A. nordenskioeldii: d = 10.75 µm, L = 25.44 µm (aspect 2.37)
  - A. alaskanum: d = 8.95 µm, L = 13.08 µm (aspect 1.46)
- **Phase 3:** the cell-volume PDF is N(2320, 542) µm³, the mean and SD of per-sample biovolume per cell across 180 S6 2016 surface-ice samples (`williamson2020_count_biomass_S6_2016.csv`). The aspect ratio is uniform between the two species means.

### Intracellular pigment concentration (formerly 50 kg m⁻³)
- **Phenolics:** c_i = (phenolic mass per cell) / (biovolume per cell).
  - 0.04322 ± 0.01716 ng cell⁻¹ (Williamson et al. 2020, n = 53).
  - Pooled S6 biovolume 1962 µm³ cell⁻¹, from the same campaign.
  - Result: **22.0 kg m⁻³**; Phase 3 PDF N(22.0, 8.8), from the per-cell SD.
- **Chlorophyll a, chlorophyll b and carotenoids (Number 3):** computed the same way from the same file.
  - Chl a: 0.00396 ng cell⁻¹.
  - Chl b: 0.00070 ng cell⁻¹.
  - Six carotenoids: 0.00456 ng cell⁻¹ in total.
  - They absorb with the in vivo MACs of `williamson2020_mass_absorption_coefficients.csv`, used exactly as tabulated. Williamson et al.'s own analysis script (`photoinhibition_main_script.Rmd`, lines 275–295) multiplies these MACs by mg cell⁻¹.
  - Carotenoids use the photoprotective ("ppc") spectrum, following Williamson et al., who treated the carotenoid pool of glacier algae as photoprotective.
- **Stated assumption:** the concentration is the same in both species, so pigment per cell scales with volume. Chevrollier et al. (2022) support the direction but not the significance:
  > "The log-transformed purpurogallin content per cell showed a positive correlation with biovolume, though not significantly (R 2 = 0.49, p = 0.074)…"

### Priors of the Sentinel-2 inversion (formerly hand-set)

**log₁₀ B ~ N(3.559, 0.778).** Fitted to the 180 natural surface-ice samples at S6 in 2016 with cells > 0 (`williamson2020_count_biomass_S6_2016.csv`):
- Habitats: h, m, l, clean, dark, dirty, snow, h20log.
- Excluded: incubation treatments ('13', '50', 'full', 'uv'), meltwater/cryoconite ('cryh20', 'discry', 'supra') and 'biofilm'.
- These samples are independent of the 2017 validation samples.

**f_n (A. nordenskioeldii fraction) ~ Beta(17.37, 11.42).** Method-of-moments fit to the three published Greenland surveys:
- Williamson et al. (2018), FEMS Microbiol. Ecol. 94:fiy025: "Ancylonema typically demonstrated the greatest relative abundance (∼65%) during the present study…" Mesotaenium is not modelled, so this is used as the Ancylonema-nordenskioeldii share of the community.
- Halbach et al. (2025), Nat. Commun. 16, doi:10.1038/s41467-025-56664-6: "…of which ~66% were filamentous Ancylonema cf. nordenskiöldii and ~34% were unicellu[lar]…"
- Chevrollier et al. (2022): "in almost equal proportions of both species", taken as 0.50.

**k (HCRF/albedo) ~ N(0.898, 0.175).** Mean and SD over 51 field spectra of the anisotropic reflectance factor, band-averaged with the ESA SRFs (`biosnicar_field_ARF.csv`).

**Bubbly-ice optical radius:** log-uniform on [300, 20 000] µm. There is no measured distribution, so the prior spans the range that reproduces the field NIR reflectances at the measured densities. Cooper et al. (2021), The Cryosphere 15:1931, doi:10.5194/tc-15-1931-2021, measured in this ice:
> "The optimal reff values are ∼ 9.3 and ∼ 10.6 mm…"

**Ice density.** Cooper et al. (2018), The Cryosphere 12:955, doi:10.5194/tc-12-955-2018:
> "thin (<0.5 m), even lower density (0.33–0.56 g cm−3, µ = 0.45 g cm−3) unsaturated weathering crust"

> "low-density (0.43–0.91 g cm−3, µ = 0.69 g cm−3) ice to at least 1.1 m depth"

The model uses a surface layer of 450 kg m⁻³ (Phase 3: U(330, 560)) over 690 kg m⁻³.

**Melt diagnostics (mapped only, not used in the prior).** PROMICE KAN_L/M/U daily air temperature, interpolated linearly in elevation day by day. The fitted 2019 Jun–Jul lapse rate is 5.16 K km⁻¹.

### Field validation data
Cook et al. (2020), The Cryosphere 14:309, doi:10.5194/tc-14-309-2020:
> "…the measurements presented were all made during constant conditions of clear skies at solar noon ±2 h."

> "Samples were vortexed thoroughly before 20 µL was pipetted into a Fuchs–Rosenthal haemocytometer."

`COUNT_RESOLUTION_CELLS_ML = 62.5`: 24 of the 25 non-zero field counts are (integer-rounded) multiples of 62.5 cells mL⁻¹ (the exception is 44 861), so 62.5 is the counting step and zero counts mean below it.

## What is still NOT empirical (disclosed in the READMEs)
1. **TD-DFT error model in Phase 3:**
   - ΔE ~ N(0, 0.075 eV), f-scale ~ logN(0, 0.2), FWHM ~ U(0.25, 0.40) eV.
   - These are user-specified. No experimental spectrum of the purpurogallin glucoside exists to calibrate them.
2. **Tier D Fe(III)–phenolic LMCT surrogate** (ε, λmax):
   - Procházková et al. (2025) report Fe–phenolic complexation but no numeric spectrum.
   - Tier D stays PROVISIONAL until Level 3/4 TD-DFT output replaces it.
3. **Equal intracellular concentration in both species** (see above).
4. **Bounds of the ice-radius prior** (data-spanning, not a measured distribution).
5. **Fixed modelling choices:**
   - Cell asymmetry parameter g = 0.96.
   - Clear-sky SW↓ transmissivity of 0.75 (used for forcing in W m⁻² in Phases 2–3).
   - Algae confined to the upper 2 cm. This follows the BioSNICAR convention and matches Cook et al.'s sampling of the upper ice surface.
6. **Phase 3 bubbly-radius PDF:** a log-normal fitted to the radii retrieved by our own Phase 4 inversion from the 31 field spectra, so it is model-derived from data rather than directly measured.

Data that could not be obtained: the Stibal et al. (2017) cell counts and spectra. Their supplementary data host blocked automated download.
