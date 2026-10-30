# Independent optical-variability data: what exists, and what it can support

Searched: the Cook 2020 data archive (Zenodo 10.5281/zenodo.3564501, downloaded and inspected), the Dark Zone paper text and supplement (Sci. Adv. 12, eady9482) and archive (Zenodo 18826013), and the SNICAR-fx particle files. No literature search beyond these sources was completed in this milestone; the absence of a dataset below means "not found here", not "does not exist".

## Inventory

| Item | Organism / material, surface | What was measured | Normalisation | Conditions | Uncertainty given | Variation across samples available? |
|---|---|---|---|---|---|---|
| Ice-algal pigment content (Cook archive `Cell_Counts/Pigment_concentrations.xlsx`, sheet `pigment_concs_per_cell`) | glacier ice algae (S6, 2016-2017 campaigns; species not stated in the sheet) | per-cell pigment amounts: chl a 0.0040 +- 0.0025 ng/cell (CV 64 %), carotenoids, phenolic total 4.1 % dry weight | ng per cell; % dry weight | extraction method not in the sheet (not read) | sd and se over samples; n not stated in the sheet | summary statistics only; **no per-sample values** in the sheet |
| Phenolic mass extinction coefficient (same file, sheet `mass_extinction_coefs`) | glacier-algal phenolic pigment | regression slope vs concentration, r2 0.82-0.84 at 250-260 nm, with error | m2 per mg | not read | slope err, r2 | one spectrum with regression error; not sample-resolved |
| Cell dimensions (Cook archive `Cell_Counts/Cell_dimensions.xlsx`) | Mesotaenium and Ancylonema, per plot | width, length, biovolume per sample | biovolume 290-1680 um3 across the 40 count-bearing plots (factor 5.8; CV 0.32) | microscopy | not stated | yes, per plot |
| Ice-algal in vivo absorption (SNICAR-fx `ice_algae_1nm_resolution.nc`; BioSNICAR `ice_algae_empirical_Chevrollier2023`) | glacier ice algae, S21 (and for the 2023 spectrum SW Greenland) | in vivo absorption cross-section; scattering neglected; no absorption > 830 nm assumed | per cell | not read | none in the file | **one spectrum** |
| Greenland mineral dust (Dark Zone paper; archive `mineral_dust_1nm_resolution.nc`) | mineral dust from snow, cryoconite holes and surface ice, 3 sites | filter-pad absorption of H2O2-cleaned, freeze-dried samples (n = 10); XRD mineralogy (n = 17, Table S1) | m2 per kg inorganic dry mass | resuspension, filter-pad with half null-point correction at 750 nm | **per-sample spectra are not in the archive**; only the mean file | mineralogy yes (hornblende 0.7-20 wt %, K-feldspar 4-36 wt %, hematite <= 0.1 wt %); absorption reported as similar |
| Greenland dust reflectance (Cook archive `Mineral_Dust_Data/Mineral_dust_reflectance.csv`) | dust from two samples (GRIS1617, GRIS1628) | integrating-sphere reflectance 350-2500 nm | reflectance factor | sample preparation / grain size dependent | none | **two spectra**, differ by about 0.07 at 350 nm |
| Cryoconite optics (Dark Zone archive `dark_cryoconite_1nm_resolution.nc`) | cryoconite, S21 | optics inverted from two reflectance measurements of optically thick layers (diffuse); g fixed at 0.8 | per mass | outdoor/diffuse | none | **two measurements, one site**, then merged |
| Snow-algae optics (`snow_algae.nc`, Chevrollier 2025) | red snow algae | in vivo absorption | per cell | snow | - | not substitutable for glacier ice algae (different pigment assemblage, medium, environment); not used |

## Is the variation chemical, or instrumental / preparation?

* **Pigment per cell:** the CV of 64 % for chl a is real between-sample variation in the measured quantity, but the sheet does not allow separating biology from extraction and counting error, and it is not paired with spectra. A per-cell absorption change cannot be derived from it without a pigment-to-absorption model that this repository does not have (and that would need its own validation).
* **Dust reflectance pair (GRIS1617 vs GRIS1628):** two powders measured in an integrating sphere; a 0.07 reflectance offset is as consistent with grain size, packing or sample mass as with mineralogy. **Two spectra do not establish chemical variability.** The matching mineral abundances in `Mineral_Abundances.csv` could be used to test whether the offset tracks Fe-bearing phases; this was not done.
* **Dust absorption in the Dark Zone data:** the authors report similar mass absorption across 10 samples despite varying mineralogy, so for this region the measured chemical variability maps to little optical variability. That argues against a dust-chemistry extension for Greenland.
* **Cryoconite:** the two merged measurements cannot show variation; the competitor itself names cryoconite optics as the least constrained input.

## Exploratory check (development data, retrospective)
If per-cell absorption scaled with measured cell biovolume, the Baseline A cell-count error would correlate with biovolume. Spearman correlation of (retrieved - counted) log10 cells with log10 biovolume across 40 count-bearing plots: 0.04 (p = 0.80). No support; not a proof of absence given the crude biovolume weighting of the two species (their relative proportions were not available).

## What is missing
Paired, sample-resolved measurements: absorption (or reflectance of an optically thick layer) **and** pigment composition for the same glacier-algal samples, repeated across at least several surfaces and sites with preparation replicates, to separate chemistry from preparation. Equivalently for dust: per-sample absorption spectra with mineralogy, which the Dark Zone authors measured (n = 10) but did not archive in the files we obtained. Contacting the authors for those per-sample spectra is the cheapest next step; none of this requires a computational chemistry calculation.
