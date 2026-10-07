# Phase 4: Physics-informed Bayesian inversion of Sentinel-2 imagery

Phase 4 retrieves glacier-algal abundance, *Ancylonema* community fraction and ice structure from Sentinel-2 L2A surface reflectance (B2 490, B3 560, B4 665, B8 842 nm). The forward model is the Phase 1–2 physics through BioSNICAR. The posterior is propagated into albedo and radiative-forcing maps and compared with standard empirical retrievals.

## Quick start (Colab: `phase4_colab.ipynb`)

```bash
git clone --depth 1 https://github.com/jmcook1186/biosnicar-py.git        # next to this repo
pip install numpy scipy pandas matplotlib seaborn pyyaml rasterio pyproj pystac-client cartopy emcee miepython

python phase4/tests/test_phase4.py
# 1) synthetic-truth experiment on the real Dark Zone grid (error reduction):
python phase4/run_phase4.py --source synthetic --phase1-l2 phase1/results/level2
python phase4/run_phase4.py --source synthetic --truth-model tierA --phase1-l2 ...   # reverse test
# 2) real Sentinel-2 scene (SW Greenland Dark Zone near PROMICE S6, 23 Jul 2019, 2.6 % cloud):
python phase4/run_phase4.py --source s2 --scene-id S2A_22WEV_20190723_0_L2A \
       --phase1-l2 ... --resolution 40 --f-step 0.2
# the measured in vivo phenolic MAC (Williamson et al. 2020) instead of the TD-DFT spectrum:
python phase4/run_phase4.py --source s2 ... --phenol williamson2020
```

**All default inputs are empirical.** Every number and its source is listed in `data/empirical/SOURCES.md` and loaded through `phase2/empirical_data.py`. Every run first validates against 31 field samples (spectra plus cell counts; see "Field validation" below).

| Run | Time (4 CPU cores) |
|---|---|
| Synthetic, 3 × 3 km at 20 m (2.3 × 10⁴ pixels), field validation included | about 4 min |
| Real scene, 6 × 6 km at 40 m (2.2 × 10⁴ pixels), field validation included | about 2.5 min (MCMC at 3000 steps) |

## Files

| File | Purpose |
|---|---|
| `s2_io.py` | STAC search (earth-search) or direct S3 item, windowed COG reads, reflectance scale/offset, SCL snow/ice mask, Copernicus GLO-30 DEM and slope, GeoTIFF export |
| `emulator.py` | Two-species BioSNICAR emulator (and Tier A); Sentinel-2 SRF convolution, albedo, forcing, pigment mass per node; PCHIP refinement |
| `priors.py` | Empirical priors: abundance, community fraction, ice radius, reflectance factor k (see table below) |
| `field_validation.py` | Leave-one-out validation on the Cook et al. (2020) field spectra and haemocytometer counts; picks the reflectance noise σ by maximum marginal likelihood |
| `inversion.py` | Exact grid-quadrature posterior with analytic k marginalization (all pixels); emcee MCMC with R-hat/ESS (selected pixels) |
| `synthetic.py` | Synthetic ground truth computed with direct BioSNICAR runs (avoids the inverse crime) |
| `empirical.py` | Band-ratio index baseline, fitted to the 25 counted field samples |
| `validation.py` | Bias, RMSE, MAE, R², 95 % coverage, error-reduction tables, field-point matching |
| `maps4.py` | Cartopy UTM maps with graticule, scale bar, north arrow; corner plots; validation scatter |
| `run_phase4.py` | Driver (`-h` for all options) |

## Data and preprocessing

- **Imagery.** Sentinel-2 **L2A** (Sen2Cor bottom-of-atmosphere surface reflectance) Cloud-Optimised GeoTIFFs from the public AWS archive. Reads are windowed over HTTP, so only the area of interest is downloaded.
- **Reflectance scaling.** Reflectance = DN × 10⁻⁴ + offset. The offset of −0.1 is applied for processing baseline ≥ 04.00 (from 25 January 2022).
- **Mask.** The scene classification layer keeps only class 11 (snow/ice).
- **Level-1C input.** Correct it with ESA Sen2Cor first. Sen2Cor's aerosol retrieval over bright ice is weakly constrained; the multiplicative factor k and the per-band σ absorb errors of a few percent.
- **Default scene.** `S2A_22WEV_20190723_0_L2A` (SZA 47°), 6 × 6 km around 67.08 °N, 49.35 °W.
- **Spectral response.** ESA's official Sentinel-2A/B/C SRFs (TN-15-0007 v4.0), averaged into BioSNICAR's 10 nm bins. The spacecraft is read from the scene ID.
- **Elevation.** Copernicus GLO-30 DEM, reprojected to the scene grid.
- **Melt diagnostics (mapped only).** Positive degree days from PROMICE KAN_L/M/U daily air temperature, interpolated in elevation day by day. These are available for 2019 scenes. They do not enter the prior: no published calibration links melt to abundance at the pixel scale.
- **Searching.** Use `--search START END` (STAC API) to find other dates or glaciers, or `--source local` for your own GeoTIFFs.

## Model

### State and fixed choices

State z = [log₁₀ B, f_n, r (, dust)], plus a nuisance factor k.
- **Fractions:** f_a = 1 − f_n, so the requested [B, f_alaskanum, f_nordenskioeldii, grain] lies on the simplex.
- **Units:** B is in cells mL⁻¹ in the 2 cm surface layer (BioSNICAR's convention).
- **Species optics:** each species is a Phase 2 packaged cell (tier C, or tier D with `--tier D`).
  - **Geometry:** Greenland mean volumes (Chevrollier et al. 2022) with measured length:width ratios (Procházková et al. 2021):
    - A. nordenskioeldii: 10.75 × 25.44 µm.
    - A. alaskanum: 8.95 × 13.08 µm.
  - **Phenolic concentration:** 22.0 kg m⁻³ = measured phenolics per cell / measured biovolume per cell (Williamson et al. 2020, S6).
  - **Photosynthetic pigments:** cells also hold chl a, chl b and carotenoids at their measured per-cell concentrations, with Williamson et al.'s in vivo MACs (`--no-photosynthetic` removes them).
  - **Phenolic MAC:** the Phase 1 TD-DFT spectrum (`--phenol tddft`, the default) or the measured in vivo phenolic MAC (`--phenol williamson2020`).
- **Ice:** bubbly ice (granular ice cannot reach the field NIR reflectance). A 2 cm surface layer at 450 kg m⁻³ (weathering crust) lies over 690 kg m⁻³ (Cooper et al. 2018). The optical radius is 0.3–20 mm.

### Likelihood

R_b = k·F_b(z) + ε_b, with ε_b ~ N(0, σ²).
- **k** ~ N(0.898, 0.175²): the measured HCRF/albedo ratio of 51 field spectra, band-averaged with the ESA SRFs.
- **σ** is chosen by maximum marginal likelihood on the field samples: 0.020 for our model and 0.010 for Tier A. Pass `--sigma` to override it.

### Priors (all flags)

| Parameter | Prior | Source |
|---|---|---|
| log₁₀ B | N(3.559, 0.778) | 180 S6 surface-ice counts, 2016 (Williamson et al. 2020); independent of the 2017 validation samples |
| f_n | Beta(17.37, 11.42) | Moment fit to 3 Greenland surveys: 0.65, 0.66, 0.50 |
| r | Log-uniform on [300, 20 000] µm | Range spanned by the field NIR at the measured densities. Not a measured distribution; Cooper et al. (2021) measured about 9.3–10.6 mm |
| k | N(0.898, 0.175) | 51 field ARF spectra |
| dust (`--dust` only) | Log-uniform over its nodes | Off by default: S6 dust was found weakly absorbing |

Use `--prior-scale` to test sensitivity to the priors. `--spatial-pooling` adds an empirical-Bayes second pass that pools neighboring pixels.

### Inference

- **Grid quadrature.** With three or four state dimensions, the posterior is evaluated *exactly* on the emulator grid for every pixel. k is integrated analytically using the matrix determinant lemma and Sherman–Morrison.
- **Outputs.** Posterior mean, SD and 2.5/50/97.5 % quantiles; MAP; posterior mean ± SD of pigment mass, albedo and forcing; per-pixel log evidence (Bayes-factor maps between models); and a χ² goodness-of-fit at the MAP.
- **MCMC cross-check.** emcee (differential-evolution moves) on the interpolated emulator with explicit k, for two representative pixels. It reports split R-hat, ESS and autocorrelation time, and overlays the MCMC samples on the exact grid marginals (Fig. 4C).

## Validation (`tests/test_phase4.py`, all passing)

| Check | Result |
|---|---|
| Analytic k-marginalized likelihood vs brute-force integration | Agrees, for k ~ N(1, 0.1) and for the empirical N(0.898, 0.175) |
| Emulator vs direct BioSNICAR at off-grid states (default empirical configuration) | < 0.002 reflectance |
| Grid posterior vs MCMC | Means and SDs within 0.08 dex; R-hat < 1.05 |
| Posterior calibration (truth drawn from the empirical priors) | 95 % intervals cover the true log B in ≥ 90 % of 120 synthetic pixels |
| Sentinel-2 offset rule and DEM slope | Correct |

## Field validation (runs before every inversion)

**Data.** Cook et al. (2020), S6, July 2017: 31 surface samples, each with:
- a nadir HCRF spectrum (ASD FieldSpec, clear sky, solar noon ± 2 h);
- a haemocytometer count of the cells in the ice inside the spectrometer footprint.

**Method.**
- Spectra are band-averaged with the ESA S2A SRFs, and each sample uses the emulator at its own noon SZA.
- Everything is out of sample: σ is chosen leave-one-out, and the band-ratio regression is also leave-one-out.
- 6 samples have zero counts (below the 62.5 cells mL⁻¹ counting step) and are reported separately.

**Results for the 25 counted samples** (log₁₀ cells mL⁻¹; `--phenol williamson2020`, with photosynthetic pigments):

| Method | Bias | RMSE | Spearman ρ | 95 % coverage |
|---|---|---|---|---|
| **Physics-informed Bayesian (ours)** | **+0.36** | **0.47** | **0.96** | 0.76 |
| Tier A Bayesian (BioSNICAR empirical algae) | +0.69 | 0.75 | 0.94 | 0.40 |
| Empirical band-ratio regression (LOO) | −0.01 | 0.51 | 0.74 | – |
| Cook et al. (2020) published inversion | −0.13 | 0.24 | 0.75 | – |

The Cook row is scored on only 13 samples: 12 of the 25 counted samples were retrieved as 0 cells and cannot be scored on a log scale.

**What this shows**
- **Ranking:** our model ranks samples best (ρ = 0.96) and beats Tier A on bias, RMSE and coverage.
- **Bias:** it still **overestimates abundance by 0.36 dex** (about 2.3×), and the 95 % intervals cover only 76 %. Do not quote it as calibrated.
- **Zero-count samples:** the two dark ones (22_7_SB6/7) are retrieved at about 10⁵ cells mL⁻¹ with large radii. Their darkness is not explained by counted algae.
- **TD-DFT default:** with `--phenol tddft`, the results depend on the Phase 1 spectrum. Rerun with production Level 2 output.

## Results so far: read before quoting

### 1. Synthetic truth generated with our optics (3 × 3 km, 2.3 × 10⁴ pixels; empirical priors, σ = 0.02)

| Method | log B bias | log B RMSE | Forcing bias | Forcing RMSE | 95 % coverage (log B / RF) |
|---|---|---|---|---|---|
| Empirical index (field-fitted) | −0.26 | 0.84 | – | – | – |
| Tier A Bayesian | +0.06 | 0.69 | −2.5 W m⁻² | 16.3 W m⁻² | 0.83 / 0.86 |
| **Physics-informed Bayesian** | +0.19 | **0.64** | **+1.4 W m⁻²** | **10.9 W m⁻²** | **0.93 / 1.00** |

- Forcing RMSE is 33 % lower than Tier A and albedo RMSE 39 % lower.
- With the wide empirical abundance prior, the log B RMSE gain is small (6 %).

### 2. Community fraction is not identifiable

The posterior SD of f_n is 0.85–1.00 of its prior SD, so the f_n maps show the empirical prior. What *is* constrained is the total pigment mass per mL.

### 3. Real scene (S2A, 23 July 2019, 40 m, `--phenol williamson2020`)

- **Fit:** both models fit every pixel (0 % χ² failures) with bubbly ice and no dust axis.
- **Model evidence:** median ln Bayes factor −1.6, weakly favouring Tier A.
- **Forcing:** median algal forcing is 36 W m⁻² (ours) vs 58 W m⁻² (Tier A).
- **Optics tension:** with the measured phenolic MAC at the measured concentration, modelled cells absorb almost all visible light incident on them. BioSNICAR's measured in vivo cell absorption (Tier A) is lower towards 700 nm. This disagreement between two empirical sources is reported, not tuned away.

## Outputs (`phase4/results/`)

- `geotiff/phase4_maps.tif`: one band per map. Bands include our posterior mean, SD and 2.5/97.5 % quantiles of log B, f_n, r, pigment, albedo (mean, SD) and forcing (mean, SD), plus χ², k, the Tier A and empirical maps, ΔRF, the log Bayes factor, elevation, slope and melt stage.
- `tables/validation_metrics.csv`, `error_reduction_vs_tierA.csv`, `error_reduction_vs_empirical.csv`, `mcmc_diagnostics.csv`, `field_validation_*.csv`
- `figures/Fig4A_biomass_maps`, `Fig4B_albedo_forcing_maps`, `Fig4C_corner_pixel{1,2}`, `FigS4_validation_scatter` (each as `.pdf` and 600-dpi `.png`)
- `summary.json`, plus emulator caches in `cache/`.

## Known limitations (non-empirical items are listed in `data/empirical/SOURCES.md`)

- **Equal concentration.** Intracellular pigment concentration is assumed equal in both species.
- **Radius prior.** The ice-radius prior bounds span the data; they are not a measured distribution.
- **Validation scale.** Field validation is at plot scale (spectrometer footprint), not at the 10–20 m pixel scale.
- **Reflectance geometry.** Satellite HCRF is compared with modeled albedo, with anisotropy absorbed by k.
- **Pixel independence.** Pixels are independent apart from the optional pooling pass.
- **Fixed structure.** Algae are confined to a 2 cm layer, and surface density is fixed per run (`--rho`).
