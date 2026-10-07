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
       --phase1-l2 ... --ice-mode bubbly --dust --resolution 40 --f-step 0.2
```

| Run | Time (4 CPU cores) |
|---|---|
| Synthetic, 2 × 2 km at 20 m (10⁴ pixels) | about 1.5 min |
| Real scene, 6 × 6 km at 40 m with dust (2.2 × 10⁴ pixels) | about 3.5 min |

## Files

| File | Purpose |
|---|---|
| `s2_io.py` | STAC search (earth-search) or direct S3 item, windowed COG reads, reflectance scale/offset, SCL snow/ice mask, Copernicus GLO-30 DEM and slope, GeoTIFF export |
| `emulator.py` | Two-species BioSNICAR emulator (and Tier A); Sentinel-2 SRF convolution, albedo, forcing, pigment mass per node; PCHIP refinement |
| `priors.py` | Environmental priors: melt stage from lapse-rate degree-days, slope, community-ratio Beta prior, grain prior, reflectance factor k |
| `inversion.py` | Exact grid-quadrature posterior with analytic k marginalization (all pixels); emcee MCMC with R-hat/ESS (selected pixels) |
| `synthetic.py` | Synthetic ground truth computed with direct BioSNICAR runs (avoids the inverse crime) |
| `empirical.py` | Band-ratio index baseline (calibrated on Tier A simulations, or your field coefficients) |
| `validation.py` | Bias, RMSE, MAE, R², 95 % coverage, error-reduction tables, field-point matching |
| `maps4.py` | Cartopy UTM maps with graticule, scale bar, north arrow; corner plots; validation scatter |
| `run_phase4.py` | Driver (`-h` for all options) |

## Data and preprocessing

- **Imagery.** Sentinel-2 **L2A** (Sen2Cor bottom-of-atmosphere surface reflectance) Cloud-Optimised GeoTIFFs from the public AWS archive. Reads are windowed over HTTP, so only the area of interest is downloaded.
- **Reflectance scaling.** Reflectance = DN × 10⁻⁴ + offset. The offset of −0.1 is applied for processing baseline ≥ 04.00 (from 25 January 2022).
- **Mask.** The scene classification layer keeps only class 11 (snow/ice).
- **Level-1C input.** Correct it with ESA Sen2Cor first. Sen2Cor's aerosol retrieval over bright ice is weakly constrained; the multiplicative factor k and the per-band σ absorb errors of a few percent.
- **Default scene.** `S2A_22WEV_20190723_0_L2A` (SZA 47°), 6 × 6 km around 67.08 °N, 49.35 °W.
- **Elevation.** Copernicus GLO-30 DEM, reprojected to the scene grid.
- **Melt stage.** Positive degree-days from a station temperature record (`--tref-csv` with columns date and t_air, e.g. PROMICE S6 or KAN_M), lapsed at 6.5 K km⁻¹. Without a record, a documented climatological curve is used.
- **Searching.** Use `--search START END` (STAC API) to find other dates or glaciers, or `--source local` for your own GeoTIFFs.

## Model

### State and fixed choices

State z = [log₁₀ B, f_n, r (, dust)], plus a nuisance factor k.
- **Fractions:** f_a = 1 − f_n, so the requested [B, f_alaskanum, f_nordenskioeldii, grain] lies on the simplex.
- **Units:** B is in cells mL⁻¹ in the 2 cm surface layer (BioSNICAR's convention).
- **Species optics:** each species is a Phase 2 packaged cell (tier C, or tier D with `--tier D`) with its own geometry (`emulator.DEFAULT_SPECIES`). The defaults are placeholders, so set them from microscopy.

### Likelihood

R_b = k·F_b(z) + ε_b, with ε_b ~ N(0, σ_b²), σ = 0.02–0.025 (sensor + atmosphere + model), and k ~ N(1, 0.1²) for illumination, HCRF-versus-albedo and residual atmosphere.

### Priors (all flags)

| Parameter | Prior |
|---|---|
| log₁₀ B | N(3.8 + 0.8(s − ½) − 0.3·slope/10°, 0.8), where s is the melt stage |
| f_n | Beta(2, 2) |
| r | N(1200 + 1400·s, 600) µm |
| dust | Log-uniform over its nodes |

Use `--prior-scale` to test sensitivity to the priors. `--spatial-pooling` adds an empirical-Bayes second pass that pools neighboring pixels.

### Inference

- **Grid quadrature.** With three or four state dimensions, the posterior is evaluated *exactly* on the emulator grid for every pixel. k is integrated analytically using the matrix determinant lemma and Sherman–Morrison.
- **Outputs.** Posterior mean, SD and 2.5/50/97.5 % quantiles; MAP; posterior mean ± SD of pigment mass, albedo and forcing; per-pixel log evidence (Bayes-factor maps between models); and a χ² goodness-of-fit at the MAP.
- **MCMC cross-check.** emcee (differential-evolution moves) on the interpolated emulator with explicit k, for two representative pixels. It reports split R-hat, ESS and autocorrelation time, and overlays the MCMC samples on the exact grid marginals (Fig. 4C).

## Validation (`tests/test_phase4.py`, all passing)

| Check | Result |
|---|---|
| Analytic k-marginalized likelihood vs brute-force integration | Agrees |
| Emulator vs direct BioSNICAR at off-grid states | < 0.002 reflectance (measured 0.0004) |
| Grid posterior vs MCMC | Means and SDs within 0.08 dex; R-hat < 1.05 |
| Posterior calibration | 95 % intervals cover the true log B in ≥ 90 % of 120 synthetic pixels |
| Sentinel-2 offset rule and DEM slope | Correct |

## Results so far: read before quoting

These were obtained with a placeholder minimal-basis Phase 1 spectrum. Rerun with your production Level 2 spectra.

### 1. Synthetic truth generated with our optics (2 × 2 km, 10⁴ pixels)

| Method | log B bias | log B RMSE | Forcing bias | Forcing RMSE |
|---|---|---|---|---|
| Empirical index | −0.19 | 0.41 | – | – |
| Tier A Bayesian | −0.13 | 0.44 | +7.2 W m⁻² | 7.9 W m⁻² |
| **Physics-informed Bayesian** | **+0.01** | **0.34** | **+0.7 W m⁻²** | **1.9 W m⁻²** |

That is a 93 % lower |bias| and 23 % lower RMSE in log B than Tier A, and 76 % lower forcing RMSE.

### 2. Reverse test, truth generated with Tier A optics

- Our model now underestimates forcing by 9.5 W m⁻², and its forcing intervals cover only 56 % of truths.
- Biomass retrieval stays robust in both directions (bias ≤ 0.04 dex).
- **The forcing advantage is conditional on our optics being correct.** The four bands alone cannot decide between the two models: synthetic Bayes factors are ≈ 0 in both directions.

### 3. Community fraction is not identifiable

The posterior SD of f_n is 0.92–0.98 of its prior SD. Four broad bands cannot separate two species with the same pigment and similar cells. The maps show the prior, and Fig. 4A f says so on the panel.

What *is* constrained is the total pigment mass per mL. In the corner plots the B–f trade-off collapses onto it.

### 4. Real scene (23 July 2019)

- **Ice structure.** In granular mode the observed NIR (≈ 0.30) is unreachable even with 5 mm grains, so the real-data run uses solid bubbly ice plus mineral dust. Both models then fit (χ² failures < 1 %).
- **Model evidence.** The data **favor BioSNICAR's empirical algae** (median ln Bayes factor −13). Physically, real cells also hold chlorophyll and carotenoids that absorb across the visible, while our model represents only the phenolic pigment.
- **Next step.** Add photosynthetic pigments to tier C/D cells. The phenolic-only forcing (median 17 W m⁻²) should be read as the *phenolic contribution*, not total algal forcing; Tier A gives 96 W m⁻².

### 5. Field validation decides the question

Pass cell counts with `--field-csv` (columns lon, lat, cells_per_ml). The tables then report bias and RMSE against observations for all three methods.

## Outputs (`phase4/results/`)

- `geotiff/phase4_maps.tif`: one band per map. Bands include our posterior mean, SD and 2.5/97.5 % quantiles of log B, f_n, r, pigment, albedo (mean, SD) and forcing (mean, SD), plus χ², k, the Tier A and empirical maps, ΔRF, the log Bayes factor, elevation, slope and melt stage.
- `tables/validation_metrics.csv`, `error_reduction_vs_tierA.csv`, `error_reduction_vs_empirical.csv`, `mcmc_diagnostics.csv`, `field_validation_*.csv`
- `figures/Fig4A_biomass_maps`, `Fig4B_albedo_forcing_maps`, `Fig4C_corner_pixel{1,2}`, `FigS4_validation_scatter` (each as `.pdf` and 600-dpi `.png`)
- `summary.json`, plus emulator caches in `cache/`.

## Known limitations

- **Spectral response.** BioSNICAR's Sentinel-2 responses are tophat approximations; replace `data/band_srfs/sentinel2_msi.csv` with ESA's official ones for final numbers.
- **Reflectance geometry.** Satellite HCRF is compared with modeled albedo, with anisotropy absorbed by k.
- **Pixel independence.** Pixels are independent apart from the optional pooling pass.
- **Fixed structure.** Algae are confined to a 2 cm layer, and surface density is fixed per run (`--rho`).
