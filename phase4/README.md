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

**All default inputs are empirical.** Every number and its source is listed in `data/empirical/SOURCES.md` and loaded through `phase2/empirical_data.py`. Every run first validates against 64 field samples from two sites (spectra plus cell counts; see "Field validation" below).

| Run | Time (4 CPU cores) |
|---|---|
| Synthetic, 3 × 3 km at 20 m (2.3 × 10⁴ pixels), field validation included | about 7 min |
| Real scene, 6 × 6 km at 40 m (2.2 × 10⁴ pixels), field validation included | about 4 min |

## Files

| File | Purpose |
|---|---|
| `s2_io.py` | STAC search (earth-search) or direct S3 item, windowed COG reads, reflectance scale/offset, SCL snow/ice mask, Copernicus GLO-30 DEM and slope, GeoTIFF export |
| `emulator.py` | Two-species BioSNICAR emulator (and Tier A); Sentinel-2 SRF convolution, albedo, forcing, pigment mass per node; PCHIP refinement |
| `priors.py` | Empirical priors: abundance, community fraction, ice radius, reflectance factor k (see table below) |
| `field_validation.py` | Leave-one-out validation on the Cook et al. (2020) field spectra and haemocytometer counts; picks the reflectance noise σ by maximum marginal likelihood |
| `inversion.py` | Exact grid-quadrature posterior with analytic k marginalization (all pixels); emcee MCMC with R-hat/ESS (selected pixels) |
| `synthetic.py` | Synthetic ground truth computed with direct BioSNICAR runs (avoids the inverse crime) |
| `empirical.py` | Band-ratio index baseline, fitted to the 59 counted field samples |
| `validation.py` | Bias, RMSE, MAE, R², 95 % coverage, error-reduction tables, field-point matching |
| `maps4.py` | Cartopy UTM maps with graticule, scale bar, north arrow; corner plots; validation scatter |
| `run_phase4.py` | Driver (`-h` for all options) |
| `bias_study.py` | Field-site bias variants (dust, film geometry, site calibration), scored by the same LOO validation |
| `dust_sensitivity.py` | Dust prior at the independent site: S6 prior vs broad, 10× cleaner and no dust |
| `satellite_validation.py` | Sentinel-2 pixels vs counted S6 2017 samples (pixel, 3 × 3, site mean; within-pixel variability) |
| `multi_scene.py` | Seasonal series of 2019 scenes (six dates); per-scene medians and daily-mean forcing scaled by measured PROMICE SW |

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

State z = [log₁₀ B, f_n, r, dust], plus a nuisance factor k.
- **Fractions:** f_a = 1 − f_n, so the requested [B, f_alaskanum, f_nordenskioeldii, grain] lies on the simplex.
- **Units:** B is in cells mL⁻¹ in the 2 cm surface layer (BioSNICAR's convention).
- **Species optics:** each species is a Phase 2 packaged cell (tier C, or tier D with `--tier D`).
  - **Geometry:** Greenland mean volumes (Halbach et al. 2022) with measured length:width ratios (Procházková et al. 2021):
    - A. nordenskioeldii: 10.75 × 25.44 µm.
    - A. alaskanum: 8.95 × 13.08 µm.
  - **Phenolic concentration:** measured phenolics per cell / measured biovolume per cell (Williamson et al. 2020, S6), scaled per species by the measured size dependence c ∝ V^γ (γ = −0.73 ± 0.73, 64 samples): 19.6 kg m⁻³ (A. nordenskioeldii) and 41.6 kg m⁻³ (A. alaskanum).
  - **Photosynthetic pigments:** cells also hold chl a, chl b and carotenoids at their measured per-cell concentrations, with Williamson et al.'s in vivo MACs (`--no-photosynthetic` removes them).
  - **Phenolic MAC:** one of
    - the Phase 1 TD-DFT spectrum **calibrated against measured spectra of the pigment** (`--phenol tddft`, the default; `phase2/tddft_calibration.py`);
    - the uncalibrated spectrum (`tddft_raw`);
    - the measured extract MAC (`--phenol williamson2020`).

    Tier D adds the measured Fe-purpurogallin absorption (Procházková et al. 2025) at the fitted complexed fraction.
  - **Scattering:** g(λ) from Mie theory with the measured cell refractive index (Phase 2).
- **Ice:** bubbly ice (granular ice cannot reach the field NIR reflectance). A 2 cm surface layer at 450 kg m⁻³ (weathering crust) lies over 690 kg m⁻³ (Cooper et al. 2018). The 2 cm is the sampling depth that defines the measured cells mL⁻¹. The bubble radius grid spans 0.3–20 mm.
- **Forcing:** SW↓ uses the clear-sky transmissivity measured at PROMICE KAN_M (0.919).

### Likelihood

R_b = k·F_b(z) + ε_b, with ε_b ~ N(0, σ²).
- **k** ~ N(0.898, 0.175²): the measured HCRF/albedo ratio of 51 field spectra, band-averaged with the ESA SRFs.
- **σ** is chosen by maximum marginal likelihood on the field samples: 0.020 for our model and 0.010 for Tier A. Pass `--sigma` to override it.

### Priors (all flags)

| Parameter | Prior | Source |
|---|---|---|
| log₁₀ B | N(3.559, 0.778) | 180 S6 surface-ice counts, 2016 (Williamson et al. 2020); independent of the 2017 validation samples |
| f_n | Beta(17.37, 11.42) | Moment fit to 3 Greenland surveys: 0.65, 0.66, 0.50 |
| r | Log-normal, estimated from the 64 field spectra by empirical Bayes (ours: median 600 µm, ln-SD 2.0) | Candidates: the measured bubbly-ice SSA prior (Cooper et al. 2021; Dadic et al. 2013; median 2.8 mm at 690 kg m⁻³) and a grid of log-normals. The field data prefer the wide distribution (log evidence 230 vs 183): the SSA measurements are of ice 0.1–1 m deep, not the surface crust |
| k | N(0.898, 0.175) | 51 field ARF spectra |
| dust | Log-normal, median 3.1 × 10⁵ ppb, ln-SD 0.46 (on by default; `--no-dust` removes it) | Moment-matched to the measured S6 loading, 342 µg g⁻¹ mean, relative SD 0.49 (Cook et al. 2020). The field spectra strongly prefer it (field-site bias study below) |

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
| Grid posterior vs MCMC (MCMC in ln r) | Means and SDs within 0.08 dex; R-hat < 1.05 |
| Posterior calibration (truth drawn from the empirical priors) | 95 % intervals cover the true log B in ≥ 90 % of 120 synthetic pixels |
| Sentinel-2 offset rule and DEM slope | Correct |

## Field validation (runs before every inversion)

**Data.** 64 co-located samples, each with a nadir HCRF spectrum and a haemocytometer count of the algae in the ice under the spectrometer (`data/empirical/SOURCES.md`):

| Dataset | Samples | Notes |
|---|---|---|
| S6, SW Greenland, 13–24 July 2017 (Cook et al. 2020 archive) | 46 (5 with 0 cells) | Counts from the primary workbook, with Poisson errors from the number of cells counted |
| Southern Greenland ice sheet, 5–6 Aug 2021 (Chevrollier et al. 2023) | 18 | An independent site, year, team and instrument; surface scraped 1–6 cm |
| S6, 17 June – 11 Aug 2014 (Stibal et al. 2017) | only if added by hand | Wiley/AGU block automated downloads (HTTP 403). Download the supporting information in a browser into `data/empirical/stibal_2017/` and fill in `mapping.json` (see that folder's README and `phase2/stibal2017.py`); the samples then enter the validation as `s6_2014` |

The earlier 31-sample version used biosnicar-py's copy of the S6 metadata. Two of its "0 cells" samples (22_7_SB6/7, the darkest misfits) were never counted: they are absent from the count workbook and are now excluded.

**Method.**
- Spectra are band-averaged with the ESA S2A SRFs. Each sample uses the emulator at its own solar zenith.
- All hyperparameters are chosen leave-one-out by maximum marginal likelihood of the other samples: the reflectance noise σ, and the ice-radius prior (the measured-SSA prior vs a grid of log-normals).
- The band-ratio regression is also leave-one-out.
- `coverage95_obs` widens the observation by its counting error.

**Results for the 59 counted samples** (log₁₀ cells mL⁻¹; `--phenol williamson2020`, with photosynthetic pigments):

| Method | Data | n | Bias | RMSE | Spearman ρ | 95 % coverage |
|---|---|---|---|---|---|---|
| **Physics-informed Bayesian (ours)** | all | 59 | **+0.38** | **0.57** | **0.79** | **0.78** |
| Tier A Bayesian | all | 59 | +0.72 | 0.84 | 0.77 | 0.37 |
| Empirical band-ratio (LOO) | all | 59 | 0.00 | 0.58 | 0.43 | – |
| **Ours** | S Greenland 2021 | 18 | **+0.17** | **0.41** | 0.51 | **0.89** |
| Tier A | S Greenland 2021 | 18 | +0.41 | 0.50 | 0.48 | 0.67 |
| Empirical band-ratio (LOO) | S Greenland 2021 | 18 | −0.49 | 0.64 | 0.10 | – |
| Ours | S6 2017 | 41 | +0.47 | 0.63 | 0.82 | 0.73 |
| Cook et al. (2020) published inversion | S6 2017 | 13 | −0.13 | 0.24 | 0.75 | – |

The Cook row covers only the 13 samples that have a published non-zero retrieval.

**What this shows**
- **Overall:** our model is the most accurate and best calibrated of the methods scored on all samples. At the independent site its bias is +0.17 dex with 89 % coverage. There the empirical regression, fitted mostly to S6, fails (bias −0.49, ρ 0.10), but the physics transfers.
- **Bias:** a positive bias remains (+0.47 dex at S6). The intervals are still somewhat too narrow (78 % coverage overall). Do not quote S6 abundances as calibrated.
- **TD-DFT path:** with the placeholder Phase 1 spectrum (`--phenol tddft`), the bias is +0.78 dex. The calibration fixes units and band positions, but a minimal-basis spectrum cannot reproduce the measured band shape (HPLC-shape R² 0.91). Rerun with production Level 2 output.

## Field-site bias study (`bias_study.py`)

Before this study, abundance at S6 was overestimated by +0.47 dex. Each candidate explanation was scored with the same leave-one-out field validation, with σ and the radius prior re-selected per variant (measured phenolic MAC, f_n step 0.2):

| Variant | S6 bias | S6 coverage | All 59: bias / RMSE / coverage | S Greenland: bias / coverage | Log evidence |
|---|---|---|---|---|---|
| Baseline (uniform 2 cm, no dust; meltwater-unit fix) | +0.49 | 0.73 | +0.40 / 0.59 / 0.78 | +0.19 / 0.89 | 227 |
| **+ measured S6 mineral dust** | **+0.17** | **0.81** | **+0.05 / 0.49 / 0.83** | −0.21 / 0.89 | **284** |
| Algae in the top 2 mm (same cells m⁻²) | +0.56 | 0.81 | +0.47 / 0.61 / 0.83 | +0.26 / 0.89 | 208 |
| Its control: same 3-layer column, uniform | +0.49 | 0.73 | +0.40 / 0.59 / 0.78 | +0.20 / 0.89 | 227 |
| Dust + film | +0.20 | 0.85 | +0.09 / 0.46 / 0.86 | −0.16 / 0.89 | 264 |

**Findings**
- **Unit audit.** Field counts are per mL of meltwater, while BioSNICAR's cells/mL are per mL of solid ice. The correction (× 0.917) is in the bridge and moves abundances by only 0.04 dex. The phenol-equivalent units cancel exactly between the measured MAC and the per-cell content (`data/empirical/SOURCES.md`, unit audit).
- **Sampling depth.** Cook et al.'s archived manuscript states that the S6 counts are of "the upper 2 cm of the ice surface", so the 2 cm layer matches the measurement. A thin surface film is disfavoured by the data (log evidence 208 vs 227): concentrating the cells changes darkening by only 4–8 %.
- **Mineral dust.** Without a dust term, all visible darkening is attributed to algae. With the measured S6 dust loading as prior, the S6 bias falls by two-thirds, the overall bias to +0.05 dex, and the evidence rises by 57. **Dust is therefore on by default.**
  - At southern Greenland the S6 dust prior is a transfer assumption, as no dust was measured there. The bias becomes −0.21 dex, and coverage stays 0.89.
- **Site calibration factor (rejected).** A multiplicative factor estimated at S6 (leave-one-out) removes the S6 bias by construction. Transferred to the independent site, it drives the bias there to −0.30 to −0.37 dex and coverage to 0.72. That is the test a satellite application faces, so no such factor is applied.

Run `python phase4/bias_study.py --phenol tddft|williamson2020` to reproduce; tables go to `results/bias_study*/`.

## Results so far: read before quoting

### 1. Synthetic truth generated with our optics (3 × 3 km, 2.3 × 10⁴ pixels; field-estimated priors, σ = 0.02)

| Method | log B bias | log B RMSE | Forcing bias | Forcing RMSE | 95 % coverage (log B / RF) |
|---|---|---|---|---|---|
| Empirical index (field-fitted) | +0.75 | 1.23 | – | – | – |
| Tier A Bayesian | −0.01 | 0.63 | −2.9 W m⁻² | 22.6 W m⁻² | 0.83 / 0.88 |
| **Physics-informed Bayesian** | +0.09 | **0.56** | **+1.0 W m⁻²** | **16.4 W m⁻²** | **0.93 / 0.99** |

Compared with Tier A, forcing RMSE is 27 % lower, albedo RMSE 30 % lower and log B RMSE 11 % lower. The MCMC cross-check samples ln r and starts from draws of the grid posterior. It converges (R-hat ≤ 1.04 at 6000 steps) and agrees with the grid means within 0.05 dex.

### 2. Community fraction is not identifiable

The posterior SD of f_n is 0.85–1.00 of its prior SD, so the f_n maps show the empirical prior. What *is* constrained is the total pigment mass per mL.

### 3. Real scene (S2A, 23 July 2019, 40 m, `--phenol williamson2020`)

- **Fit:** both models fit every pixel (0 % χ² failures) with bubbly ice (numbers below are from the run before dust became the default; they are regenerated by `scripts/reproduce_all.sh`).
- **Model evidence:** median ln Bayes factor −1.5, weakly favouring Tier A.
- **Forcing:** median algal forcing is 43 W m⁻² (ours) vs 83 W m⁻² (Tier A). Both are higher than in the previous revision because SW↓ now uses the measured clear-sky transmissivity (0.92 instead of 0.75).
- **Optics tension:** with the measured phenolic MAC at the measured concentration, modelled cells absorb almost all visible light incident on them, and BioSNICAR's measured in vivo cell absorption (Tier A) is lower towards 700 nm. This disagreement between two empirical sources is reported, not tuned away.

## Structural model error τ (calibrated intervals)

The posterior is conditional on the forward model, so it is too narrow when the model is wrong. `field_validation.fit_model_error` therefore fits a model-error term τ (dex, Gaussian in log₁₀ B):
- **Likelihood.** y ~ N(posterior mean, posterior SD² + τ² + Poisson count SD²), maximised over the counted samples, leave-one-out.
- **Interval.** The calibrated 95 % interval is mean ± 1.96 √(SD² + τ²).
- **Where it is used.** τ is fitted on the field samples and applied unchanged to the maps (`ours_log_b_q025_cal`/`q975_cal` bands, Fig. 4A panel e) and the satellite check.

With the measured phenolic MAC and dust on, **τ = 0.37 dex**. 95 % coverage:

| | Posterior only | With τ |
|---|---|---|
| S6 2017 (41) | 0.80 | 0.98 |
| S Greenland 2021 (18) | 0.89 | 1.00 |

The calibrated intervals are slightly conservative: the Gaussian τ also absorbs the remaining S6 bias of +0.17 dex.

## Satellite-scale validation (`satellite_validation.py`)

**Question.** The field validation is at plot scale (spectrometer footprint ≈ 0.5 m). Does the same retrieval, applied to Sentinel-2 L2A pixels with the field-calibrated σ, radius prior and τ, agree with the counts?

**Data.**
- **Samples:** 20 counted S6 samples with GPS positions (15 and 21–23 July 2017). Positions come from Tedstone et al. (2020), UK PDC, in **UTM 23N**; see `data/empirical/SOURCES.md`.
- **Primary scenes:** the nearest clear scene. That is 11 Jul for the 15 Jul samples (−4 d) and 21 Jul for the 21–23 Jul samples (0 to −2 d).
- **Sensitivity scenes:** the same-day scenes of 15 Jul (34 % tile cloud) and 23 Jul (thin cirrus nearby, 32 pixels masked).
- **Grid:** a 300 m window at native 10 m. Every pixel at the site is Sen2Cor class 11.

**Results** (measured phenolic MAC, dust on; `records/satellite_validation_williamson2020/`):

| Comparison (primary scenes) | n | Bias (dex) | RMSE | Spearman ρ | 95 % coverage (with τ) |
|---|---|---|---|---|---|
| Pixel containing the sample | 20 | −0.14 | 0.74 | −0.07 | 1.00 |
| 3 × 3-pixel mean | 20 | −0.19 | 0.75 | −0.17 | – |
| Same samples, ground spectra (plot scale) | 20 | +0.16 | 0.47 | 0.85 | – |
| **Site mean per day** (mean of samples vs mean of their pixels) | 4 days | **+0.04, −0.34, −0.38, +0.07** | – | – | count SE 0.24–0.42 |

**What this shows**
- **Representativeness limit.** Counts of samples that share one 10 m pixel differ with an SD of **0.83 dex** (11 df). No pixel retrieval can match individual 1 m² samples better than that, and the pixel RMSE (0.74) is at that floor. Per-sample rank agreement is therefore not expected at 10 m.
- **The pixels agree with the site.** The day-mean abundance under the samples is reproduced within 0.04–0.38 dex, comparable to the sampling error of the 5-sample day means (0.24–0.42 dex).
- **Not a loss of information.** Per-spectrum posterior SD is the same at both scales (pixel 0.55–0.65 dex; plot 0.62 dex). Pixel retrievals span a narrower range (log B 3.3–3.8) because a 10 m pixel averages clean and heavily colonised patches (pixel B2 0.46–0.63; ground plots 0.11–0.73).
- **Scene choice matters as much as model choice.**
  - The 23 Jul scene, with thin cirrus nearby, is 0.3–0.4 dex lower than the clear 21 Jul scene.
  - The same-day 15 Jul scene (cloudy tile) is 0.17 dex lower than the clear 11 Jul scene.
  - Hence clear-scene selection (tile cloud < 10 %) for the maps.

## Outputs (`phase4/results/`)

- `geotiff/phase4_maps.tif`: one band per map. Bands include our posterior mean, SD and 2.5/97.5 % quantiles of log B, f_n, r, pigment, albedo (mean, SD) and forcing (mean, SD), plus χ², k, the Tier A and empirical maps, ΔRF, the log Bayes factor, elevation, slope and melt stage.
- `tables/validation_metrics.csv`, `error_reduction_vs_tierA.csv`, `error_reduction_vs_empirical.csv`, `mcmc_diagnostics.csv`, `field_validation_*.csv`
- `figures/Fig4A_biomass_maps`, `Fig4B_albedo_forcing_maps`, `Fig4C_corner_pixel{1,2}`, `FigS4_validation_scatter` (each as `.pdf` and 600-dpi `.png`)
- `summary.json`, plus emulator caches in `cache/`.

## Known limitations (remaining stated assumptions are listed in `data/empirical/SOURCES.md`)

- **Species concentrations.** These use the measured size scaling (γ = −0.73 ± 0.73), which is weakly constrained and consistent with equal concentration.
- **Radius prior.** The surface ice-radius prior is estimated from the field spectra (empirical Bayes). It is conditional on this forward model, and the measured-SSA alternative is reported alongside it.
- **Residual bias.** With dust on (default), abundance is overestimated at S6 by +0.17 dex and underestimated at southern Greenland by −0.21 dex (that site borrows the S6 dust prior: `dust_sensitivity.py`). Calibrated intervals (τ = 0.37 dex) cover 98–100 % of samples.
- **Validation scale.** At pixel scale, Sentinel-2 retrievals agree with site-mean counts but cannot be checked sample by sample: plot-to-plot variability inside a 10 m pixel is 0.83 dex (`satellite_validation.py`). Satellite products should be read as pixel or area means, not as the abundance of individual patches.
- **Reflectance geometry.** Satellite HCRF is compared with modeled albedo, with anisotropy absorbed by k.
- **Pixel independence.** Pixels are independent apart from the optional pooling pass.
- **Fixed structure.**
  - Algae fill a uniform 2 cm layer, which is the sampling depth of the counts.
  - Surface density is fixed per run (`--rho`).
  - Both ice layers share one bubble radius.
