# Phase 3: Monte Carlo uncertainty propagation and Sobol' sensitivity analysis

Phase 3 asks how micro-scale uncertainty (calibrated pigment spectrum, cell geometry, pigment loading) compares with environmental variability (ice structure, density, cell abundance, illumination) in controlling the macro-scale forcing of glacier algae.

## Quick start

```bash
git clone --depth 1 https://github.com/jmcook1186/biosnicar-py.git        # once, next to this repo
pip install numpy scipy pandas matplotlib seaborn pyyaml tqdm miepython SALib

python phase3/tests/test_phase3.py                                  # validation suite
python phase3/run_phase3.py --phase1-l2 phase1/results/level2       # full analysis (~5-10 min)
python phase3/run_phase3.py --demo --n-mc 400 --n-sobol 256         # quick pipeline check (DEMO)
python phase3/run_phase3.py --phase1-l2 ... --reuse                 # re-analyse / re-plot cached runs
```

## Files

| File | Purpose |
|---|---|
| `param_space.py` | Parameter PDFs (frozen `scipy.stats`), Latin Hypercube and Saltelli designs, parameter table |
| `forward_model.py` | Sample → tiers A–D albedo/forcing using the Phase 2 physics, with parallel batch evaluation |
| `qstar_table.py` | Packaging-factor lookup table (Q* vs optical depth and aspect ratio), < 0.3 % from direct chords |
| `stats_tools.py` | Descriptive statistics and bootstrap CIs, Sobol' wrapper (SALib), convergence, tornado analysis |
| `figures3.py` | Fig. 3A–3C and S3 (Phase 2 style, PDF and 600-dpi PNG) |
| `run_phase3.py` | Driver (all options: `-h`) |
| `tests/test_phase3.py` | Validation suite |

## Uncertainty model (`param_space.default_parameters`)

| Scale | Parameter | PDF |
|---|---|---|
| molecular | ΔE, band shift of the Phase 1 spectrum | posterior of the fit to the HPLC spectra of the isolated pigment (`tddft_calibration`); placeholder spectrum: N(+0.07, 0.03) eV |
| molecular | f, scale (oscillator strength × glucoside → phenol-equivalent units) | log-normal posterior of the fit to the measured S6 extract MAC; placeholder: median 38 |
| molecular | band FWHM | posterior of the HPLC-shape fit; placeholder: N(1.17, 0.02) eV |
| molecular | φ, Fe-complexed fraction (tier D) | posterior of the fit of the measured Fe-purpurogallin increment to the extract MAC; placeholder: N(0.27, 0.03) on [0, 1] |
| cellular | cell volume V | N(2320, 542) µm³, truncated at 400: per-sample biovolume per cell, 180 S6 samples (Williamson et al. 2020) |
| cellular | aspect L/d | U(1.46, 2.37): between the two species' measured means (Procházková et al. 2021) |
| cellular | intracellular phenolics c_i | N(22.0, 8.8) kg m⁻³, truncated at 0.5: phenolics per cell (53 samples) / S6 biovolume |
| cellular | size exponent γ of the concentration (c ∝ V^γ) | N(−0.73, 0.73): maximum-likelihood fit to 64 S6 samples, jackknife SE (γ = 0 is equal concentration) |
| environmental | ice specific surface area | ln SSA ~ N(−0.97, 0.35) m² kg⁻¹: 19 measurements (Cooper et al. 2021; Dadic et al. 2013), converted to BioSNICAR's bubble radius at 690 kg m⁻³ |
| environmental | surface density ρ | U(330, 560) kg m⁻³: measured weathering-crust range (Cooper et al. 2018) |
| environmental | cell abundance | log₁₀ B ~ N(3.56, 0.78): 180 S6 surface-ice counts |
| environmental | clear-sky transmissivity T | N(0.919, 0.036): 2301 clear-sky hours, PROMICE KAN_M |

**Every PDF comes from data** (`data/empirical/SOURCES.md`); `parameters.csv` lists each one with its rationale.
- The molecular PDFs are the marginal posteriors of the empirical calibration of your Phase 1 spectrum, so they change when you rerun with production output. The run writes them to `tables/tddft_calibration.json`.
- Their correlations are dropped, because Sobol' analysis needs independent inputs; the calibration summary reports them.

Fixed:
- SZA 45° (about solar noon at S6 in July).
- Sub-Arctic-summer spectrum.
- Cell asymmetry parameter from Mie theory with the measured refractive index (Phase 2).
- Bubbly ice over 690 kg m⁻³ ice, with a 2 cm algal layer.
- Chlorophyll a/b and carotenoids at their measured per-cell concentrations in tiers B–D. All inputs are treated as independent, as classical Sobol' analysis requires. If you later find correlated inputs (for example cell size with c_i), switch to Shapley effects.

## Outputs analysed

| Output | Meaning |
|---|---|
| `rf_A`–`rf_D` | Instantaneous SW forcing relative to clean ice (W m⁻²), 300–2500 nm |
| `bba_*` | Broadband albedo (300–2500 nm) |
| `eff_*` | Forcing per 10⁴ cells mL⁻¹. Removes the trivial scaling with abundance so micro-scale controls become visible |
| `d_CB` | RF_C − RF_B: the packaging effect |
| `d_DC` | RF_D − RF_C: the Fe-complexation and aggregation effect |
| `d_XA` | RF_X − RF_A: anomaly relative to BioSNICAR's empirical algae |

## Methods

### 1. Monte Carlo

- **Design:** Latin Hypercube, N = 1000 (`scipy.stats.qmc`, random-CD optimized). It maps to physical values through each PDF's inverse CDF.
- **Statistics:** mean, SD, median, IQR, the P2.5–P97.5 spread, and bootstrap 95% confidence intervals of the mean and median (in `mc_summary.csv`).
- **Wording for the paper:** the P2.5–P97.5 spread is an *uncertainty interval* of the forcing. The bootstrap intervals are *confidence intervals* of a statistic. Don't mix the two terms.

### 2. Sobol'

- **Design:** Saltelli (2010) with a scrambled Sobol' sequence and N = 1024, giving N(2D + 2) = 24,576 runs.
- **Estimators:** SALib, with first-order S₁, total-effect S_T and second-order S_ij indices, each with 1000-sample bootstrap 95% CIs.
- **Interactions:** S_T − S₁ is the share of variance due to interactions.

### 3. Sobol' by scale

- **Design:** a separate grouped Saltelli design (molecular, cellular and environmental as three groups), N(G + 2) = 5,120 runs (first-order and total-effect indices only).
- **What it answers:** "which *scale* drives the variance" directly, without summing parameter-level indices that overlap.

### 4. Supporting analyses

- **Convergence (Fig. S3):** S_T recomputed on Sobol'-sequence prefixes N = 64…1024. Report the indices only once these curves have flattened.
- **Tornado:** each parameter moved from P5 to P95 with all others at their medians (local, one-at-a-time). It complements the global Sobol' indices and does not replace them.

## Speed-ups (all validated against the Phase 2 reference)

- **Q\* lookup table:** < 0.3 % from direct chord averaging.
- **Radius snapping:** sampled radii snap to the nearest BioSNICAR lookup-table radius. The bubbly table has 10 µm steps up to 5 mm and 500 µm steps above.
- **Clean-ice optics:** granular-ice optics are cached per radius.
- **In-memory lookup tables:** BioSNICAR's tables are read fully into memory before workers fork. Its lazy `.npz` handles are not fork-safe, and its grain-table accessor decompresses the whole array on every call. Results are bit-identical (tested in Phase 2).

Each model run then takes about 17 ms. Forked workers inherit the warmed caches.

## Validation (`tests/test_phase3.py`, all passing)

| Check | Result |
|---|---|
| Ishigami benchmark (analytic S₁, S_T, S₁₃) | Recovered within 0.03 at N = 4096 |
| LHS stratification | Exactly one sample per stratum in every dimension |
| Q\* table vs direct Monte Carlo chords | < 0.5 % |
| Fast forward model vs Phase 2 reference (direct Q\*, fresh BioSNICAR runner) | < 0.5 % |
| Bootstrap 95 % CI coverage, log-normal test | ≥ 50/60 |

## Figures

| Figure | Content |
|---|---|
| Fig. 3A | Kernel-density PDFs of (a) forcing and (b) forcing per 10⁴ cells for tiers A–D, with median and 95 % interval bars |
| Fig. 3B | (a, b) Parameter-level S₁ (filled) and S_T (outline) with 95 % CIs, colored by scale; (c) scale-grouped indices for total and per-cell forcing |
| Fig. 3C | (a) Tornado for RF_D; (b) Monte Carlo binned mean per-cell forcing over cell volume × c_i; (c) pairwise S_ij among cell volume, aspect, c_i, ice radius, density and abundance (full matrix in `tables/sobol_S2_*.csv`) |
| Fig. S3 | Convergence of S_T with Saltelli sample size |

## Interpreting the results

- **Total forcing.** Abundance spans two orders of magnitude, so it dominates the variance of total forcing. With the empirical PDFs (log₁₀ B SD 0.78), S_T ≈ 0.9–1 in a small development run. That's physically expected, but it isn't the scientifically interesting result.
- **Per-cell forcing.** The per-cell outputs (`eff_*`) answer the micro-scale question. With abundance factored out, cell geometry and c_i, the parameters that control packaging, carry most of the variance.
- **Molecular uncertainty.** After the empirical calibration, the molecular parameters are pinned by measured spectra. They contribute S_T < 0.003 in the development run, against about 0.3 (cellular) and 0.7 (environmental) for per-cell forcing. The calibration removes most of the TD-DFT error rather than propagating it. Re-check this with your production Phase 1 spectrum.
- **Tier D.** Tier D adds the measured Fe-purpurogallin absorption at the fitted complexed fraction φ; its uncertainty is the φ posterior.
