# Phase 3: Monte Carlo uncertainty propagation and Sobol' sensitivity analysis

Phase 3 asks how micro-scale uncertainty (TD-DFT errors, cell geometry, pigment loading) compares with environmental variability (ice grains, density, cell abundance) in controlling the macro-scale forcing of glacier algae.

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
| molecular | ΔE, systematic TD-DFT excitation shift | N(0, 0.075 eV), truncated at ±0.15 eV |
| molecular | f scale, oscillator-strength error | log-normal, median 1, σ_ln 0.20 |
| molecular | band FWHM | U(0.25, 0.40) eV |
| molecular | LMCT ε_max (tier D) | U(3000, 5000) M⁻¹ cm⁻¹ per Fe |
| molecular | LMCT λ_max (tier D) | U(520, 620) nm |
| cellular | cell length L | U(10, 30) µm |
| cellular | cell diameter d | U(5, 15) µm |
| cellular | intracellular pigment c_i | log-uniform 10–200 kg m⁻³ |
| environmental | ice grain radius | U(1, 3) mm |
| environmental | surface density ρ | U(500, 800) kg m⁻³ |
| environmental | cell abundance | log-uniform 10³–10⁵ cells mL⁻¹ |

Fixed: SZA 55°, clear-sky SW↓ (Phase 2 parameterization), sub-Arctic-summer spectrum, 2 cm algal layer. The PDFs are priors, so state them in the paper and justify them; `parameters.csv` lists each one with its rationale. All inputs are treated as independent, as classical Sobol' analysis requires. If you later find correlated inputs (for example cell size with c_i), switch to Shapley effects.

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
- **Clean-ice optics:** cached per BioSNICAR lookup-table radius. BioSNICAR only tabulates radii in 20 µm steps between 1 and 5 mm, so sampled grain radii are snapped to the nearest one, at most 10 µm off.
- **In-memory lookup table:** an in-memory copy of BioSNICAR's grain table, because its accessor decompresses the whole array on every call. Bit-identical results, tested in Phase 2.

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
| Fig. 3C | (a) Tornado for RF_D; (b) Monte Carlo binned mean per-cell forcing over cell diameter × c_i; (c) pairwise S_ij among cell size, c_i, grain radius, density and abundance (full matrix in `tables/sobol_S2_*.csv`) |
| Fig. S3 | Convergence of S_T with Saltelli sample size |

## Interpreting the results

- **Total forcing.** Abundance spans two orders of magnitude, so it dominates the variance of total forcing (S_T ≈ 0.7–0.9 on the development input). That's physically expected, but it isn't the scientifically interesting result.
- **Per-cell forcing.** The per-cell outputs (`eff_*`) answer the micro-scale question. With abundance factored out, cell geometry and c_i, the parameters that control packaging, carry most of the variance.
- **Molecular uncertainty.** On the development input, TD-DFT uncertainty (ΔE, f, FWHM) contributes only a few percent. Re-check this with your real Phase 1 spectra: if a strong band sits near 400–500 nm, the energy shift moves it in and out of the solar peak, and the molecular share will rise.
- **Tier D.** Tier D still uses the provisional Fe(III)-phenolic surrogate from Phase 2, so its LMCT parameters stand in for real Level 3/4 uncertainty.
