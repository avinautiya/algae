# Phase 2: Cellular optics, pigment packaging and BioSNICAR radiative forcing

Phase 2 connects the Phase 1 molecular spectra to ice-surface radiative transfer:

1. It starts from the Phase 1 TD-DFT MAC(λ).
2. A Duysens packaging correction inside Ancylonema-like cells gives MAC_vivo(λ).
3. That defines per-cell extinction, single-scattering albedo and asymmetry (ext / SSA / g) on BioSNICAR's 480-band grid.
4. BioSNICAR's adding-doubling solver then computes albedo for bare weathering-crust ice.
5. The outputs are broadband albedo (300–2500 nm) and instantaneous forcing for model tiers A–D.

## Quick start

```bash
git clone --depth 1 https://github.com/jmcook1186/biosnicar-py.git   # next to this repo, or set BIOSNICAR_PATH
pip install numpy scipy pandas matplotlib seaborn pyyaml tqdm miepython

python phase2/run_phase2.py --demo                       # pipeline check, outputs marked DEMO
python phase2/run_phase2.py \
    --phase1-l1 phase1/results/level1 --phase1-l2 phase1/results/level2   # real run, ~1 min
python phase2/tests/test_phase2.py                       # validation suite
```

BioSNICAR must be used from a **git checkout**: it locates its `data/` folder relative to the repository. A pip wheel alone will not find its optical-property files.

## Files

| File | Purpose |
|---|---|
| `pigment_packaging.py` | Duysens/Morel–Bricaud sphere formula, a Monte Carlo chord-length method for any convex cell (cylinders), `mac_vivo`, `mac_vivo_grid` |
| `cell_optics.py` | Phase 1 loaders, 480-band mapping, per-cell optics for tiers B/C/D (Mie asymmetry parameter), the tier A loader |
| `tddft_calibration.py` | Empirical calibration of the Phase 1 spectrum against measured spectra of the pigment; tier D (Fe-complexed) MAC from the measured Fe-purpurogallin spectrum; Fig. S0 |
| `empirical_data.py` | Every measured input (see `data/empirical/SOURCES.md`) |
| `biosnicar_bridge.py` | Locates BioSNICAR, injects in-memory impurities into its mixing step and solver, clear-sky SW↓, broadband albedo and forcing, lap.npz export |
| `figures.py` | Publication style; Fig. 2A–2C plus supplementary S1–S2 (PDF and 600-dpi PNG) |
| `run_phase2.py` | Command-line driver (all parameters are flags; run with `-h`) |
| `tests/test_phase2.py` | Validation suite (see below) |
| `phase2_colab.ipynb` | Colab notebook |

## Physics

### 1. Packaging (Duysens flattening)

- **Absorption coefficient inside the cell.** With pigment at concentration c_i, the intracellular absorption coefficient is a_i(λ) = MAC(λ) · c_i.
- **Absorption by one cell.** For a randomly oriented convex cell (ray optics, relative refractive index ≈ 1):

  σ_abs = ⟨A⟩ · ⟨1 − e^(−a_i l)⟩

  Here l is the chord length through the cell and ⟨A⟩ = S/4 is the mean projected area.
- **Packaging factor.** Because the mean chord is ⟨l⟩ = 4V/S, the packaging factor is

  Q\* = ⟨1 − e^(−a_i l)⟩ / (a_i ⟨l⟩), and MAC_vivo = Q\* · MAC.
- **Spheres:** this reduces to the closed-form Morel & Bricaud (1981) result, Q\* = 3Q_a(ρ′)/(2ρ′) with ρ′ = a_i·d.
- **Cylinders (Ancylonema):** Q\* is evaluated with 2×10⁵ isotropic, uniformly random chords.

### 2. Cell optics for BioSNICAR

- **Extinction:** σ_ext = 2⟨A⟩, the geometric-optics limit (size parameter 40–600).
- **Scattering:** σ_sca = σ_ext − σ_abs,packaged.
- **Tier B** keeps tier C's scattering and swaps in the unpackaged absorption, MAC·m_pig. So B vs C isolates packaging alone.
- **Cell water:** absorption by intracellular water (60% water by volume, k of water) is included and packaged the same way.
- **Asymmetry parameter:** g(λ) from Mie theory for the equal-volume sphere (`cell_optics.mie_g`). It uses the measured cell refractive index of 1.38 (Chevrollier et al. 2023), the ice host (n = 1.31) and k from the cell's own absorption, giving g ≈ 0.98–0.99. Measured g of green microalgae is > 0.95 (Pilon & Kandilian 2016). BioSNICAR's 0.96 is available as `g_mode="fixed"`, and the van Diedenhoven parameterization (`vd2014`) treats n as relative to air. The van Diedenhoven SSA is always written out as an independent cross-check of the packaging-based SSA.

### 3. Spectral window

- The molecular MAC is used for 350–800 nm (`--window`).
- Pigment absorption is zero above 800 nm. If the input MAC still has a tail at 800 nm, as the demo `ppg.csv` does, this produces a visible step in Fig. S1; Phase 1 Gaussian tails are much smaller there. Widen `--window` if your spectra need it.
- Below 350 nm the MAC is held constant (`--uv-mode hold`); the alternatives are `molecular` and `zero`.

### 4. Radiative transfer

- BioSNICAR v2 adding-doubling solver, direct beam, sub-Arctic-summer spectrum (`--incoming 3`).
- Two-layer column (defaults from measurements; `data/empirical/SOURCES.md`): a 2 cm algae-bearing **bubbly** weathering crust over 2 m of clean ice.
  - The crust's density is 330/450/560 kg m⁻³: Cooper et al. (2018) range and mean.
  - The ice below is 690 kg m⁻³ (Cooper et al. 2018).
  - The bubble-radius sweep is the 2.5/16/50/84/97.5 % quantiles of the **measured specific surface area** of bubbly ice (Cooper et al. 2021; Dadic et al. 2013; ln SSA ~ N(−0.97, 0.35)), converted at 690 kg m⁻³: about 1.45–5.65 mm. The reference is the median, 2.85 mm. (Cooper et al.'s 9.3–10.6 mm are ice-sphere radii, not bubble radii; only their SSA carries over.)
  - Granular ice (`--ice-mode grains`) cannot reproduce the field NIR reflectance (see Phase 4).
- **Cell counts are per mL of meltwater** (1 mL = 1 g), while BioSNICAR's input is per mL of solid ice (it divides by 917 kg m⁻³). The bridge converts (× 0.917, `biosnicar_bridge.MELTWATER_TO_BIOSNICAR`), so the column number of cells equals the measured count × ρ·dz.
- `IceSpec(film_dz=…, film_only=…)` puts the same cells per m² into a thin surface film (a 3-layer column), for tests of the vertical distribution. Compare at the same `film_dz`, because splitting a layer changes the delta-Eddington solution by about 0.005 albedo.
- Algae are in the top 2 cm, in cells mL⁻¹. This is the sampling depth that defines the measured abundances: "the top 2 cm collected" (Williamson et al. 2018); "scraping off the top ~2 cm" (Halbach et al. 2025).

### 5. Forcing

RF = SW↓ Σ_λ f(λ)[α_clean(λ) − α(λ)] over 300–2500 nm. Here f is BioSNICAR's normalized irradiance spectrum, and SW↓ = S₀ cos θ · T^(1/cos θ) for clear sky, with **T = 0.919** fitted to 2301 clear-sky hours of PROMICE KAN_M radiation (`empirical_data.clear_sky_transmissivity`; the former 0.75 gave about 35 % less SW↓); pass `--sw-down` to use measured fluxes instead (e.g. PROMICE). Broadband albedo uses the same weighting.

## Model tiers

| Tier | Absorption | Source |
|---|---|---|
| A | BioSNICAR empirical glacier algae | `ice_algae_empirical_Chevrollier2023` (whole-cell measured optics) |
| B | Calibrated Level 2 MAC, unpackaged | Phase 1 TD-DFT, calibrated (below) |
| C | Calibrated Level 2 MAC × Q\* | Phase 1 TD-DFT, calibrated + `pigment_packaging.py` |
| D | C + measured Fe-complex absorption × Q\* | Procházková et al. (2025) Fig. 4 at the complexed fraction fitted to the S6 extract (or `--level34-csv`) |

### Empirical calibration of the TD-DFT spectrum (`tddft_calibration.py`, Fig. S0)

The pigment of the field algae has measured spectra (Williamson et al. 2020 deposit). The Phase 1 spectrum is fitted to them in two stages:
1. **Band shift ΔE and Gaussian FWHM** from the diode-array spectra of the HPLC-isolated (uncomplexed) chromophore.
2. **Scale f and Fe-complexed fraction φ** from the measured MAC of whole S6 extracts. The fit is in log space, so the visible weighs as much as the UV peak.

Two details matter:
- **Units.** Phenolics were quantified as **phenol equivalents** (US EPA 420.1), in both the per-cell content and the extract MAC. f converts the TD-DFT MAC (per kg glucoside) to those units, which the uncalibrated model did not do. `--raw-tddft` turns the calibration off.
- **Tier D.** Tier D = f·[M + φ·I_M·D], where D is the measured Fe-induced absorbance change of purpurogallin (Procházková et al. 2025, Fig. 4; digitised).

With the placeholder Phase 1 spectrum (sto-3g):
- **Fit:** ΔE = +0.07 eV, FWHM = 1.17 eV, f = 38, φ = 0.27.
- **Effect:** calibrated tier C forcing is within about 6 % of BioSNICAR's measured tier A at the reference state (43 vs 41 W m⁻² at 10⁴ cells mL⁻¹; 232 vs 220 at 10⁵). Without calibration, tiers B–D combined a MAC per kg of glucoside with a concentration in phenol equivalents.

Rerun with the production Phase 1 output; `tables/tddft_calibration.json` records the fit and its diagnostics.

## Validation (`tests/test_phase2.py`, all passing)

| Check | Result |
|---|---|
| Monte Carlo chords vs analytic Duysens sphere | < 0.3 % |
| Cylinder mean chord vs Cauchy 4V/S | < 0.5 % |
| Exact Mie (miepython) vs Duysens, 10 µm cell, n_rel 1.03 | Shape < 2.5 %; magnitude +4–6 %, which is the n_rel² refractive enhancement that anomalous diffraction neglects |
| Energy conservation | σ_abs ≤ ⟨A⟩; 0 < SSA < 1; packaged ≤ unpackaged |
| Bridge vs BioSNICAR `run_model()` | Bit-identical albedo |

Two independent comparisons are printed or plotted on every run:
- **SSA cross-check:** the packaging-based SSA vs van Diedenhoven's (within about 0.07 for the demo input).
- **Packaging cross-check:** our Q\* vs BioSNICAR's measured glacier-algae packaging factor `pckg_GA.csv`, overlaid in Fig. 2A b.

## Outputs (`phase2/results/`)

- `tables/mac_vivo_grid.{npz,csv}`: MAC_vivo and Q\* over cell size × c_i × λ (task 1).
- `tables/cell_optics_480band.csv`: ext / abs / SSA / g for each tier on the 480-band grid.
- `tables/biosnicar_sweep.csv`: one row per run, with columns tier, grain_um, rho_top, sza, conc, bba (300–2500 nm), bba_biosnicar, rf, sw_down.
- `tables/summary_reference_state.csv`, `tables/spectral_albedo_reference.npz`
- `lap_entries/tier{B,C,D}_glacier_algae.npz`: our optics in BioSNICAR's `lap.npz` key format (see biosnicar-py `ADDING_DATA.md`).
- `figures/Fig2A_MAC_packaging`, `Fig2B_albedo_vs_concentration`, `Fig2C_forcing_anomaly_vs_grain`, `FigS1_cell_optics`, `FigS2_spectral_albedo` (each as `.pdf` and `.png`)
- `run_config.json`: every parameter, plus provenance of the inputs.

## Defaults to justify or replace in the paper

- **Cell size (empirical).** The reference cell is *A. nordenskioeldii*, 10.75 × 25.44 µm: Greenland mean volume (Halbach et al. 2022) with the measured length:width ratio (Procházková et al. 2021). Fig. 2A shows both species; *A. alaskanum* is 8.95 × 13.08 µm. `--sizes` overrides this.
- **Pigment concentration (empirical).** The phenolic concentration is c_i = 22.0 kg m⁻³: measured phenolics per cell / measured biovolume per cell at S6 (Williamson et al. 2020). The packaging grid uses c_i and c_i × (mean ± 1 SD)/mean.
  - **Stated assumption:** the concentration is the same in both species.
- **Photosynthetic pigments.** Chlorophyll a, chlorophyll b and carotenoids are included in tiers B–D at their measured per-cell concentrations, with in vivo MACs (Williamson et al. 2020). They absorb inside the same packaged cell. `--no-photosynthetic` gives the phenolic-only cell.
- **Species concentrations.** The size dependence of the measured intracellular concentration (`empirical_data.phenolic_size_scaling`, γ = −0.73 ± 0.73, consistent with 0) sets the per-species values: 19.6 (A. nordenskioeldii) and 41.6 kg m⁻³ (A. alaskanum). `--ref-species` selects the reference cell.
- **Remaining stated assumptions** (see `data/empirical/SOURCES.md`):
  - equal PG amounts in the two solutions of Procházková Fig. 4;
  - the equal-volume sphere for g;
  - uniform algae within the 2 cm layer.
- **Concentration range.** 10⁷ cells mL⁻¹ is far above observed blooms, which reach about 10⁴–10⁵ cells mL⁻¹.
- **Tier A cell size.** Tier A's extinction (7.1×10⁻¹⁰ m² per cell) reflects BioSNICAR's empirical cell size, which differs from our reference cell. Fig. S1 shows the per-cell absorption so the difference stays visible.
- **Forcing scope.** Forcing is instantaneous, clear-sky and direct-beam, with no melt feedbacks.
