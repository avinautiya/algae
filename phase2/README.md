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
| `cell_optics.py` | Phase 1 loaders, 480-band mapping, per-cell optics for tiers B/C/D, the provisional tier D Fe(III)-phenolic model, the tier A loader |
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
- **Asymmetry parameter:** g = 0.96 for every tier, which is BioSNICAR's empirical value. `--g-mode vd2014` uses the van Diedenhoven parameterization instead. The van Diedenhoven SSA is always written out as an independent cross-check of the packaging-based SSA.

### 3. Spectral window

- The molecular MAC is used for 350–800 nm (`--window`).
- Pigment absorption is zero above 800 nm. If the input MAC still has a tail at 800 nm, as the demo `ppg.csv` does, this produces a visible step in Fig. S1; Phase 1 Gaussian tails are much smaller there. Widen `--window` if your spectra need it.
- Below 350 nm the MAC is held constant (`--uv-mode hold`); the alternatives are `molecular` and `zero`.

### 4. Radiative transfer

- BioSNICAR v2 adding-doubling solver, direct beam, sub-Arctic-summer spectrum (`--incoming 3`).
- Two-layer column (defaults from measurements; `data/empirical/SOURCES.md`): a 2 cm algae-bearing **bubbly** weathering crust over 2 m of clean ice.
  - The crust's density is 330/450/560 kg m⁻³: Cooper et al. (2018) range and mean.
  - The ice below is 690 kg m⁻³ (Cooper et al. 2018).
  - The optical radius sweep is 1–15 mm; the reference of 10 mm follows Cooper et al. (2021), who measured about 9.3–10.6 mm.
  - Granular ice (`--ice-mode grains`) cannot reproduce the field NIR reflectance (see Phase 4).
- Algae are in the top layer, in cells mL⁻¹, using BioSNICAR's convention.

### 5. Forcing

RF = SW↓ Σ_λ f(λ)[α_clean(λ) − α(λ)] over 300–2500 nm. Here f is BioSNICAR's normalized irradiance spectrum, and SW↓ = S₀ cos θ · 0.75^(1/cos θ) for clear sky; pass `--sw-down` to use measured fluxes instead (e.g. PROMICE). Broadband albedo uses the same weighting.

## Model tiers

| Tier | Absorption | Source |
|---|---|---|
| A | BioSNICAR empirical glacier algae | `ice_algae_empirical_Chevrollier2023` (whole-cell measured optics) |
| B | Level 2 MAC, unpackaged | Phase 1 TD-DFT |
| C | Level 2 MAC × Q\* | Phase 1 TD-DFT + `pigment_packaging.py` |
| D | Fe(III)-complexed + aggregated MAC × Q\* | **Provisional surrogate** until Level 3/4 TD-DFT exists (`--level34-csv`) |

The tier D surrogate (`cell_optics.FePhenolicSurrogate`) works on the Level 2 bands:
- It red-shifts them by 0.10 eV and adds 0.25 eV of extra broadening to represent aggregation.
- It adds an Fe(III)←phenolate LMCT band at 570 nm (FWHM 0.70 eV, ε = 4000 M⁻¹ cm⁻¹ per Fe, one Fe per two ligands). These are typical literature ranges for Fe(III)-catecholate/galloyl complexes, and they **must be cited and varied** before quoting tier D numbers.
- MAC is expressed per kg of ligand, so tiers B–D contain the same number of pigment molecules per cell.

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

- **Cell size (empirical).** The reference cell is *A. nordenskioeldii*, 10.75 × 25.44 µm: Greenland mean volume (Chevrollier et al. 2022) with the measured length:width ratio (Procházková et al. 2021). Fig. 2A shows both species; *A. alaskanum* is 8.95 × 13.08 µm. `--sizes` overrides this.
- **Pigment concentration (empirical).** The phenolic concentration is c_i = 22.0 kg m⁻³: measured phenolics per cell / measured biovolume per cell at S6 (Williamson et al. 2020). The packaging grid uses c_i and c_i × (mean ± 1 SD)/mean.
  - **Stated assumption:** the concentration is the same in both species.
- **Photosynthetic pigments.** Chlorophyll a, chlorophyll b and carotenoids are included in tiers B–D at their measured per-cell concentrations, with in vivo MACs (Williamson et al. 2020). They absorb inside the same packaged cell. `--no-photosynthetic` gives the phenolic-only cell.
- **Not empirical:**
  - g = 0.96.
  - Clear-sky transmissivity 0.75.
  - The Tier D surrogate (provisional until Level 3/4 TD-DFT output exists).
- **Concentration range.** 10⁷ cells mL⁻¹ is far above observed blooms, which reach about 10⁴–10⁵ cells mL⁻¹.
- **Tier A cell size.** Tier A's extinction (7.1×10⁻¹⁰ m² per cell) reflects BioSNICAR's empirical cell size, which differs from our reference cell. Fig. S1 shows the per-cell absorption so the difference stays visible.
- **Forcing scope.** Forcing is instantaneous, clear-sky and direct-beam, with no melt feedbacks.
