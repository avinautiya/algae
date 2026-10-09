# Quantity / unit / normalization audit: TD-DFT → glacier radiation (2026-10-09, from current code)

**Status codes:**
- **M**: directly measured;
- **I**: inferred from measurements;
- **C**: calibrated (fitted);
- **A**: assumed;
- **Q**: quantum-chemical.

| # | Source → destination | Units | Convention / equation (code) | Constants, mass, grid | Status | Assumptions | Tests |
|---|---|---|---|---|---|---|---|
| 1 | TD-DFT states → stick spectrum | E_i (eV), f_i (dimensionless) | PySCF TD/TDA, PCM water ε 78.4 (`phase1/run_phase1.py`) | B3LYP: 30 roots full TD; CAM: 15 roots TDA; 6-31G* (Cartesian d); neutral glucoside C18H18O12 | Q | neutral species; one conformer; B3LYP geometry provisional (optimisation not converged, `level2/COMPLETE.json`) | `phase1/tests` (root coverage, checkpoint validity) |
| 2 | sticks → ε(E) | L mol⁻¹ cm⁻¹ | Hilborn: ∫ε dν̃ = f/4.319×10⁻⁹; ε(E) = K Σ f_i g_i(E), K = 2.8707×10⁴ (`spectra.gaussian_broaden`) | area-normalised Gaussian in ENERGY, FWHM 0.3 eV raw (calibrated w later) | Q (+A width) | ε decadic; symmetric energy broadening; ε is a coefficient, so ε(λ) = ε(E(λ)) with no Jacobian (correct for a cross-section) | `phase1/tests` normalisation integral |
| 3 | ε → MAC | m² kg⁻¹ (Napierian) | MAC = ln10 · 0.1 · ε / M × 1000 (`spectra.epsilon_to_mac`) | M = 426.33 g mol⁻¹ (glucoside, not aglycone) | Q | Napierian, as BioSNICAR consumes | unit test |
| 4 | raw MAC → calibrated shape | dE (eV), w (eV FWHM) | stage 1: unit-area shape vs HPLC isolated chromophore 265–600 nm, AR(1) residuals (`tddft_calibration._shape_loglik`) | 1-nm grid; trapezoid area | C | the isolated HPLC peaks 2–4 are the uncomplexed chromophore | `phase2/tests` calibration |
| 5 | calibrated shape → per-mass MAC | f (dimensionless) | stage 2: ln E_extract = ln f + ln[M(dE, w) + φ I_M D] (log-space, 260–750 nm) | E_extract per kg **phenol equivalents** (EPA 420.1); f absorbs the glucoside → phenol-equivalent conversion AND the oscillator-strength error | C | not separable: an assay normalisation and an f-error are confounded in f | — |
| 6 | Fe complex term | φ ∈ [0, 1], D (per unit integrated PG absorbance) | φ · I_M · D(λ), D = measured PG-Fe minus PG absorbance (Procházková 2025, digitised) | I_M = ∫ M over 265–600 nm | M (D), C (φ) | φ is an effective mixture coefficient, **not** a measured complexed fraction; D measured for the aglycone purpurogallin, applied to the glucoside | — |
| 7 | MAC function → 480-band grid | m² kg⁻¹ | `cell_optics.to_480`: evaluated 350–800 nm; **held at the 350 nm value below 350 nm**; 0 above 800 nm | BioSNICAR 10-nm band centres | A | the UV TD-DFT bands (B3LYP 296/331 nm) never reach the glacier model as computed; the measured-MAC treatment is 0 outside 250–750 nm | — |
| 8 | MAC → intracellular absorption | m⁻¹ | a = MAC × c_int / vacuole fraction | c_int = phenolic mass (phenol eq.) per cell / pooled S6 2016 biovolume (`ED.intracellular_concentration_kg_m3`) | I | equal intracellular concentration across species, scaled by measured size dependence γ | `phase2/tests` |
| 9 | photosynthetic background | m⁻¹ | chl a, chl b, carotenoid in-vivo MACs × per-cell mass / biovolume | Williamson 2020 / Dauchet 2015 | M/I | shares self-shading with phenolics | — |
| 10 | cell packaging and scattering | m² per cell | Q* by ray-chord (cylinders); scattering = 2A − absorption; g by Mie of the equal-volume sphere (n 1.38 / 1.31) (`CellModel.optics`) | geometry: Halbach 2022 volumes, Procházková 2021 L:W | A/M | geometric-optics extinction 2A; homogeneous cell | packaging tests |
| 11 | cells mL⁻¹ meltwater → column loading | cells m⁻² | BioSNICAR conc × 0.917 so column = conc[cells/g] × ρ × dz (`biosnicar_bridge._layer_concs`) | crust 2 cm at 450 kg m⁻³ over 690 | M (counts), A (crust) | the count sample represents the top 2 cm of crust | bridge tests |
| 12 | column RT → spectral albedo | 1 | BioSNICAR adding-doubling, bubbly ice, sub-Arctic summer clear-sky direct beam | 480 bands, 200–5000 nm | A (ice state) | 2-layer column; dust uniform in the crust | — |
| 13 | albedo → observable | 1 | broadband = irradiance-weighted 300–2500 nm; S2 bands = SRF × flux weighting; HCRF = k × albedo (scalar k) | k ~ N(0.898, 0.175) from S6 ARF (D3) | A | **scalar k cannot represent wavelength-dependent anisotropy** | — |
| 14 | albedo → absorbed SW / forcing | W m⁻² | SW × Σ flux (α_clean − α) | clear-sky SW parameterisation (T = 0.919 fitted to KAN_M) or measured SW↓ | M/C | clear sky for emulator RF | SEB tests |

## Answers to specific checks

- **Oscillator-strength normalisation and integration:** Hilborn relation, energy-domain area-normalised Gaussians. **Correct.**
- **Wavelength vs energy broadening, Jacobian:** broadening in energy; ε is evaluated pointwise, so no Jacobian is needed. A Jacobian would be needed only for a density per nm, which is not used. **Correct.**
- **Napierian vs decadic:** Napierian throughout; the measured extract MAC is taken as Napierian per the source loader. Verify against the Williamson deposit: listed as m² mg⁻¹.
- **Molar vs mass normalisation; phenol-equivalent vs pigment mass; glycoside vs aglycone:**
  - the raw MAC is per glucoside mass;
  - every downstream magnitude is per phenol-equivalent mass via the fitted f;
  - the magnitude of the TD-DFT oscillator strengths never reaches the glacier model.
  - The assay conversion and the f-error are **confounded** in f, which is reported as such, never interpreted as an oscillator-strength error.
- **Pigment mixture:** the calibration shape uses HPLC peaks 2–4 (purpurogallin-type); the extract contains complexes and other phenolics, carried by φ·D and by f.
- **Extract vs isolated vs in vivo:** shape comes from the isolated chromophore; magnitude from the extract (in solution); packaging is applied in the cell model. In-vivo cell absorption is never fitted, so there is no double packaging.
- **cells mL⁻¹ meltwater vs solid ice; scrape depth; crust density:** handled in `_layer_concs` (row 11). The column loading depends on the assumed 2 cm crust at 450 kg m⁻³; the field scrape depth is about 2 cm (S6) or 2–5 cm (PROMBIO).
- **Interpolation/extrapolation:** row 7 (UV hold, red zero). The calibration fit window is 260–750 nm; the Fe increment is held at its 255 nm value below 255 nm.
- **Root coverage:** CAM covers ≥ 195 nm (15 roots); B3LYP has 30 roots (highest about 5.25 eV ≈ 236 nm). The consumed window starts at 350 nm (glacier) or 260 nm (calibration). Root-drop sensitivity: `tddft_calibration.root_count_sensitivity`.
