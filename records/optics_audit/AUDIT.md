# Optics audit: the 2–7× low-biomass forcing gap vs published estimates

**Question.** Our daily algal forcing and melt at low and medium abundance are 2–7× below Cook et al. (2020) Lbio and Williamson et al. (2020). Is there a bug or an improper assumption in the molecule → cell → BioSNICAR chain?

**Answer.** No bug was found. Our per-cell absorption matches laboratory-derived in vivo glacier-algal cross-sections to within 7 % (band mean, 350–700 nm). At the published abundances, three independent optics give the same darkening to within 5 %. The gap comes from how the published numbers were derived (Section 5).

Reproduce:

    python phase2/audit_optics.py

Outputs: `optics_audit.json`, `optics_audit_curves.csv`, `optics_audit_consequence.csv`, Fig. S8.

## 1. Assay-factor integration (f = 39.5 = 8.7 × the 1:1 glucoside→phenol value)

**Code path.**
- `phase2/tddft_calibration.py`, `complexed_mac()`: `return f * (m(wl) + phi * i_m * ...)`. Here f multiplies the TD-DFT MAC per kg glucoside, giving a MAC **per kg phenol equivalent**. That is the unit of the measured extract MAC it is fitted to.
- `phase4/emulator.py`, `_Builder.__init__`: `mac = co.to_480(cal.mac_D if cfg.tier == "D" else cal.mac_C)`.
- `phase2/cell_optics.py`, `CellModel.optics`: `a_pig = mac480 * c_comp`. Here `c_comp` = `empirical_phenolic_concentration()` = 22.0 kg m⁻³ of cell volume, **also in phenol equivalents** (same 4-AAP assay, Williamson et al. 2020).

**Numerical check.**

| Quantity | Value |
|---|---|
| Tier D MAC / measured extract MAC, 350–700 nm | median 0.82 (range 0.50–1.52) |
| Intracellular absorption a_i at 400 / 550 / 650 nm | 4.9 × 10⁶ / 7.7 × 10⁵ / 6.1 × 10⁵ m⁻¹ |

The f factor is applied once, in the right direction, and in consistent units. It cancels in a_i: a ninefold error in the absolute pigment mass would scale MAC and c_i inversely, leaving their product unchanged.

## 2. Packaging (Duysens / Morel–Bricaud, chord-averaged cylinder)

**Code path.**
- `phase2/pigment_packaging.py`, `q_star()`: Q* = ⟨1 − e^(−a l)⟩ / (a ⟨l⟩) over isotropic chords.
- `cell_optics.CellModel.optics`: σ_abs = MAC × m_pig × Q*.

**Checks.**

| Check | Result |
|---|---|
| Density independence | `CellModel.optics` has no concentration argument. Abundance enters only BioSNICAR's `conc` (`emulator._Builder.node`), so Q* cannot depend on cell density |
| Q* over 350–700 nm | Mean 0.12 (*A. nordenskioeldii*), 0.08 (*A. alaskanum*); range 0.008–1.0 |
| σ_abs (packaged) / projected area | 0.96 and 0.97: the cells are essentially black over 350–700 nm, the physical limit |
| σ_abs unpackaged / projected area | 16–25× (band mean; 4–88× by wavelength) |

The unpackaged (dissolved-pigment) cross-section is physically impossible: a particle cannot absorb more light than falls on it in ray optics. Packaging is therefore not "over-suppressing". It enforces the geometric limit.

## 3. Units into BioSNICAR

**Code path.**
- `phase2/biosnicar_bridge.py`: `CustomImpurity` (`unit = 1`, `mac` = ext_xsc in m² cell⁻¹) and `_layer_concs` (× 0.917).
- BioSNICAR `column_OPs.mix_in_impurities`:
  - cells kg⁻¹ = `conc / 917 * 10**6`;
  - τ = ρ·dz × cells kg⁻¹ × `mac` (comment in the source: "cells m-2 * m2 cells-1").

**Numerical check.**

| Quantity | Value |
|---|---|
| Column cells at 10⁴ cells mL⁻¹ (meltwater), ρ = 450 kg m⁻³, dz = 2 cm | 9.0 × 10⁷ cells m⁻² = count [g⁻¹] × ρ [g m⁻³] × dz, as measured |
| Our ext_xsc at 500 nm | 5.2 × 10⁻¹⁰ m² cell⁻¹ |
| BioSNICAR's own glacier-algae entry at 500 nm | 6.4 × 10⁻¹⁰ m² cell⁻¹ (same unit and magnitude) |
| Algal optical depth at 500 nm | 0.047 |
| Area covered by cells | 2.3 % |

Units are consistent.

## 4. Absolute per-cell absorption vs measured cells (linear scale, community 60 % *A. nordenskioeldii*)

Reference: `ice_algae_empirical_Chevrollier2023` (BioSNICAR's default glacier-algae optics, derived from measured in vivo absorption; our tier A).

| Band (nm) | Ours, tier D (10⁻¹⁰ m²) | Measured cells (10⁻¹⁰ m²) | Ratio |
|---|---|---|---|
| 350–400 | 2.04 | 2.75 | 0.74 |
| 400–450 | 2.03 | 2.58 | 0.79 |
| 450–500 | 1.99 | 2.36 | 0.84 |
| 500–550 | 1.96 | 2.17 | 0.90 |
| 550–600 | 1.95 | 2.07 | 0.94 |
| 600–650 | 1.94 | 1.70 | 1.14 |
| 650–700 | 1.91 | 1.20 | 1.59 |
| **350–700 mean** | **1.98** | **2.12** | **0.93** |

The measured extract MAC put through the same cell model also gives 0.93.

Two physical (not coding) differences:
- **Below 600 nm** the measured cells absorb up to 40 % more than the projected area of our Halbach-volume cells. They must be larger, and our cells are already at the black limit.
- **Above 620 nm** the measured cells become transparent, while the Fe-complex tail of tier D keeps ours dark.

## 5. Consequence at the published abundances

Same ice column (bubbly, r = 600 µm, ρ = 450/690 kg m⁻³) and S6 dust; SZA 47°; clear sky.

| Abundance (cells mL⁻¹) | Optics | Δ broadband albedo | Noon RF (W m⁻²) |
|---|---|---|---|
| 2.9 × 10⁴ (Cook Hbio) | ours, tier D | 0.093 | 73 |
| | measured extract MAC | 0.088 | 69 |
| | Chevrollier 2023 measured cells | 0.091 | 71 |
| 4.73 × 10³ (Cook Lbio) | ours, tier D | 0.019 | 14.8 |
| | measured extract MAC | 0.018 | 14.1 |
| | Chevrollier 2023 | 0.019 | 14.9 |
| 8989 (Williamson high) | ours, tier D | 0.034 | 27.0 |
| | measured extract MAC | 0.033 | 25.7 |
| | Chevrollier 2023 | 0.034 | 27.0 |

**Crust density (the one assumption with leverage).** Cells per m² = count × ρ × dz.

| ρ_top (kg m⁻³) | Lbio noon RF | Hbio noon RF |
|---|---|---|
| 450 | 14.8 | 73 |
| 690 | 21.0 | 99 |
| 850 | 24.2 | 114 |

The Phase 3 prior is U(330, 560) kg m⁻³ (measured S6 weathering crust, Cooper et al. 2018). Even ρ = 850 leaves Lbio about 5× below Cook's.

**Where the published numbers come from.**

- **Cook et al. (2020).** RF is the **measured albedo difference between algal sites and clean ice**. UAV broadband albedo (Table 2): Hbio 0.25, Lbio 0.44, clean ice 0.53, so Δ = 0.28 and 0.09 (0.28 × 372 W m⁻² ≈ 104 ≈ 116 W m⁻²). This attributes all the extra darkness of algal sites to algae, including dust, cryoconite and weathering-crust structure.
  - Our field validation reproduces those same spectra with algae + dust + ice structure. The counted abundance is recovered with bias +0.25 dex (S6) and −0.06 dex (independent site).
  - Had our per-cell darkening been 2–7× too low, the retrieval would have needed 2–7× more cells: bias +0.3 to +0.85 dex. That rules out a gap of that size in our optics. At S6, at most a factor 10^0.25 ≈ 1.8 could be missing, which may equally be non-algal darkening.
- **Williamson et al. (2020).** Per-cell light interception is "half the lateral surface area" of a 29.6 × 12.0 µm cell (560 µm²). At 8989 cells mL⁻¹ sampled over 20 × 20 × 2 cm (≈ 1 × 10⁸ cells m⁻²), perfectly black cells would intercept at most about 6 % of the surface, and only below 750 nm (≈ 55 % of SW).
  - Upper bound: about 11 W m⁻² daily mean with 329 W m⁻² mean irradiance. Ours is 9.4.
  - Their 1.86 cm w.e. d⁻¹ needs about 72 W m⁻², i.e. about 6× more cells per m² than the stated sampling implies.
  - The volume-to-area "correction factor of ref. 23" is not stated in either paper. It is the most likely source.

## Conclusion for the paper

- **Code.** The molecule → cell → BioSNICAR chain is unit-consistent and reproduces measured in vivo glacier-algal absorption cross-sections (ratio 0.93).
- **Comparison.** The literature comparison should present ours as **algae-only forcing**: the cells' own absorption, physically capped by their projected area.
- **Cook et al. 2020.** Site-differencing values include non-algal darkening of algal sites.
- **Williamson et al. 2020.** That melt estimate is not reproducible from its stated cell geometry and sampling. It exceeds the black-cell geometric limit by about 6×.

No code change is warranted. The one assumption that moves the result is crust density, already sampled in Phase 3 (S_T 0.16 for forcing efficiency).
