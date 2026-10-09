# Molecular-contribution study: results (protocol `docs/molecular_contribution_protocol.md`)

**Inputs and code:**
- `manifest.json`: TD-DFT state/spectrum file sha256, empirical-data fingerprint, code hashes;
- calibrations cached by content fingerprint in `phase4/results/cache/calibrations`.

**Command:** `python3 phase4/molecular_contribution.py --stages raw calibrate compare expA` (budgeted), then `--stages summarise` (the expA summary was re-run from the saved states after a pandas dtype bug; no state was recomputed).

## 1. Raw quantum chemistry vs measured spectra (Q1; measured data were later calibration targets, so this is not independent)

**File:** `raw_vs_measured.csv`, `measured_features.json`.

| | Band maxima | Visible (400–700) share | Unit-area shape rel. RMSE vs isolated chromophore |
|---|---|---|---|
| Isolated chromophore, HPLC peaks 2–4 | 303, 393 nm | 0.137 | — |
| TD-B3LYP raw (FWHM 0.3 eV) | 296, 331, 382 nm | 0.124 | 0.41 |
| TD-CAM-B3LYP raw | 278, 333 nm | 0.0006 | 1.32 |
| Whole extract | — | 0.30 | — |

- Raw B3LYP reproduces the two-band structure of the uncomplexed chromophore, about 0.05–0.1 eV too blue.
- Raw CAM-B3LYP is about 0.5–0.8 eV too blue.
- **Neither represents the extract's additional visible absorption.** Its share is 2.2× that of the isolated chromophore; it comes from complexes or other phenolics.

## 2. What survives calibration (§5)

**File:** `calibration_comparison.csv`, `calibrated_mac_differences.json`, `calibrated_macs.csv`. Identical procedure and freedom: AR(1), same data and bounds.

| | dE (eV) | w (eV) | φ | ln f |
|---|---|---|---|---|
| B3LYP | +0.043 ± 0.008 | 0.61 ± 0.01 | 0.10 ± 0.02 | 4.12 ± 0.76 |
| CAM-B3LYP | **−0.786 ± 0.019** | 0.76 ± 0.02 | 0.06 ± 0.02 | 3.70 ± 0.70 |

- **The calibration absorbs a 0.8 eV functional error** through dE.
- After calibration, the tier-D MACs still differ by 74 % on average in the visible, but only −34 % solar-weighted.
- The B3LYP tier-D posterior itself spans ±124 % (relative SD) in the visible. The visible magnitude is poorly constrained: it is carried by φ × the measured Fe increment and by f.
- Both calibrated tier-D MACs fit the extract only to a log-RMSE of about 0.72–0.77 (a factor of about 2).
- **Root coverage:** CAM's 15 roots move the calibrated MAC by up to 4.2 % when the top 5 are removed (flagged by the code).

## 3. Experiment A, controlled pigment-level substitution (Q2)

**Files:** `expA_states.csv` (2424 states), `expA_pair_differences.csv`, `expA_posterior_spread.csv`.

- **Held identical:** cells, packaging, photosynthetic background, species mix, column, bubbly ice, dust, illumination (SZA 50°, clear sky), observation quantity.
- **States:** log B 3–5, r 1000–8000 µm, dust 3×10⁴ / 3×10⁵ ppb.
- **Uncertainty:** 24 joint calibration draws propagated per TD treatment.

| Pair | max \|Δ broadband albedo\| | max \|Δ B2\| | max \|Δ absorbed SW\| |
|---|---|---|---|
| B3LYP-D − CAM-D (functional, after calibration) | **0.0034** | 0.0074 | **2.5 W m⁻²** |
| B3LYP-D − measured MAC | 0.0051 | 0.0046 | 3.7 |
| CAM-D − measured MAC | 0.0080 | 0.0028 | 5.9 |
| B3LYP-C − CAM-C (no Fe term) | 0.063 | 0.017 | 46.5 |
| B3LYP-D − B3LYP-C (Fe/visible term) | 0.080 | 0.010 | 58.5 |

- **Calibration posterior spread** (max SD over states): B3LYP-D 0.021 broadband (15.8 W m⁻²); CAM-D 0.008 (5.6 W m⁻²).
- **Readings:**
  - After calibration, the choice of functional and the choice "TD-DFT vs measured MAC" change glacier albedo by ≤ 0.008 and absorbed SW by ≤ 6 W m⁻².
  - That is below the 0.01 / 10 W m⁻² meaningful thresholds, and smaller than the calibration's own posterior spread.
  - The only consequential component is the **visible absorber** (tier D vs C: up to 0.08 albedo, 58 W m⁻²), and in the current model it is the **measured** Fe increment scaled by a fitted φ, not a quantum-chemical prediction.

## 4. Forward-error diagnosis, abundance slopes (§7)

**Files:** `albedo_vs_abundance_slopes.csv`, `albedo_slope_day_bootstrap.csv`; `records/forward_diagnostics/`.

Observed vs modelled (TD-DFT D, reference nuisance state) change in hemispherical albedo per decade of counted cells. 41 positive-count S6 2017 plots, 8 sampling days, day-block bootstrap:

| Band | Observed slope (95 %) | Unexplained by the pigment model (95 %) |
|---|---|---|
| B2 490 nm | −0.150 (−0.201, −0.118) | −0.037 (−0.107, +0.011): **not distinguishable from 0** |
| B4 665 nm | −0.134 (−0.183, −0.105) | −0.063 (−0.124, −0.023) |
| B8 833 nm | −0.095 (−0.143, −0.065) | **−0.087 (−0.138, −0.057)** |

- **Every pigment treatment** (Tier A, measured MAC, TD-DFT C/D) predicts essentially no NIR response.
- **About 40 % of the broadband abundance response** (observed −0.106 per decade vs modelled −0.060) comes from darkening that co-varies with algae at wavelengths where algal pigments do not absorb.
- **Consistent with** non-biological surface change co-located with blooms: crust water content, ice micro-structure, cryoconite/dust. It is **not** consistent with an error in the molecular absorption spectrum.
- **Zero-count plots:** −0.14 bias for every model (model too dark): the reference non-algal state is wrong in the other direction for clean surfaces.
- **Geometry** (HCRF/albedo nearly flat spectrally, 0.83–0.88) and **illumination** (SZA ± 5°, spectrum) do not explain the deficit (P7-FWD-1).

## 5. Decision gates

| Gate | Classification | Evidence |
|---|---|---|
| **A** molecular predictions independently supported | **Unvalidated.** No independent laboratory spectrum of the pigment is available. Raw B3LYP agrees with the (later calibration) isolated-chromophore bands within about 0.1 eV; CAM does not. | §1. Required: Remias et al. 2012 spectrum (FEMS Microbiol Ecol 79:638; publisher 403 here) or a new measurement. |
| **B** molecular treatment changes glacier optics meaningfully | **No, after calibration** (≤ 0.008 albedo, ≤ 6 W m⁻² between functionals and vs measured MAC). Yes only for the visible-absorber term, which is empirical. | §3 |
| **C** controlled substitution improves independent observations | **Not supported / underdetermined.** H1 shows no gain. The predicted differences (≤ 0.008) are below the observation floor (0.01), so no glacier dataset can resolve them. | §3, P6-H1-1 |
| **D** transfer without recalibration | **Untested.** No untouched campaign with co-located biology and albedo. | `docs/new_data_inventory.md` |
| **E** melt benefit | **Untested** (H3 untested). | P5-H3-1 |

## 6. Repair decision (§11)

- **No molecular repair is justified.**
  - The molecular treatment is not the binding error: pigment-level substitution moves albedo by ≤ 0.008, while the model–observation gap is 0.10 broadband and 0.24 in B2.
  - Expanding calibration freedom would only manufacture fit.
- **The supported diagnosis is non-molecular:** algae-co-varying red/NIR darkening and a mis-centred non-algal (clean-ice) state.
  - A repair would have to be constrained independently: for example, crust water/structure measured at the plots, or fitted on zero/low-algae controls of a development campaign.
  - It would then need testing on an untouched campaign. No such campaign exists in hand, so **no repaired model is frozen**.

## 7. Where molecular information could still matter (prospective discriminating experiment, §8)

- **Condition.** Molecular treatments differ materially only in the visible absorber (400–650 nm) at high loading (log B ≥ 4.5), where tier C vs D and functionals differ by 0.06–0.08 broadband albedo before an empirical visible term is fitted.
- **Discriminating measurement: laboratory, not satellite.**
  - Absorption of glacier-algal extracts with known pigment mass (absolute, not phenol-equivalent) and known Fe content;
  - the isolated-pigment spectrum (Remias 2012);
  - the computed Fe(III)–purpurogallin complex (`FE_CAT`, now prioritised) vs the measured PG–Fe band (~593 nm, Procházková 2025), after state-character analysis.
- **Glacier observations.** Four-band S2 and broadband albedo cannot discriminate (differences ≤ 0.008 < 0.01 floor).
  - Hyperspectral albedo with co-located HPLC-quantified pigments could, by resolving the 450–650 nm shape at known loading. No such labelled dataset is identified.

## Outcome category

**"Quantum chemistry provides mechanistic explanation (band assignment of the uncomplexed chromophore; B3LYP within about 0.1 eV) while empirical optics predict as well; current data cannot distinguish the molecular treatments in glacier radiation."**

The forward failure that limits prediction is non-molecular.
