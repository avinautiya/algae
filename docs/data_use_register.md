# Data-use register (authoritative; supersedes the "Uses so far" column of `docs/data_manifest.md`)

Audited 2026-10-09 by tracing code paths (`grep` of every loader and its callers), not from memory. Keep this file current: any new use of a dataset is added here **before** the run.

## Column codes

| Code | Use |
|---|---|
| CAL | pigment/cell calibration |
| BIO | biological priors (abundance, community) |
| GEO | anisotropy / geometry prior (k = HCRF/albedo) |
| ICE | ice / dust priors and structural ice choices (mode, density, radius grid, dust nodes) |
| DISC | discrepancy / observation-error fitting (σ, τ, albedo discrepancy SD) |
| SEL | hyper-parameter or model / variant selection |
| VAL | validation target |
| INT | final interpretation (maps, effects) |

✔ = used · (✔) = indirect or campaign-level use · — = not used.

| Dataset | CAL | BIO | GEO | ICE | DISC | SEL | VAL | INT |
|---|---|---|---|---|---|---|---|---|
| Williamson et al. 2020, S6 **2016** (counts, pigments/cell, MACs, HPLC) | ✔ TD-DFT calibration targets; measured-MAC optics | ✔ abundance prior N(3.56, 0.78) | — | — | — | (✔) calibration residual model AR(1) vs iid | literature prior baseline in `heldout_v2` | via optics |
| Cook 2020 archive, S6 **2017** counts (47 plots) | — | — | — | — | — | — | ✔ `heldout_v2` targets; H1 inputs (known state); forward diagnostics inputs | — |
| Cook 2020 archive, S6 **2017** HCRF (same plots) | — | — | (✔) numerators of ARF (29 of 51 ARF columns are these plots) | ✔ bubbly-ice choice ("field NIR requires bubbly ice") and radius-grid range set on these spectra; EB radius prior in `field_validation` | ✔ σ (max marginal likelihood), τ; primary-fold calibration of `heldout_v2` | ✔ | ✔ `heldout_v2` reverse-fold target | scene product uses the primary-fold σ, τ and radius prior trained on them |
| biosnicar `Albedo_master.csv` (S6 2017 hemispherical albedo) | — | — | ✔ denominators of ARF → k prior N(0.898, 0.175) | — | — | — | ✔ **H1 target** | — |
| biosnicar `ARF_master.csv` = HCRF/albedo, S6 2017 | — | — | ✔ k prior (all physics retrievals: `heldout_v2`, scene product, H2/H4 retrievals) | — | — | — | — | via retrievals |
| Chevrollier et al. 2023, S Greenland 2021 (counts + HCRF) | — | — | — | — | ✔ σ, τ for the reverse fold | ✔ | ✔ `heldout_v2` primary target; ML comparison | — |
| McCutcheon/Cook S6 2017 dust loading (campaign level) | — | — | — | ✔ dust prior (log-normal) | — | — | — | via priors |
| Ice SSA (Cooper 2021; Dadic 2013) | — | — | — | ✔ radius prior (measured-SSA option) | — | — | — | — |
| Cooper 2018 densities (crust 450, near-surface 690) | — | — | — | ✔ fixed structure | — | — | — | — |
| Community surveys (Williamson 2018; Halbach 2022, 2025) | — | ✔ f_n prior | — | — | — | — | — | — |
| Procházková 2025 Fe absorbance (digitised) | ✔ tier D Fe increment | — | — | — | — | — | — | — |
| PROMICE KAN_M hourly radiation 2016–2019 | — | — | — | — | — | ✔ clear-sky transmissivity T = 0.919 (scales modelled SW and RF only; BBA and Δα = RF/SW do not depend on it) | ✔ SEB validation; H2/H4 station (targets SW↑/albedo only) | multi-scene legacy |
| PROMICE KAN_M heights 2016–2019 | — | — | — | — | ✔ χ = 0.3 fitted on 2016–2018 (scenario, not validation) | — | ✔ SEB/H3 prerequisite checks | — |
| PROMICE KAN_L/KAN_M L3 hourly 2016–2023 (raw) | — | — | — | — | — | — | ✔ H2/H4 targets; H3 prerequisites | — |
| Sentinel-2 L2A (scenes over S6, KAN_L, KAN_M) | — | — | — | — | — | — | inputs for H2/H4 (never targets) | ✔ scene product |
| ESA S2 SRFs | band integration | — | — | — | — | — | — | — |
| PROMBIO 2024 / 2021–2023 (GEUS) | — | — | — | — | — | — | **reserved; roles to be frozen** (`docs/new_data_inventory.md`) | — |
| UPE_U 2018 counts and UAS (BAS PDC) | — | — | — | — | — | — | candidate (not acquired) | — |
| `cook2020_field_metadata.csv` (BioSNICAR-GO retrievals) | — | — | — | — | — | — | reference only, never a truth target | — |
| ML comparison (`records/ml_comparison`) | — | — | — | — | — | — | descriptive baselines on `heldout_v2` rows (retrospective) | — |

## k-prior overlap with S6 test plots (deviation D3), stated precisely

- **Construction.** k ~ N(0.898, 0.175) is the SRF-averaged mean/SD of 51 ARF spectra. ARF = HCRF/albedo of S6 2017 plots, verified numerically for shared plots. 29 of the 51 are counted plots scored in the `heldout_v2` reverse fold (test = S6 2017).
- **Not model-symmetric.** In each retrieval k multiplies that model's forward reflectance. The physics models over-predict broadband HCRF by different amounts (+0.10 to +0.26 at S6), so each model needs a different k to reconcile its own brightness. A k prior centred on the test plots' observed anisotropy helps a model in proportion to how well "model albedo × observed ARF" matches the observed HCRF. That differs between models.
- **Consequence.** The leakage does not cancel in a contrast such as tddft_D − tierA. Its sign and size per model are unknown until a leakage-free rerun exists.
- **The original `heldout_v2` result is retained unchanged**, with this qualification.
- **Planned sensitivity.** A separately named `heldout_v2_kfree` sensitivity, not a replacement:
  1. reverse fold with a k prior built from the training site only. The S Greenland data have no albedo, so no training-only ARF exists. Fall back to an independently justified prior: N(1, 0.25) as used for Stibal 2017 (instrument assumption).
  2. Prior-sensitivity range: k SD × {0.5, 1, 2} and mean ∈ {0.8, 0.9, 1.0}.
  3. Report each model's reverse-fold score across the range.
  - Status: not run (needs the physics tables, about 2 GB; can run under the budget).

## H1: every path by which earlier use of the same plots could reach the implemented prediction

H1 predicts 300–2500 nm broadband albedo of S6 2017 plots from the measured count, by day-held-out empirical Bayes over (radius population, discrepancy SD), with frozen emulators. Paths audited:

| Path | Reaches H1? | Why |
|---|---|---|
| k prior (built from the same albedo) | **No** | `h1_albedo.py` never calls `prior_logpdfs` for k; albedo is predicted directly. |
| σ, τ (fitted on these plots' HCRF) | **No** | not used; H1 fits its own discrepancy SD on training days only. |
| Radius prior | **No** for parameters; **yes** for the grid | EB on H1 training days replaces the HCRF-fitted prior. The EB *grid* (`MU_LNR`, `SD_LNR`) and the emulator radius range 300–20 000 µm were set during earlier work on these plots' HCRF. Structural, campaign-informed. |
| Bubbly-ice mode, crust density, 2 cm crust | **Yes (structural)** | chosen because the field NIR of these S6 spectra "requires bubbly ice". A model-selection decision informed by the same plots' HCRF. It applies identically to M0, M1, M1b and M3 (same ice), so it does not favour one optics model through ice, but H1 is not a test of the ice representation. |
| Dust prior | **Yes (campaign-level)** | from S6 2017 heavy-biomass samples (McCutcheon/Cook), same campaign, not the plots' albedo. Identical across models. |
| Emulators (`heldout_v2/cache`) | no albedo use | built from Phase 1/2 optics; posterior-mean calibration (D1) from S6 2016 data. |
| Counts | input | known state, as designed. |

**Conclusion.** H1 is "new target quantity (albedo) on previously seen sites, with campaign-informed ice structure and dust prior". It is not external evidence. A favourable H1 contrast between optics models is interpretable because ice and dust treatment is identical across models. An absolute H1 error is not interpretable as transfer accuracy.

## Rules

- `heldout_v2` remains immutable. Any leakage-free rerun gets a new name and manifest.
- A dataset listed under SEL, DISC, GEO, ICE or CAL is never relabelled "untouched".
- New datasets (PROMBIO, UPE_U) get development/test roles frozen in `docs/new_data_inventory.md` before any fit.
