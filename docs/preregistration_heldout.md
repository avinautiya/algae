# Pre-registration: held-out evaluation (frozen before any held-out score was computed)

Frozen on 2026-10-09 at ~12:10 UTC, in response to the review of commit 453e751.

State at freezing:
- No held-out score, summary or per-sample table existed.
- `phase4/results/heldout_v1/` contained only the evidence/marginal table for `tierA_empirical`, which holds no scores.
- The run then in progress (an earlier code version) is NOT the confirmatory run; its outputs are not used for any conclusion.

## Primary model (the only model whose held-out performance is a confirmatory result)

**`tddft_D`:**
- Phase 1 Level 2 TD-B3LYP spectrum (`phase1/results/level2`).
- Calibrated on laboratory spectra only (`phase2/tddft_calibration.py`):
  - modular two-stage posterior;
  - residual model AR(1), with the residual correlation ρ marginalised over its grid (not fixed at the boundary MLE);
  - stage 2 by grid quadrature without boundary clipping;
  - tier D (uncomplexed band model + measured Fe increment).
- Optics:
  - Molecular parameters from joint posterior draws.
  - Two measured species mixed by f_n, chlorophylls and carotenoids included.
  - Bubbly ice: densities 450/690 kg m⁻³; dust axis with the measured S6 prior.
- Retrieval:
  - Grid posterior with log B step 0.05.
  - σ and the ice-radius prior by training-site empirical Bayes.
  - Discrepancy τ by maximum likelihood with inner LOO.
- Emulator settings: the `phase4/emulator.py` defaults at this commit, plus `dust_ppb = DUST_NODES_PPB`.

## Reported alternatives (descriptive only; no winner is selected among them)

- `tddft_C` (no Fe term)
- `tddft_D_iid` (independent-residual calibration)
- `measured_mac_C` (measured in vivo MACs)
- `tierA_empirical` (BioSNICAR empirical cell optics; the reference "empirical algal optics")
- Statistical baselines: climatology, literature prior, band ratio, ridge

If an alternative scores better than the primary, that is reported as such. It is not then called the validated model: it would have been selected using the test data.

## Primary comparison and decision rule

- **Metric:** mean log predictive score of the observed counts. The observation model is a Poisson likelihood for S6 (counted volume known; zeros kept) and a log-normal with SD 0.10 dex for S Greenland. Sensitivity runs use 0.05 and 0.20 dex.
- **Primary fold:** train S6 2017 → test S Greenland 2021. The secondary fold is the reverse.
- **Contrast:** `tddft_D` − `tierA_empirical`, paired by sample. The uncertainty is a cluster (block) bootstrap over sampling days, because samples taken on the same day are not independent.
- **"Improvement demonstrated"** requires both of the following:
  - the block-bootstrap 95 % interval on the primary fold excludes 0 in favour of `tddft_D`;
  - the secondary fold has the same sign.
- Otherwise the result is reported as **"improvement not demonstrated"**, whatever the point estimates. With only 2 sampling days at S Greenland, the block interval may be too wide to demonstrate anything; that outcome is reported as such.
- **Secondary metrics** (reported, not used for the decision):
  - RMSE, bias, Spearman ρ, 95 % coverage and CRPS on positive counts;
  - the forward 4-band HCRF log density given the measured abundance;
  - broadband HCRF error.

## Scope

This evaluation tests abundance retrieval and reflectance prediction at an unseen site. It does **not** test melt prediction. No melt observation co-located with these samples exists, so no statement about better glacier-melt prediction follows from it.
