# Figure audit and paper-figure suite (P8-FIG-1)

## Paper-figure suite

- **Entry point:** `python3 paper_figures/make_figures.py [--only fig2 fig5 ...]`, run under the shared budget. It only reads saved artifacts: no chemistry or BioSNICAR runs.
- **Outputs** go to `paper_figures/output/`:
  - per figure: a vector PDF, a 300-dpi PNG, `<name>_data.csv` (the plotted values);
  - `manifest.json`: per figure the caption, validation status, source files with sha256, config, and a build record (code commit, script sha256, whether the script had uncommitted changes, matplotlib version).
- **Style:**
  - Palette checked with the colour-vision validator: no FAIL; a contrast WARN on light surfaces, so every series is also identified by a legend or label, and the values are in the CSV.
  - Fixed series colours: Tier A blue, TD-DFT B3LYP orange, CAM-B3LYP yellow, measured MAC green, M2 pink, no-algae grey.
  - Tier C is dashed or hatched.
  - Log axes are labelled "log scale" and explained in captions; any axis-range limit is disclosed in the caption.
- **Visual inspection:** every PNG was rendered and inspected. Fixed after inspection:
  - Fig1: text overflow and arrows crossing boxes;
  - Fig2: shading label colliding with the legend; panel letters colliding with titles;
  - Fig3: legend over the data;
  - Fig4: y-range (−2…1.2) hiding the residuals, clipped B2 interval, clipped title;
  - Fig5: a crash on the H5 split, clipped titles, inconsistent comparator colours;
  - FigS2: two bars sharing a colour (now hatched);
  - FigS3: overlapping direct labels (now a legend, with marker overlap disclosed);
  - FigS5: raw run keys; interpolation lines now dashed and explained.

| File | Content | Evidence role |
|---|---|---|
| Fig1_workflow | Workflow; computed / measured / model / data-role boxes | schematic |
| Fig2_spectra | Raw TD-B3LYP and TD-CAM-B3LYP spectra with sticks; 350 nm hold-region marked; calibrated MAC (median, 95 % of 200 joint draws) vs extract; dE posterior (CAM shift −0.79 eV) | calibration-data comparison, not validation |
| Fig3_substitution | Controlled substitution: plug-in vs posterior contrasts, practical thresholds, D − C visible-absorber ablation, 24 states | model experiment |
| Fig4_forward_diagnosis | Observed vs predicted at the reference nuisance state; spectral residuals; B2/B4/B8 slopes with day-block intervals; zero-count controls; unexplained fraction (no causal label) | development diagnostics |
| Fig5_benchmark | H1, H2/H4 against all baselines; counts, thresholds, simplified M2, H5 (2 test days), deviations D1/D3 | independent tests per protocol (H1 on previously seen sites) |
| Fig6_satellite | **BLOCKED**: no completed, validated scene product exists | — |
| FigS1_chemistry_progress | Durable TD checkpoints (roots converged vs cycles; chunk wall time vs shortest uptime) | reliability |
| FigS2_calibration | Calibration posteriors (dE–w, ln f–φ); visible log-R² by residual model | calibration diagnostics |
| FigS3_coverage_ml | Coverage vs width (H1 vs nuisance-prior ensemble); physics vs ML forward RMSE | development diagnostics |
| FigS4_sensitivity | Reference-state sensitivity of broadband/B2 bias; Monte Carlo stability of the posterior contrasts | sensitivity / numerical convergence |
| FigS5_host_seb | Conditional host SEB: Δα series and modelled melt increments (Modes A/B), not validated melt | conditional scenario |

## Previously committed figures

All predate the 2026-10-09 repairs and are kept for provenance only. Do not cite them in the paper.

| File | Class | Reason |
|---|---|---|
| `phase1/results/level2/level2_B3LYP_spectrum.png` | diagnostic | Quick-look of one TD run; replaced for the paper by Fig2a. |
| `records/phase1_production/FigS0_tddft_calibration.*` | stale | Calibration before the cut-posterior / residual-model repairs; current posteriors are in Fig2 and FigS2. |
| `records/phase3_full/figures/Fig3A_forcing_distributions.*`, `Fig3B_sobol_indices.*`, `Fig3C_interactions.*`, `FigS3_sobol_convergence.*` | stale | Phase 3 global sensitivity with pre-repair optics and calibration; not regenerated. Numbers must not be quoted without a rerun. |
| `records/optics_audit/FigS8_optics_audit.*` | stale | Optics audit before the unit/quantity audit (`docs/quantity_unit_audit.md`). |
| `records/satellite_validation_*/FigS6_satellite_validation.*` (2 files) | superseded | Station comparisons superseded by the frozen-protocol H2/H4 benchmark (Fig5b). |
| `records/multi_scene_2019_tddft/FigS7_seasonal_series.*` | superseded | Development scene series. The satellite interpretation figure (Fig6) is blocked until a validated scene product exists. |

Generators `phase2/figures.py`, `phase3/figures3.py` and `phase4/maps4.py` produced the files above. They are not used by the paper suite.
