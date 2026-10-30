# Baseline A results (what actually executed)

Command: `python3 common/run_budgeted.py --background --threads 1 --mem-mb 1500 --name m1_fit -- python3 inversion_m1/baseline_fixed_optics.py fit`, then `... evaluate`.
Data: S6 2017 hemispherical albedo, 47 plots with colocated counts, 9 sampling days. Config sha256 `7c9fe755ccff99a633f5418bb6188f1e5ab888593573a81c022da91985850ac4`; roles hash `3bce82ee9756a21317ec7b7c5ff60562eeafa60917ea416e5b6e51c6762a28c9`; BioSNICAR commit `fe74eeef2999fd285f3dc06c47b0f0e1a36a17a6`; library versions in `records/inversion_m1/baseline_A/summary.json`.
Observation operator: 10 nm box average of the 1 nm albedo; window 400-1300 nm; bounds-only (no informative prior).

| Quantity | Value |
|---|---|
| Full-window RMSE (median over 47) | 0.0097 (IQR 0.0075-0.0112) |
| Held-out 1000-1300 nm MAE (fit on 400-1000), median | 0.0116; mean bias -0.0075 |
| 95 % band coverage of the held-out window, leave-one-day-out | 0.957, but with sigma_total 0.016-0.059 inflated by kappa 5.8-6.0 (the clip is 6): coverage is reached by inflation, so it is not evidence of a well-specified error model |
| Algae vs counts (47 plots) | Spearman 0.78; bias +0.69 dex; RMSE 1.0 dex |
| Profile interval for log10 cells | median width 4.2 dex; **open at an end for 100 %**; wider than 1 dex for 100 % |
| Interval coverage of counted cells | 1.00 (uninformative because intervals span most of the bounds) |
| Plots with zero counted cells (n = 5) | retrieved log10 cells median 4.1 |
| Parameters at a bound | dust 21 %, dark-particle proxy 66 %, cells 2 %, radius 4 % |
| Median abs. correlation | cells-dust 0.91, cells-dark 0.95, dust-dark 0.82; median condition number 3.2e6 |
| Broadband albedo (H1 convention) | bias -0.0067, MAE 0.0084 against measured broadband (n = 47) |
| vs. published Cook 2020 retrievals (13 shared) | Spearman 0.55 between retrievals; RMSE vs counts 0.48 dex (ours) vs 0.24 dex (published) |

## Reading it without overclaiming
* The spectral shape is reproduced to ~0.01 albedo, yet the algal amount is **not identifiable from albedo alone** with free mineral and free dark-particle terms: profile intervals are essentially the whole bound range. The rank ordering (Spearman 0.78) carries some information; the absolute amount does not, and the bias is positive (+0.7 dex), with spurious retrievals on count-zero plots.
* This is a property of Baseline A **as specified here** (three free loadings, bounds only, albedo only). It is not a statement about the Dark Zone method, which uses a different forward model, informative Beta priors and 5 nm hyperspectral data to 1200 nm.
* The published retrievals are model output, used as a reference only. Our lower agreement with counts than the published values (0.48 vs 0.24 dex) is on 13 plots only and should not be generalised.
* The pre-declared informative-prior variant has **not** been run; it is the natural next test of whether priors, rather than chemistry, close the gap.
* None of this addresses melt or absorbed energy.
