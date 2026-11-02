# Observed-boundary field-temperature benchmark

All 11 eligible Kennicott debris profiles completed the frozen test: 27 interior sensor series and 9,099 hourly withheld predictions. Homogeneous diffusion improves on depth-linear quasi-steady interpolation in every profile, but does not outperform the strongest statistical comparator. This is component validation, not field validation of cliff melt, hydraulic recovery or a forecast driven by predicted weather.

| Model | Equal-profile MAE (°C) | Equal-profile RMSE (°C) |
|---|---:|---:|
| Homogeneous diffusion | 0.7552 | 1.0198 |
| Quasi-steady interpolation | 1.1933 | 1.7124 |
| Training-only lag filter | 0.7534 | 1.0161 |
| Training-only boundary-lag ridge | 0.2863 | 0.5286 |

Diffusion reduces mean MAE by 36.7% relative to quasi-steady interpolation. It wins against that baseline in 11/11 profiles, against the lag filter in 7/11, and against ridge in 0/11. Consistent superiority is therefore **not established**. All profiles and comparator scores are retained, rather than selecting favorable depths or years.

## Design and interpretation

The original protocol, raw-source hashes, sensor depths, eligible cohort and splits were fixed before scoring. A later append-only audit-metadata update is recorded separately in `results/PREFIT_AUDIT_UPDATE.json`; original freeze files remain unchanged. Each model is fitted only on elapsed days 1–14, initialized once at the first observation, and evaluated on days 14–28. Measured shallowest and deepest probe temperatures are imposed throughout. Interior test temperatures are never assimilated. Training-only ridge uses causal boundary lags and a training-internal penalty selection; no held-out target sets its coefficients.

The 2020b profile is excluded because its apparent deepest temperature column is an author-constructed zero boundary rather than an observation. Negative observed temperatures in other profiles are retained. This linear sensible-heat subcolumn calculation does not assert a freezing-phase model. Every eligible profile, including failures of the physical approximation, is included.

Profiles belong to one glacier, and sensors, hours and profiles share weather and physical context. The 9,099 predictions are not 9,099 independent field replications. No statistical confidence interval or cross-glacier superiority follows. The fitted diffusivity identifies k/C only: jointly rescaling conductivity and volumetric capacity leaves temperature predictions unchanged while changing flux. Consequently this experiment cannot validate melt rates.

## Numerical and reproducibility checks

All fitted diffusivities are interior to the 1e-8–1e-5 m²/s search interval. The maximum 64-versus-128-mode temperature difference is 0.000183°C, far below the predeclared 0.02°C gate. Poorer field accuracy is therefore not explained by this modal truncation or an optimizer trapped at the search boundary.

An independent saved-result audit recomputed all metrics and checked all 17,523 train/test prediction rows against raw timestamps, depths and values. Maximum arithmetic disagreement was 1.78e-15. A training-only ridge replay reproduced predictions exactly, with fitted weights saved separately. Two synthetic leakage tests verify that modifying withheld interior values cannot change parameters/predictions, and changing future boundary inputs cannot change fitted parameters.

The executed solver source is preserved as `../solver/dirichlet_slab_temperature_v1_executed.py`, SHA256 `ef965a927157529152c068b5a5b21bee0226bfec2727e71c0d806dd5b0dcb04a`. Later gradient-API additions to the main solver do not rewrite this execution freeze. The executed runner is retained unchanged.

The physical model's failure against ridge motivates a separate, prospectively frozen two-layer test; it does not demonstrate that heterogeneity is the unique missing mechanism. Advection, sensor representativeness and geometry uncertainty remain possible explanations.

## Figures

`figures/field_temperature_all_profiles.pdf` reports every profile and all four models. `figures/field_temperature_example.pdf` shows a deterministic illustration: the first sorted eligible profile, its deepest interior probe and the first 73 hourly test samples. That example is not selected for agreement.
