# Prospective reserved temporal test of heterogeneous diffusion

The fixed-midpoint two-layer thermal model lowers equal-profile temperature MAE from0.748°C to0.522°C relative to homogeneous diffusion (30.2%) on the reserved days28–42. It wins6/9profiles, while a training-only boundary-lag ridge remains better in all9. Thus the predeclared consistent-superiority gate fails. This is a useful bounded physical-model improvement, not an established new environmental mechanism or melt-prediction improvement.

| Model | Equal-profile MAE°C | Equal-profile RMSE°C |
|---|---:|---:|
| Two-layer diffusion |0.5217|0.7239|
| Homogeneous diffusion |0.7477|1.0071|
| Quasi-steady interpolation |1.1823|1.7724|
| Training-only lag filter |0.7664|1.0414|
| Training-only ridge |0.2314|0.3821|
| Homogeneous finite-volume control |0.7477|1.0070|

All9eligible profiles,23interior sensor series and7,751reserved hourly predictions are retained. Layered diffusion beats quasi-steady9/9, homogeneous6/9, lag6/9, ridge0/9. Four profiles select equal top/bottom diffusivities, so the data do not consistently demand heterogeneity under this coarse grid and common-capacity assumption. No selected parameter reaches a grid endpoint.

The initial days14–28 benchmark motivated this specific hypothesis; it is therefore development evidence. The present protocol and executable source were committed before scoring the reserved interval, at `4549f6c92b8ea0bca32b58d70f5d28f3007f994c` on `codex/cliff-physics-v2-20261101`. `FREEZE.json` preserves the original scientific source hashes. A separately added prescore arithmetic-audit source hash is recorded in `PRESCORE_AUDIT_SOURCE.json`; it changes no physical method or inputs.

The diffusion equation uses constant common volumetric capacity, two conductivities, a fixed midpoint interface, harmonic face conductances and exact-in-time modal recurrence for linearly interpolated endpoint forcing. Only the first13-day objective selects the81grid candidates. Comparator parameters are preserved from their original early fits. Interior observations initialize the first state once and never reset it. Boundary probes are observed contemporaneously, so these are retrospective component predictions, not autonomous weather forecasts.

The maximum40-to80cell change is0.00515°C, below the0.02°C gate. The homogeneous finite-volume control agrees closely with the original modal solver, excluding a solver discretization change as the principal explanation of the aggregate improvement. Independent checks match every saved truth value to the raw archive and recompute metrics with maximum disagreement1.78e-15. Independent backend tests cover homogeneous analytical decay, steady serial resistance, positivity/conservation and invalid gaps.

Conductivity and volumetric-capacity scales remain unidentifiable from temperature alone. The two fitted ratios do not establish uniquely measured physical strata, and performance could reflect unmodeled advection or sensor/geometry errors. These dependent profiles are from one glacier; no cross-glacier, cliff-lifetime, wetting-cooling, melt-volume or professional-pipeline superiority follows. The result is retained even though the stronger predeclared target failed. No additional complexity is fitted after examining these reserved errors.

Reproduce with `python3 run.py` on a fresh versioned copy (existing final outputs are protected), then `python3 audit.py` and `python3 make_figures.py`. `results/RESULTS.json` contains every training-grid objective, chosen pair, profile score and precision check. The all-profile figure reports all five substantive comparators; the homogeneous finite-volume control is tabulated numerically.
