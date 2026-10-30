# Sensitivity analyses (declared as sensitivity, not headline)

## Sky mixture (`config_sensitivity_sky.json`, 47 plots, full refit with a diffuse-fraction parameter)
* The fitted diffuse fraction is at its lower end for 68 % of plots (median 7e-10) and at 1 for 4 %; median RMSE 0.0097 is **unchanged** (Baseline A 0.0097); bias against counts +0.70 dex (A: +0.69); dark-particle proxy at bound in 64 % (A: 66 %); dust 21 % (same).
* Reading: the extra physical freedom did not improve the fit or reduce the count bias. The +0.053 forward difference between the two sky end members is therefore not what limits Baseline A; f is weakly constrained (a flat albedo offset is traded against the other terms) and should not be interpreted as a retrieved sky condition. Per-plot cloud values exist in the logsheets and were not used.

## Informative prior on log10 cells (from stored profiles, no refit; `prior_sensitivity.json`)
Prior N(3.75, sd) with 3.75 the midpoint of the imposed bounds (not derived from any count), sd 2, 1, 0.5 dex.
* Interval contraction relative to no prior: median ratio 1.0, 1.03, 1.15 (>= 1: the prior never narrowed the interval, because it pulls against data that prefer higher values); mode shift 0.005-0.09 dex. The data, not the prior, set the width (likelihood width 0.14 dex at sigma 0.015, 0.78 dex at sigma_total, versus prior 95 % widths 2-8 dex).
* Bias against counts 0.53-0.69 dex across priors; RMSE 0.72-0.96 dex. Counts were used only to report outcomes, not to set priors.
* A tighter interval would not have counted as better inference; none was obtained.
