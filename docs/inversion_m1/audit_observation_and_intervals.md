# Audit 1: the observation quantity, and what the "4 dex" intervals were

## 1. Which quantity is the S6 spectrum, and in what quantity is the ~0.01 residual?

**Evidence (read, not assumed).** The Cook 2020 archive (Zenodo 10.5281/zenodo.3564501, downloaded; `Albedo_Reflectance_Metadata/13_7_log.docx`) is a bio-albedo logsheet for an ASD FieldSpec with natural illumination, "solar noon +/- 2 hr", fore-optic **RCR/8degree** (remote cosine receiver for the up/down pairs), sensor 30-55 cm above the surface, and for each plot a sequence `ALUP1, ALDN1, ALUP2, ALDN2` (up-facing then down-facing) followed by a separate `HCRF` scan at nadir (bare fibre, Spectralon reference). The `Albedo_master.csv` quantity is therefore **ALDN/ALUP with a cosine receiver: hemispherical spectral albedo under the actual sky**; the `HCRF_master.csv` quantity is nadir **HCRF** (directional). They are different measurements of the same plot and are not interchangeable: their ratio is the anisotropy factor k that the repository previously used as a prior. Baseline A fits the albedo file only and never fits HCRF.

**The 0.0097 median residual (full-window RMSE) is a difference in hemispherical albedo** (dimensionless 0-1) between the model's direct-beam spectral albedo and the 10-nm box-averaged ALDN/ALUP spectrum. The transfer check on S Greenland 2021 and the Dark Zone reproduction use HCRF and are a different quantity.

**Approximations documented.**
1. *Sky condition.* The measurement is under the actual sky (the logsheet records clear sky through thin cirrus to thicker cloud; Dark Zone tables S2/S4 record the same for their spectra). The primary model assumes direct beam at noon SZA. At the fitted parameters the same surface has a **diffuse-minus-direct albedo difference of +0.053 (median RMS 0.053, p90 0.056) over 400-1300 nm** (`records/inversion_m1/baseline_A/sky_sensitivity.json`; BioSNICAR `direct=0` at the fitted parameters, no refit). That is **five times the reported residual**, so the fit's parameters absorb an unmodelled illumination effect and the residual is not a clean measure of fit quality for the physical quantity. Most of the difference is the Fresnel reflection at the air-ice interface (about 3 % for direct beam at 45 degrees against about 6 % hemispherically diffuse), which is not spectrally flat in its effect on retrieved darkening.
2. *Cosine-receiver footprint.* A cosine receptor at 30-55 cm integrates a finite surface patch with a cosine response; the model is plane-parallel and laterally homogeneous. Not quantified.
3. *Band operator.* 1-nm albedo averaged into 10-nm boxes centred on the model bands (`operators.band_average`); tested in `tests/test_operators.py`.
4. *Quality rule.* Bands with fewer than half of the 1-nm samples finite are dropped; albedo must lie in [-0.02, 1.02].
5. *No HCRF-to-albedo conversion is used in Baseline A.* Where HCRF is fitted (S Greenland 2021, Dark Zone), the model albedo is multiplied by a free anisotropy factor k (bounded [0.5, 1.4]), which is an approximation that treats k as wavelength-independent; it does not validate the conversion.
6. *Uncertainty of the observation.* sigma_obs = 0.015 comes from repository documentation and was not re-checked against the source paper; the logsheet gives no per-spectrum uncertainty.

**Sky-mixture sensitivity fit (pre-declared variant, `config_sensitivity_sky.json`).** A fifth parameter, the diffuse fraction f in [0, 1], mixes direct and diffuse albedo linearly. The mixture is an approximation (exact only if both components shared the same spectral weighting). Results are in section 4 below once complete.

## 2. The intervals: construction, and why they were 4 dex

Construction (from `baseline_fixed_optics.py` and `interval_audit.py`):
* profile parameter log10 cells per mL meltwater, grid 1.0, 1.5, ..., 6.5 (the grid ends are the imposed bounds);
* at each grid value the remaining three parameters (log10 dust, log10 dark-particle proxy, ln bubble radius) are re-optimised by bounded least squares from **one warm start** (the global optimum), inside their own bounds. A non-global nuisance optimum can only raise SSE and narrow the interval; the consistency check min(profile SSE) - SSE_opt is >= 0 for all 47 plots (median 1.4e-3, max 1.4e-2), so the optimum is not violated, but no multi-start check of the profile points was run;
* threshold: (SSE(g) - SSE_min) / (sigma x kappa)^2 <= 3.84, sigma and kappa from leave-one-day-out training residuals.

Outcome of the audit (`records/inversion_m1/baseline_A/interval_audit.json`; width on a 0.5-dex grid, linear interpolation of the crossing):

| Error scale | median width (dex) | p10-p90 | lower end bound-limited | closed on both sides |
|---|---|---|---|---|
| as reported, sigma_total x kappa | 4.14 | 3.65-4.59 | 96 % | 4 % |
| sigma_total, no kappa | 0.51 | 0.15-3.51 | 26 % | 74 % |
| sigma_obs = 0.015, no inflation | 0.5 (grid-limited) | 0.02-0.5 | 2 % | 98 % |
| sigma = 0.003 (Dark Zone value) | 0.5 (grid-limited) | 0.31-0.5 | 2 % | 98 % |

* **The 4-dex width was produced by the error model, not the data.** kappa hit its clip of 6 for every day because the AR(1) rule applied to smooth residuals returns rho = 0.95, and several days have sigma_total = 0.056-0.059 (driven by a few large held-out-window residuals). The product sigma x kappa is 0.1-0.35, an order of magnitude above the 0.01 spectral residual. My earlier statement that the algal amount was "not identifiable from albedo alone" was wrong and has been retracted in `baseline_results.md`.
* **Bound-limited vs data-closed.** Under the as-reported scale the interval's lower side reaches the lower bound (log10 cells = 1) in 96 % of plots: bound-limited, i.e. not closed by data. Under the narrower scales almost all intervals close inside the bounds.
* **Narrow does not mean right.** With sigma_obs only, 2 % of counted values lie within their interval (n = 42); with sigma_total (no kappa) 50 %; with the inflated scale 100 % (trivially). The retrieved cells are biased +0.65 dex relative to counts and the five zero-count plots retrieve about 10^4. So intervals conditional on fixed optics, direct-beam sky and a perfect model are far too narrow relative to the systematic discrepancy: **model discrepancy dominates**, and neither the inflated nor the narrow scale is a calibrated uncertainty. No scale tested is calibrated simultaneously for held-out spectra (coverage 0.96 at sigma_total, 0.90 at 0.015, 0.27 at 0.003) and for counts.

## 3. Does ambiguity persist under the competitor's observation model?
See `competitor_reproduction.md` (profile over algae under the archived target). Short answer there once the run is complete.

## 4. Sky-mixture sensitivity and informative priors
See `sky_and_prior_sensitivity.md`.
