# Frozen evaluation plan (to be frozen by commit BEFORE any B result or transfer result is viewed)

Frozen identifiers are written to `records/inversion_m1/baseline_A/summary.json` (config hash, `roles_hash`, commit). Changing any rule below after results are seen requires a new, labelled plan version.

## Methods compared
* **A**: fixed-optics bounded least squares (our BioSNICAR forward).
* **E**: published Dark Zone target density (NumPy MAP). Reference for an existing method; different forward model, different state vector.
* **B**: A plus the chemistry-variable parameters of `methods_spec.md` (not runnable until end-member spectra are located).
* **C**: a PINN is **not justified** at this stage: with 47 labelled development plots in one campaign, a trainable network adds flexibility the data cannot constrain, and the physics here is already a differentiable-free solver evaluated directly. C is not run; it may not be called a PINN unless its physical constraints are implemented and documented.
* Equal treatment: A, B (and E where applicable) get the same observations, the same wavelength window, the same error-scale rule, the same bounds-or-prior policy and **no calibration information from the test group**.

## Data roles (from `records/inversion_m1/sample_groups.csv`, unit of independence = sampling day at a site)
| Role | Content | Use |
|---|---|---|
| development | S6 2017, 9 sampling days, 47 count-bearing plots (+39 unlabelled) | error scale, hyperparameters, leave-one-day-out |
| transfer check (reused) | S Greenland 2021 HCRF, 2 days, 1 site; Dark Zone archive S21/S22A HCRF | labelled **reused data**, not a final test; within the same region |
| final test | PROMBIO 2024 station-days (withheld) and `records/heldout_v2` (immutable) | **not touched** by this milestone; no model selection uses it |
| none | Himalayan data | not available |

## Endpoints (reported separately, no composite score)
1. **Spectral reconstruction** within the fit window (400-1300 nm): RMSE, band-averaged residual structure.
2. **Held-out spectral region**: fit 400-1000 nm, predict 1000-1300 nm (NLL and RMSE). Within-spectrum, so not transfer.
3. **Leave-one-day-out** parameter and prediction stability on S6 (development only).
4. **Count agreement** on the 47 colocated plots: rank correlation and log-ratio of retrieved cells (meltwater count) to counted cells, with the caveat that a count is not a retrieved apparent concentration.
5. **Transfer check** on S Greenland 2021 with a free anisotropy factor k (HCRF is not albedo): reconstruction only; labelled reused.
6. **Interval calibration**: coverage of 95 % intervals on held-out days.
7. **Albedo, absorbed energy and melt**: not evaluated; no independent colocated data of the right quantity exist in the manifest. Better attribution does not imply better melt prediction.

## Synthetic tests (separate, never mixed with real endpoints)
Forward-simulated spectra with known parameters, noise at the documented sigma, to test recovery and the degeneracy profile. Passing says the inversion is self-consistent, not that it is right on natural ice.

## Decision rule for B (pre-declared)
Adopt B as an improvement only if, on endpoints 2, 3 and 5, B beats A by a paired day-block bootstrap whose 95 % interval excludes zero, w is resolved (per spec) on a majority of spectra, and interval coverage is not worse. Otherwise report no improvement. A tie or loss is a result.
Selection among variants of B uses development days only.
