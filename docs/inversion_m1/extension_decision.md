# Decision: is the proposed chemistry-variable extension supported by measured variability and independent validation?

**Verdict: Narrow. No new training or chemistry calculation is justified now.** The extension is specifiable only in a restricted form, and its validation is limited to Greenland algal counts.

## Strongest baseline status
* **E (published Dark Zone method, MAP form):** reproduces the reported spectral error (1.79 % mean, 11.9 % max vs 1.8 %, 11.8 %) on its own archive spectra; emulator and target density agree with the original TF/TFP code to float32 precision. Posterior uncertainty was not reproduced.
* **A (our fixed-optics BioSNICAR inversion):** reproduces S6 albedo to RMSE 0.0097 (hemispherical albedo), with a +0.69 dex algal bias against counts that neither a sky mixture nor weak priors changed. Its 4-dex intervals were an error-model artefact (retracted).
* The two baselines use different forward models, observation quantities (albedo vs HCRF) and state vectors; they are not interchangeable and no head-to-head on a common dataset exists.

## The exact supported extension (and nothing larger)
1. **One bounded scalar s per algal surface** multiplying the per-cell absorption of the fixed ice-algal spectrum, with log10 s bounded by the measured between-sample range of phenolic mass per cell (Williamson 2020, 64 samples: median 0.041 ng, p10-p90 0.023-0.100 ng, i.e. about +-0.3 dex around the median; full range 0.001-0.174). Shape is unchanged.
2. Consequence that must be stated up front: **s is exactly degenerate with log10 cells** in the forward model. Only the product (an apparent absorbing load) is identifiable from a spectrum. s can be separated only where cells are counted, and the pigment per cell is measured on the same sample; no such pair exists (Williamson is 2016 habitat samples, the spectra are 2017 plots).
3. A shape variable (free vs Fe-complexed purpurogallin, the digitised Prochazkova 2025 absorbances) is **not supported**: it is a solution-phase digitised figure, not in vivo, not paired to any spectrum or sample, and two spectra cannot separate chemistry from preparation.
4. Dust and cryoconite variability extensions are **not supported**: dust absorption was reported similar across 10 samples in the competitor's own data (and per-sample spectra are not in the archive); cryoconite optics rest on two measurements at one site.

## What the data can and cannot test
| Claim | Can it be tested? | With what | Limit |
|---|---|---|---|
| Spectral prediction (held-out spectral region, held-out day) | yes | S6 2017 (47 plots, 9 days) | development site only; within-region |
| Concentration accuracy | algal term only | counts at S6 2017 (47) and S Greenland 2021 (18) | counts are sample-based, not footprint means; the Greenland 2021 spectra are reused (also the competitor's S21 optics site); no mineral, cryoconite or BC truth |
| Uncertainty calibration | partly | held-out spectral coverage, counts | no scale tried is calibrated for both (0.96 held-out at sigma_total; 50 % of counts) |
| Absorbed energy and melt | **no** | PROMBIO station-day cells with station albedo/SEB exist, with no spectra; spectra and chemistry have no thermal record | nearby is not colocated; modelled melt is not observed melt; 2024 station-days are withheld |
| Himalayan transfer | **no** | no data | - |

A better spectral fit, a narrower interval, or a better count agreement would each be reported separately and none establishes better attribution or melt prediction.

## Closest prior work (what already exists)
* Dark Zone, Chevrollier et al. 2026 (Sci. Adv. 12, eady9482; code and data Zenodo 18826013): Bayesian emulator inversion of algae, cryoconite and dust with a weathering-crust column, measured and fixed optics, 104 field spectra, no concentration truth, no leave-site-out test, no melt evaluation.
* snowlaps (Chevrollier et al. 2025, snow algae, emulator plus gradient descent, no uncertainty): not applicable to glacier-ice algae.
* Bohn et al. 2022 (optimal estimation with PRISMA), Cook et al. 2020 (BioSNICAR-GO retrieval at S6), Chevrollier et al. 2023 (ice-algal optics): earlier work that the above builds on (see `prior_art_comparison.md`; those were read at the level recorded there).
The combination "Bayesian separation of particle types on bare ice with measured optics" is therefore **not novel**. What remains untested in the literature we read is whether bounded, measurement-based variability of per-cell absorption improves **held-out-site** retrieval; that is a hypothesis, and the outcome may be no gain.

## Precise missing evidence
1. Co-located, same-sample: glacier-algal count, phenolic/pigment amount per cell (and ideally in vivo absorption) and a spectrum, on several surfaces and at least two glaciers that are independent of the optics origin.
2. Per-sample mineral absorption spectra with mineralogy tied to spectra (Dark Zone authors hold n = 10; not in the archive).
3. Two or more independently measured cryoconite optical spectra from different sites.
4. Spectra with colocated radiometer or flux and melt records, to connect optics to absorbed energy.
5. Any non-Greenland data for a Himalayan statement.
The cheapest first step is a data request to the Williamson and Chevrollier groups for item 1 and 2 files, not a calculation.

## What was not done
No neural network or PINN was trained or specified (not justified: 47 labelled plots, one site, and the degeneracy above). No Deep Research handoff was implemented. The informative-prior variant is a sensitivity analysis only. The 2.5 nm grid offset in the competitor's pipeline is a hypothesis.
