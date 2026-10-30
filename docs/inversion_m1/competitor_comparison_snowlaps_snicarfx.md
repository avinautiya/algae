# snowlaps vs SNICAR-fx / Dark Zone inversion: code-level comparison

Scope: what the two existing implementations do, so that our proposed contribution is stated against them.
It is not a statement that either is wrong. A limitation of snowlaps is not evidence that our approach is new.

## Sources and pins

| Item | Pin | How read |
|---|---|---|
| snowlaps (snow) | github.com/... snowlaps checkout, commit `01a7ec62ec3b01a08a161b253ad05456d5649a54` | code read in full for the inversion path |
| SNICAR-fx forward model | github.com/openosmia/snicar-fx, commit `c1780abdcc3126c34580377a3f0135fcfccaa25c` | code + README; **contains no inversion code** (grep for bayes/prior/likelihood/invert/retrieve: no hits outside the README citation) |
| Dark Zone inversion code and data | Zenodo `10.5281/zenodo.18826013`, `data_software.zip`, md5 `0955d9c3ddf93ff9b69229d4ac28aa17` (3,036,631 B); files `src/mcmc_pipeline.py`, `src/emulator.py`, `data/emulator/*`, `data/optical_properties/*`, `data/spectra/*` | downloaded to an isolated directory; the sklearn scaler pickle was **not** unpickled (constants read from the raw bytes) and TensorFlow was not installed or run |
| Paper | Chevrollier et al., Sci. Adv. 12, eady9482 (2026), full text via Europe PMC PMC13371923 | read: Results, Discussion, Methods. Supplementary tables S1-S4 and figs S1-S8 **not** read |
| SNICAR-fx release used in the paper | Zenodo `10.5281/zenodo.18609613` | not downloaded; whether commit c1780ab equals that release is **unverified** |

Labels used below: **[A]** assumption made by the method, **[L]** limitation acknowledged by the authors,
**[D]** failure demonstrated in the authors' own results, **[O]** observed by us in the code/files and not stated by the authors,
**[U]** not verifiable from what we could read.

## Snow algae vs glacier ice algae (not interchangeable)

* snowlaps uses **red snow algae** optics: `data/inputs/snicar_config_for_emulator.yaml` line 60 `snow_algae_empirical_Chevrollier2023.nc`
  (Chevrollier 2025 TC, snow, Norway/Svalbard-type surfaces), with dust `dust_skiles_size5.nc` (line 67, Colorado dust) and BC.
* The Dark Zone work uses **ice algae** in vivo absorption measured at S21 (`ice_algae_1nm_resolution.nc`: "scattering is neglected ... no absorption beyond 830 nm is assumed"), a different pigment assemblage (purpurogallin-type phenolics dominate ice algae; snow algae carry astaxanthin-type carotenoids) and a different environment (weathering crust, liquid water, cryoconite).
  SNICAR-fx ships both `snow_algae.nc` and `ice_algae_1nm_resolution.nc`, so the distinction is available in the software.
* Consequence for us: a snowlaps result on snow says nothing about ice-algal retrieval, and our fixed-optics Baseline A uses the ice-algal file (Chevrollier 2023) for the same reason. Pigment statements in this paragraph about compound classes are general background and were **not** verified against these papers' text here.

## Comparison on the seven requested axes

### 1. Fixed vs variable absorber spectra
| | snowlaps | SNICAR-fx / Dark Zone |
|---|---|---|
| Spectra | **Fixed** per absorber, set in the emulator's training config (algae, dust, BC files above) [A] | **Fixed** per absorber: one "representative" mineral mixture, one ice-algae spectrum from S21, one cryoconite spectrum [A] |
| Variability inferred? | No. Only concentrations, grain size, LWC and SZA are free (`snowlaps.py` ~l.152 `gradient_mask`, output columns l.~372-380) | No. Inferred: 2 densities, 2 layer depths, 3 concentrations (`mcmc_pipeline.py` `n_params = 7`, l.52) |
| Spectral flexibility | none | none (no per-wavelength correction) |

The paper states the cryoconite spectrum is **inferred from two reflectance measurements of optically thick layers at S21**, "manually inverted with SNICAR-fx assuming a constant asymmetry parameter of 0.8" (Methods, Model development) [A]; and the file note says scattering of algae and dust is neglected. Whether this is adequate across sites is not tested.

### 2. Mineral-composition and pigment variability
* **Mineral:** the paper measured **n = 10** dust samples across the three sites (filter-pad absorption, XRD), found consistent mineralogy and mass absorption, and used one mixture; dust was 2-5x less absorbing than mixtures elsewhere. So for this region **measured** mineral variability was small and was reported as such. This weakens a mineral-chemistry-variability extension for Greenland specifically; it does not speak to other regions (no Himalayan data exist in our manifest) [L/measured].
* **Pigment:** one in vivo algal spectrum "shown to be representative of different sites in the southwestern ablation area" (their ref. 35 = the 2023 J. Glaciol. paper; **ref. 35's representativeness evidence was not re-read by us**). Seasonal or per-bloom pigment change is not represented. Not tested either way [U].
* **Cryoconite:** the least constrained optical input (derived from 2 spectra, one site). The authors name cryoconite vs dense weathering crust as spectrally confusable [L].

### 3. Cryoconite, liquid water and weathering-crust confounding
* snowlaps has **no** cryoconite, no weathering crust, LWC fixed to 0 in the emulator config (`snicar_config_for_emulator.yaml` line 25 `LWC: [0, 0]`) although `lwc` appears as an output column (l.~372) [O: the column is returned although the config sets LWC to zero; whether the emulator input is live was not tested]. It is a snow tool; ice confounding is out of its scope [A].
* Dark Zone: 3-layer column (granular surface layer, weathering crust, impermeable ice), liquid water as 10 % of density in the two upper layers, particles in the top 5 mm, Fresnel interface (Methods, Model configuration) [A]. BC excluded because Greenland BC is low [A].
* **Demonstrated confounding [D]:** errors > 4 % on fully saturated surfaces "resembled liquid water"; 9 of 174 spectra had multimodal or non-converged chains (3 saturated surfaces with two modes, 2 with one bad chain, 2 single-chain, 2 oscillating, shadowing); 10 of ~104 retrievals (<10 %) deviated noticeably from expectation, half attributed to shadowing, the rest to crust/cryoconite similarity or unseen subsurface cryoconite. "the presence of cryoconite was more difficult to systematically distinguish from a dense ice column". These are the authors' own admissions, and they are the sharpest published statement of the degeneracy that our profile-likelihood diagnostic targets.

### 4. Independent concentration validation vs spectral reconstruction
* snowlaps: validated by reconstruction MAE (best-of-20 starts) and, in the original paper, compared with sample measurements. We recompute the shipped per-spectrum MAE (see `reproduce_snowlaps.py`, pending). Reconstruction quality is what the code returns (`mean_MAE`, l.~362).
* Dark Zone: reconstruction error 1.8 % mean (max 11.8 %, n = 104) and 1.6 % (max 3.7 %, n = 70) on earlier campaigns [measured]; concentrations were checked **only against qualitative surface descriptions and photographs**, not against measured cell counts or cryoconite mass for the same footprint. The spectra file (105 spectra, tags like `050821_S10_HCRF`) carries no concentration variable [O].
* Lack of independent concentration data weakens the claim that the separated amounts are correct. It does not show the inversion is wrong. A good fit also cannot show it is right, because the same fitted model is used to judge itself ("perfect model assumption", Methods).

### 5. Cross-site transfer and exclusions
* Neither implementation fits a model on some sites and tests on others. The emulators are trained on **simulated** data (10.5 M SNICAR-fx runs; train/test split 90/10 of random simulations), so no site is calibrated on, but also no leave-site-out evaluation of the retrieval exists [O].
* Dark Zone applies the same inversion to S21, S22A, earlier K-transect S6 spectra and Sentinel-2 imagery 2017-2025; the S2 retrievals are **not validated against ground truth in the text read**; the paper says satellite posteriors are "much broader" [L]. No Himalayan or non-Greenland data.
* Which sites were excluded from calibration: the optical inputs are measured at S21 (algae, cryoconite) and on dust from all three sites; S21 spectra are then inverted, so **S21 is both a source of the absorber optics and a retrieval target** [O]. This is a circularity the authors do not discuss in the text read.

### 6. Parameter uncertainty, model discrepancy, calibration
* snowlaps: gradient descent from 20 random starts, best MAE kept (`snowlaps.py` l.~340-390). **Multi-start scatter is optimizer variability, not an uncertainty estimate.** No priors, no likelihood, no discrepancy term [A].
* Dark Zone: NUTS, 4 chains x 1000 draws (`mcmc_pipeline.py` l.212 `num_results=1000, num_burnin_steps=1000, num_chains=4`), Student-t likelihood nu = 3 (l.151), Beta priors (l.64-95), HDI intervals; convergence by visual inspection for field spectra and R-hat <= 1.006 mean for satellite pixels [measured]. This is a real posterior over the model's parameters. Calibration of those intervals against known truth is **not** reported in the text read [U].
* **Observed in the archive code [O]:** `sigma = 0.003` is a constant (l.61), whereas the paper says per-spectrum replicate SD "upper bound of 0.004". The paper says it adopts the **perfect-model assumption**; no model-discrepancy term (no GP, no inflation) exists in `mcmc_pipeline.py`. With sigma = 0.003 on 171 bands and a 1.8 % mean error, residuals of order 0.01 are 3 sigma; the t-likelihood damps rather than models this. The posterior is therefore conditional on fixed optics and an understated error and can be overconfident; **this is a reasoned expectation, not a demonstrated failure.**
* Other archive observations [O]: `emulator.py` l.46 hard-codes `np.arange(305, 2806, 10)` (251 values) while the network emits 440 outputs and the inversion slices `[:, 10:181]` of a 5 nm grid (350-1200 nm, consistent with the paper; the 305-2805 attribute appears stale and unused). The scaler's constants (read from raw bytes) give training ranges: inputs 1-2 `250-924.6` and `250-924.7` (density-like, kg m-3), dz1 `1e-5-0.15` m, dz2 `1e-5-2.5` m, LAP maxima `9.98e6`, `4.9e6`, `9.85e5`. The paper text lists algae to 5e5 cells/mL; the prior scale in the code for the third LAP is 5e5, i.e. the prior's hard upper limit is about half of the scaler maximum. Which column is which absorber, and the units, follow the paper's ordering and were **not** checked against the pickle's metadata.

### 7. Independent albedo, absorbed-energy and melt evaluation
* snowlaps: none in the code; the repository includes a Finse irradiance CSV for the forward app, not a melt evaluation [O].
* Dark Zone: darkening is expressed as **modelled** broadband-albedo differences (BBA with a particle set to zero), using a subarctic irradiance spectrum at SZA 50 degrees. The paper explicitly says it quantifies "direct radiative impact" and that indirect effects are attributed to the ice column; no observed albedo (e.g. pyranometer), absorbed-energy or melt comparison appears in the text read [L]. Statements such as "81 % of darkening is cryoconite + densification" are model outputs, not observations; the lab-measured inputs are observations.
* Therefore **no existing implementation evaluates independent albedo, absorbed energy or melt**; the Dark Zone paper's HCRF is a directional (nadir, 10 degree FOV) quantity treated as albedo under a Lambertian assumption (Methods) [A]; they restrict to < 1200 nm to limit anisotropy. We must not equate HCRF with albedo either; Baseline A fits hemispherical albedo directly and treats HCRF only through a free anisotropy factor.

## What this does to the research claim

1. The broad contribution "a Bayesian, physically based separation of algae, mineral and cryoconite darkening on bare ice with weathering-crust structure" **already exists and is released** (Dark Zone paper + archive). It cannot be our contribution. The earlier note that this was a blocking, unreadable dependency is resolved: methods and code are now read.
2. Their weak points are acknowledged or visible in their own text: cryoconite/crust degeneracy, single fixed optics per absorber, the perfect-model assumption, qualitative-only concentration checks, no leave-site-out test, no observed melt or absorbed-energy test. Those are **opportunities to test**, not proof that the method fails. We have not shown any of them matters quantitatively.
3. The remaining proposed hypothesis: that **experimentally supported variability in absorber chemistry** (a few measured or literature-bounded parameters, not free spectral shape) improves transfer to **independent sites**, compared with the Dark Zone configuration. For Greenland the authors' own measurements indicate that the mineral and algal spectra are fairly stable regionally, so a gain, if any, is more plausible for cryoconite optics, which are the least constrained input. That is a hypothesis; the outcome could be no gain.

## Hard limits on testing it with data we have

* Greenland optical validation only. No Himalayan data; Greenland cannot establish Himalayan transfer.
* The Zenodo spectra (105 HCRF, S21 2021 and S22A 2022) overlap geographically with data already used in this repository (S Greenland 2021) and carry **no concentration truth**. They allow a test of reconstruction and of cross-day/cross-site consistency of retrieved parameters, not of correctness of concentrations.
* S21 spectra cannot be an independent test of cryoconite optics because the cryoconite spectrum was derived at S21.
