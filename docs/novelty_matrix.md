# Novelty matrix and ablations

What this project adds relative to the closest published work, and how each claimed addition is tested
on the same held-out data. The descriptions of other work are limited to what their abstracts state
(retrieved from Crossref on 2026-10-09). Anything beyond the abstracts is marked "not verified here".

## Closest work

| Work | What it does (from the abstract) | Optics of glacier algae | Evaluation |
|---|---|---|---|
| Cook et al. 2020, The Cryosphere 14:309, doi:10.5194/tc-14-309-2020 | Field spectroscopy + radiative transfer, UAV/satellite classification and runoff modelling; algae added 4.4–6.0 Gt of runoff from bare ice in SW Greenland in 2017 (10–13 % of the total); up to 26 % faster melting in high-biomass patches | BioSNICAR-GO (details not verified here) | Field spectra, classification |
| Williamson et al. 2020, PNAS, doi:10.1073/pnas.1918412117 | Photophysiology of glacier algae; secondary phenolic pigmentation about 11× chlorophyll a; captured UV/SW radiation repurposed for melt | Measured in vivo pigment absorption (the source of this project's `williamson2020` MACs) | Field and laboratory |
| Whicker et al. 2022, The Cryosphere 16:1197, doi:10.5194/tc-16-1197-2022 | SNICAR-ADv4: a physically based two-stream adding–doubling RT model for glacier ice (bubbles, SSA, light-absorbing constituents, overlying snow) | Algal representation not verified here | Not verified here |
| Whicker-Clarke et al. 2024, JGR Atmospheres, doi:10.1029/2023JD040241 | SNICAR-ADv4 in E3SM; bare-ice properties from satellite observations; effect on Greenland albedo and surface mass balance | Not verified here | Earth-system-model evaluation |
| Williamson & Tedstone 2026, Commun. Earth Environ., doi:10.1038/s43247-026-03758-8 | Greenland-wide glacier-algal bloom simulations, 2000–2022, quasi-Monte Carlo ensemble informed by sensitivity analysis; conditions for growth around the entire margin every year | Not verified here (a bloom/biomass model) | Not verified here |

## What is claimed as new here, and the test of each claim

| # | Claim | Test on identical held-out data | Status |
|---|---|---|---|
| N1 | Pigment optics from first principles (TD-DFT on the glucoside chromophore), calibrated only on laboratory spectra, transfer to field reflectance | `phase4/heldout.py`: `tddft_C`, `tddft_D` vs `measured_mac_C` and `tierA_empirical`, both folds, log score + forward reflectance score | see results below |
| N2 | An Fe-complexed fraction is needed to reproduce visible absorption | `tddft_D` vs `tddft_C` (ablation of the Fe term) | see results below |
| N3 | The physics retrieval beats statistical reflectance baselines at an unseen site | the physics models vs `band_ratio`, `ridge`, `climatology`, `literature_prior` | see results below |
| N4 | Algal melt increment from an SEB rather than "potential melt" | `phase4/seb.py`: paired on/off runs; SEB checked against PROMICE fluxes, surface temperature and three ablation records | conditional modelled estimates only: the increment/potential ratio is about 0.9 (χ = 0) or about 0.6 (χ = 0.3, subsurface shortwave), with no independent algal-melt observation; the 2019 ablation sensors disagree by 2.5× |
| N5 | A reusable surrogate with domain flags | `phase4/surrogate_api.py` + benchmark against direct BioSNICAR | code and unit test done; benchmark pending |

Results of N1–N3 are filled in from `phase4/results/heldout_v1/` once the run completes. If pigment-aware optics do
not improve held-out scores, that is reported as the result.
