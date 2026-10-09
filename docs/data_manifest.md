# Data manifest

sha256 prefixes (16 hex) refer to the files as stored in this repository or in the pinned BioSNICAR checkout (revision fe74eeef2999fd285f3dc06c47b0f0e1a36a17a6, MIT licence). Full provenance and citations are in `data/empirical/SOURCES.md`.

**Use codes:**
- **F** — fitted on: model parameters, priors or calibration.
- **S** — used in selecting hyper-parameters or variants.
- **E** — evaluation target.
- **—** — not yet used.

## Field and laboratory data

| File | Content | Support / time | Uncertainty and QC | Licence | sha256 | Uses so far |
|---|---|---|---|---|---|---|
| `data/empirical/cook2020_archive_cell_counts.csv` | S6 2017 algal counts, cells mL⁻¹ meltwater, cells counted | 47 plots, 13–25 Jul 2017, about 1 m², top about 2 cm | Poisson 1/√N; 5 zero counts; zero-count counted volume not recorded (assumed 0.016 mL) | Zenodo 10.5281/zenodo.3564501, CC BY | 41128a21… | E (heldout_v2 S6 fold); F (heldout_v2 training, reverse fold) |
| `data/empirical/cook2020_archive_hcrf.csv` | S6 2017 nadir HCRF spectra, 350–2499 nm | same plots | ASD; 20_7_SB4 empty | same | 521da52f… | S (σ and radius prior, field_validation and heldout_v2); E (heldout_v2 forward test) |
| biosnicar `Albedo_master.csv` | S6 2017 **hemispherical spectral albedo**, 350–2500 nm, 87 plots (47 with counts) | same campaign | cosine receptor; accuracy about ±0.01–0.02 (not stated per plot) | MIT (biosnicar-py) | b4f645cb… | **H1 target**; NB: indirectly used via the k prior (ARF = HCRF/albedo, `biosnicar_field_ARF.csv`), see deviation D3 |
| biosnicar `ARF_master.csv` / `data/empirical/biosnicar_field_ARF.csv` | HCRF/albedo ratio spectra, 51 plots | S6 2017 | derived quantity | MIT | a85768fd… / 1878bf6c… | F (k prior N(0.898, 0.175) in all retrievals) |
| biosnicar `HCRF_master_16171819.csv` | HCRF 2016–2019, 351 columns | S6 and others | not yet audited | MIT | 42fdab2e… | — |
| `data/empirical/cook2020_field_metadata.csv` = biosnicar `Spectra_Metadata.csv` | published BioSNICAR-GO retrievals for 31 plots | — | model output, not observation | MIT | bbe18a46… | E (reference only; never a truth target) |
| `data/empirical/chevrollier2023_tableS4.csv` + `chevrollier2023_hcrf.csv` | S Greenland 2021 counts (18) + ASD HCRF | 5–6 Aug 2021, 2 sampling days | number counted not published (log-normal SD 0.10 dex assumed) | CC BY 4.0 | 3f8c04d4… / 9f933e52… | E (heldout_v2 primary fold); F (heldout_v2 reverse fold) |
| `data/empirical/tedstone2020_s6_2017_sample_locations.csv` | GPS of 20 S6 plots | UTM 23N | handheld GPS, about 3–5 m | UK PDC | 67d94ae7… | E (satellite validation, legacy) |
| `data/empirical/williamson2020_*` | S6 2016 counts, biovolume, pigments per cell, MACs, HPLC spectra, extract MAC | S6 2016 | per-file in SOURCES.md | per source | see SOURCES.md | F (abundance prior, cell size/concentration, calibration of TD-DFT, measured-MAC optics) |
| `data/empirical/prochazkova2025_fig4_*` | Fe–purpurogallin absorbance (digitised) | lab | digitising error not quantified | publication | 650d62e2… | F (tier D Fe increment) |
| `data/empirical/ice_ssa_measurements.csv` | bubbly-ice SSA | Greenland and lab | — | publications | efdd290e… | F (radius prior, measured_ssa option) |
| `data/empirical/esa_s2_srf_TN-15-0007_v4.0.csv` | Sentinel-2 spectral response functions | — | official | Copernicus | 175e0123… | F (band integration) |

## PROMICE (GEUS Dataverse doi:10.22008/FK2/IW73UU, dataset version 40, CC BY 4.0)

| File | Content | Time convention | Known issues | sha256 | Uses so far |
|---|---|---|---|---|---|
| `data/empirical/promice_KAN_M_hour_JJA_radiation.csv` | KAN_M SW↓, cloud cover, T | hourly means, time stamp = start of hour (readme) | 2019: no `dsr_cor` | 4a814a7a… | F (clear-sky transmissivity T = 0.919); E (multi-scene daily SW, legacy) |
| `data/empirical/promice_KAN_M_hour_JJA_2016_2019_seb.csv` | KAN_M full SEB forcing + heights | as above | 2019 `tilt_y` missing → uncorrected radiation and albedo; pressure transducer vs stake disagree 2.5× in 2019 | e6300835… | F (χ sensitivity fitted on 2016–2018); E (SEB validation against ablation) |
| `data/promice_aux/KAN_M_hour_MJJAS_2016_2019_qc.csv` | tilt, sensor heights, QC columns | as above | — | 5462fdd9… | E (sensor audit only) |
| `data/promice_aux/AWS_data_readme.pdf` | GEUS AWS readme (2024 update) | — | — | dea361c4… | documentation |
| `data/empirical/promice_KAN_{L,M,U}_day_2019.csv` | daily T, altitude | daily | — | 2426cb0f… / 07dfe514… / e37be365… | F (PDD maps, diagnostic only) |
| `data/promice_raw/KAN_L_hour.csv`, `KAN_M_hour.csv` (untracked, read-only; sha256 in `data/promice_raw/SHA256SUMS`: KAN_L 2a6fee92…, KAN_M 4665d74e…) | PROMICE L3 hourly, GEUS THREDDS `aws/l3sites/csv/hour` (live product; last-modified 2026-10-09 16:19 UTC, downloaded 2026-10-09; Dataverse doi:10.22008/FK2/IW73UU) | KAN_L 67.095 N −49.958 E 651 m; KAN_M 67.068 N −48.844 E 1268 m (median 2016–2023) | hourly means, hour-start; `dsr_cor` coverage JJA: KAN_L 0.82–0.87 (2016–2019), 0.42–0.74 (2020–2022); KAN_M 0 in 2019 and 2021, 0.27–0.43 in 2020/2022 (tilt) | CC-BY 4.0 (GEUS) | see SHA256SUMS | **H2/H4 targets** (SW↑, albedo: targets only, never inputs); SW↓ as forcing |

## Site transfer assumptions

- **KAN_M (67.067 N, −48.836 E, 1270 m) is 20–25 km from S6 (≈ 1000 m).** Using KAN_M meteorology for the S6 surface assumes:
  - the same SW↓;
  - LW↓, T and q adjusted only through the station measurement, with no lapse-rate correction.
  - The KAN_L–KAN_M difference (670 vs 1270 m) brackets the elevation effect. Once KAN_L hourly data are added, the S6 SEB increments are re-run with each station's forcing as a transfer-uncertainty range.
- **Station radiometer footprint (tens of m²) vs Sentinel-2 pixel (100–400 m²).** These are not co-located supports.

## Satellite

- Sentinel-2 L2A from the Element 84 earth-search v1 collection `sentinel-2-l2a` (COGs in s3://sentinel-cogs). Per-scene provenance (asset hrefs, processing baseline, offset rule, tilt, sun/view geometry, grid) is written by the scene-interpretation pipeline (`docs/scene_product_schema.md`).
- Original assets are read remotely and never modified.

- **M2 conversion source.** Naegeli et al. (2017) could not be obtained verbatim on 2026-10-09: mdpi.com returned 403 and the ZORA mirror is behind bot protection. M2 is therefore the protocol's "simple empirical" fallback (`phase4/h2h4_station.py`).
- **Sentinel-2 L2A availability.** Element 84 earth-search returns 0 acquisitions over KAN_L in June–August 2016; see protocol amendment A1.
