# Additional data discovery and acquisition

Checked 2026-10-09. Public records are not automatically compatible training labels. Preserve site/year/date/footprint and biological-domain differences. None of the new records below has been used to train a model on this branch.

## Priority 1: PROMBIO + matched AWS + satellite observations

**PROMBIO 2024:** [GEUS V1, DOI 10.22008/FK2/D0JC0P](https://doi.org/10.22008/FK2/D0JC0P). CC0. Original 29 KB count CSV and dataset metadata acquired in `agency_design/data/`; provider MD5 verified (`d037bcf469ba7b935e17cac92053efc8`).

Audited 118 records: QAS-L 35, QAS-M 35, QAS-U 30, TAS-A 15, KAN-L 3. Sampling strata: randomized grid 75, under radiometer 31, dark ice 6, average-area ice 6. Dates: 2024-09-01 (70), 2024-08-27 (30), 2024-09-09 (15), 2024-08-23 (3). These are sample rows, NOT 118 independent site-days.

The table distinguishes ice/snow algae counts and qualitative mineral/cell ratios. **Ten of 118 records say "no counts" with blank numerical columns; only 108 contain numeric ice-algae observations.** Do not treat that annotation as measured zero. It has no latitude/longitude or measured spectrum columns. Preserve European decimal commas/scientific notation. `prombio.load_2024()` validates the checksum and preserves raw rows, explicit missing/censored observations and unknown positions. No synthetic coordinates or satellite features are assigned.

**PROMBIO 2021–2023:** [GEUS V2, DOI 10.22008/FK2/IJA2Y4](https://doi.org/10.22008/FK2/IJA2Y4). CC0. Original 21 KB workbook and full metadata acquired, provider MD5 verified. Workbook has station coordinates, radiometer/representative/dark strata, qualitative abiotic load and cells/mL, including `<1000` censored values and `NA`. Do not use min/max/average summary columns as additional independent observations. Coordinates may be station-level, not individual scraped-sample positions. Supporting [PROMBIO report](https://dce.au.dk/fileadmin/dce.au.dk/Udgivelser/Tekniske_rapporter_300-349/TR314.pdf) explains the program; read before ingestion.

**Validation route:** match these campaigns to the corresponding [PROMICE/GC-Net curated AWS](https://promice.org/download-data/) and original satellite acquisitions. First match station/year/time, radiometer footprint, maintenance, snow state, QC and sample coordinates. Radiometer sampling strata are especially useful but do not prove direct footprint coincidence. Randomized grids can support subpixel heterogeneity; dark-ice targeted samples must not be treated as unbiased regional averages. Keep a complete campaign/site-year withheld for final evaluation BEFORE model fitting.

## Priority 2: northwestern Greenland spatial-support experiment

[UPE_U 2018 counts](https://data.bas.ac.uk/full-record.php?id=GB/NERC/BAS/PDC/01289), DOI 10.5285/ab953cb8-8675-4a85-b561-add6ceba015f: **75 randomly distributed surface samples in a 250×250 m area**, one date (2018-07-26). Open Government Licence v3. Verified metadata; raw files not acquired on this branch. Strong for spatial heterogeneity and support matching, weak for temporal generalization.

The [Tedstone et al. 2020 data-availability section](https://tc.copernicus.org/articles/14/521/2020/) links complementary products:

- Processed UPE UAS: DOI 10.5285/2dd66461-94af-458f-a9d2-c24bb0bd0322.
- Raw UPE UAS: DOI 10.5285/a87b7897-354c-4435-a1bc-e6053e7569e0.
- Processed S6 UAS/trained classifier: DOI 10.5285/77ca631f-a3a4-4f26-bc90-57bb17baa6fc.
- Classified Sentinel-2: DOI 10.5285/8e0a573d-61a4-4a6f-9fca-fc34cbd5fb45.

These are source-linked candidates, not all independently downloaded/verified. Existing classifier outputs are model predictions, not new biological ground truth. Audit S6 overlap with the repo before evaluation. Counts and co-temporal UAS footprints may support an independent pixel-support test; exact pairing remains to be checked.

## Priority 3: hyperspectral inputs for identifiability

[EnMAP campaign portal](https://www.enmap.org/data_tools/flights/) provides airborne campaign metadata and free Creative Commons data; matching glacier scenes/field labels must be checked. [EnMAP science plan](https://www.enmap.org/data/doc/Science_Plan_EnMAP_2022_final.pdf) discusses glacier algae retrieval with imaging spectroscopy.

[NASA/JPL AVIRIS data access](https://aviris.jpl.nasa.gov/data/get_aviris_data.html) is a portal lead, not evidence that any chosen flight line has glacier-algae labels. Select exact scene IDs, dates, processing level and co-located observations before downloading large cubes. JPL's [Niklas Bohn research profile](https://www.jpl.nasa.gov/site/research/bohn/) documents existing coupled atmosphere/snow retrieval and uncertainty work: this is prior art to compare against, not a novel capability to claim first.

[NASA HLS algorithms](https://hls.gsfc.nasa.gov/algorithms/) document atmosphere, cloud/shadow, BRDF, bandpass and gridding treatment. HLS reflectance can provide an independent sensor/product sensitivity test; it is not ground truth and is not interchangeable with ground HCRF or algae labels. Check native red-edge/SWIR availability rather than restricting every experiment to four bands.

## Supporting sources, not new glacier-algae labels

[GA_BLOOM code/source data](https://zenodo.org/records/20138073): download small figure-source tables/code first, not the 1.3 GB archive indiscriminately. Reuse for prior-art replication and ecology scenarios. Model trajectories are not supervised true abundance; aggregated literature observations may duplicate training records already in this repo.

[Antarctic spectral library](https://doi.org/10.7488/ds/7720), linked by [Walshaw et al. 2024](https://www.nature.com/articles/s41561-024-01492-4): green snow algae and other vegetation/lichen spectra can help negative controls and domain-shift tests. Do not pool them with brown glacier-ice algae as the same pigment/biology.

[GFZ metabolic/taxonomic preprint record](https://zenodo.org/records/17775183): currently the visible files are PDFs, not a verified machine-readable training table. Useful biological context, not a promised new supervised dataset.

## Ingestion gate

For every candidate require: original DOI/version/license/checksum; observation versus modeled label; dates/coordinates/CRS and spatial support; target units and counting protocol; censoring/missingness; sampling design; overlapping published records; and a predefined training/validation role. Do not count repeated samples or derived summaries as independent site-days. No final score or agency-readiness claim follows from dataset acquisition alone.
