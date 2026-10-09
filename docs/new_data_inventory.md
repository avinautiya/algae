# New-data inventory and frozen roles (frozen 2026-10-09, before any model fit or prediction on these data)

## Roles

| Dataset | Role | Unit of independence | Notes |
|---|---|---|---|
| PROMBIO 2021 (GEUS doi:10.22008/FK2/IJA2Y4 V2) | **development** | station-day | 3 station-days (QAS-U/M/L) |
| PROMBIO 2023 (same file) | **development** | station-day | 12 rows. KAN-L/KAN-M are **excluded**: the provider report (Anesio 2024, DCE TR314 §2) says the KAN samples "were not deemed of good quality enough to be considered". KAN-M also has no date. |
| PROMBIO 2024 (GEUS doi:10.22008/FK2/D0JC0P V1.0) | **final test, withheld** | station-day: KAN_L 08-23, QAS_U 08-27, QAS_L 09-01, QAS_M 09-01, TAS_A 09-09 | Not used for any calibration, prior, discrepancy, threshold or selection. A test station-day without usable radiation or imagery is **excluded from the test**, never moved to development. |
| UPE_U 2018 random grid (BAS PDC doi:10.5285/ab953cb8-…) and UAS products | candidate (support experiment) | one site-day | not acquired; classifier outputs are predictions, not labels |
| AVIRIS / EnMAP / PRISMA | candidate only | — | no scene with co-located labels identified |

## What has been seen before freezing (disclosure)

- **PROMBIO 2021–2023:** the full workbook was printed during this audit (development data).
- **PROMBIO 2024** (seen while importing and auditing the file structure):
  - the CSV header and its first three data rows (QAS-L randomized-grid samples 2024_0–2024_2: 4.00e4, 4.03e4, 6.78e4 ice-algae cells/mL);
  - counts by station/date/stratum/observation kind;
  - the notes of the ten "no counts" rows.
  - No other 2024 count value has been viewed; none has been used.

## Observation states (preserved separately, never merged)

| State | Meaning | Source examples |
|---|---|---|
| `observed` | finite count, including an exact 0 if reported as 0 | numeric cells/mL |
| `left_censored` | below the detection limit | `<1000` (detection limit 1000 cells/mL: 2 µl counted) |
| `missing` | no information | `NA`, blank |
| `missing_sample` | provider note says the tube had no or insufficient sample | 2024 "no counts" + "no sample in the tube", "almost no liquid", "missing" |
| `unresolved_annotation` | "no counts" without an explanatory note | 2024_69 |

## Counting and support facts (from the provider report)

- **Sampling:** top 2 cm of ice (ice screw in 2023; chisel or screw recorded per row in 2024), preserved in paraformaldehyde (pigment lost, cells countable).
- **Counting:** haemocytometer, two counts, 2 µl in total, so the detection limit is 1000 cells/mL and the Poisson count is N = 0.002 mL × concentration. Concentrations are per mL of melted sample.
- **Strata:**
  - "under radiometer": the sample is under the station's downward-looking radiometer, so it is associated with the station albedo footprint, not proven to represent it;
  - "average ice in the area" / "representative";
  - "dark ice": targeted, **not an unbiased area mean**;
  - "randomized grid" (2024): grid extent not given in the table.
- **Coordinates:** workbook positions and station GPS are **station positions**, not sample positions.

## Matching manifest

`agency_design/prombio_manifest.py` writes `records/prombio_matching/manifest.csv`: one row per biological observation, with source DOI/version/sha256, station-day, stratum, observation state, candidate Sentinel-2 acquisitions (±3 days), station radiation availability on the day, mismatch, QC and role.

## Disclosure added 2026-10-09, after the matching manifest
- The matching summary printed the **daily station albedo (Σ usr_cor / Σ dsr_cor)** of the 2024 test station-days: KAN_L 0.577, QAS_L 0.152, QAS_M 0.367; QAS_U and TAS_A have none.
- These are potential H5 targets. They were seen after roles were frozen and before any prediction or protocol choice that could depend on them.
- The H5 protocol (amendment A2) uses no albedo-based selection rule, so these values cannot influence inclusion.
