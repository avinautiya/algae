# Verification of available chemistry and matched thermal data

Verified by opening files in this repository or the downloaded archives, not from summaries. "Matched" means the same place and time as a spectrum or as the other quantity named; nearby is not matched.

## Chemistry
| Quantity | Where | Samples | Matched to a spectrum? | Matched to a thermal/energy record? |
|---|---|---|---|---|
| Algal cells mL-1 (S6 2017) | `cook2020_archive_cell_counts.csv` | 47 plots, 9 days | **yes** (same plot albedo and HCRF) | no |
| Algal cells (S Greenland 2021) | `chevrollier2023_tableS4.csv` | 18 of 20 spectra, 2 days | yes (HCRF) | no |
| Phenolic mass per cell | `williamson2020_phenolics_per_sample.csv` | 64 habitat samples, S6 2016 | **no** (habitat IDs, different year from the 2017 spectra) | no |
| Pigment per cell summary, MACs (phenol, chl a/b, carotenoids), HPLC peak spectra | `williamson2020_*` | 53 samples pooled / regressions | no | no |
| Purpurogallin vs Fe-purpurogallin absorbance | `prochazkova2025_fig4_*` (digitised) | 2 solutions | no | no |
| Cell biovolume per plot | Cook archive `Cell_dimensions.xlsx` | 40 plots with counts | yes, 2017 plots | no |
| Mineralogy, dust MAC | Dark Zone Table S1 (n = 17) and Fig. S2 (n = 10); Cook `Mineral_Abundances.csv` | not paired to spectra in files we hold | no | no |
| Cryoconite mass/concentration | none measured | - | - | - |
| Black carbon | none measured at the plots | - | - | - |

## Thermal / energy
| Record | Where | Content | Matched to spectra? |
|---|---|---|---|
| KAN_M hourly SEB forcing 2016-2019 | `data/empirical/promice_KAN_M_hour_JJA_2016_2019_seb.csv` | SW/LW up and down, albedo, air and surface temperature, turbulent fluxes, boom/stake/pressure-transducer heights | **no**: KAN_M (about 1270 m) is a different location from S6 (about 1000 m); the repository's own manifest treats S6 as one site and KAN as two other glaciers. Nearby stations are not colocated observations |
| PROMICE L3 hourly, 8 stations | `data/promice_raw/*_hour.csv` (main checkout, untracked) | as above, 2016-2023 | no spectra anywhere |
| PROMBIO station-days | `records/prombio_matching/station_days.csv` | ice-algal cells under the radiometer on a station-day, plus the station's albedo, SW and surface temperature | **cells-to-station-albedo is matched in time and place**; no spectrum, no pigment or mineral measurement. Development: 2021 (3 station-days) and 2023 (12 listed, 3 flagged for missing radiation/tilt data); about 30 sample rows. Final test (2024, 5 station-days) is withheld |
| Ablation | sensor heights in the KAN_M file | boom, stake and transducer heights disagree in 2019 per the manifest | no |

## Consequences
1. **No dataset pairs spectra, pigment/mineral chemistry and a thermal or melt record.** The only triplets that exist are pairwise: (spectrum, algal count) at S6 2017 and S Greenland 2021; (algal count, station broadband albedo, SEB) at PROMBIO.
2. A **chemistry extension therefore cannot be tested end-to-end** on melt or absorbed energy with data in hand. It can only be tested on spectra against algal counts, and only the algal term has a measured truth.
3. The measured chemistry that does vary (phenol per cell, CV 0.73) is not paired to any spectrum, so its effect on optics is not observable in these data.
4. Modelled melt, absorbed energy and observed melt must be reported as three different quantities; none is reported here.
