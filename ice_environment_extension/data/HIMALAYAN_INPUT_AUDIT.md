# Himalayan inputs and what they can support

Checked on 1 November 2026. A public catalogue is not proof that every file has
been downloaded or that variables are co-located.

| Input | Source and access evidence | Supported use | Missing link |
|---|---|---|---|
| Terrain and cliff cross profiles | [Kraaijenbrink & Immerzeel v1.2.0](https://zenodo.org/records/13073856), DOI 10.5281/zenodo.13073856; selected members actually downloaded by bounded HTTP ranges | 391,911 observed corrected points, 11 profiles, 83 profile epochs across 9 cliffs, Lirung and Langtang | No bare-ice, water or debris masks; not co-located discharge/temperature |
| Survey errors | Same archive, `accuracy/survey-error-statistics.csv`, actually downloaded | Relative registration/error audit; non-reference vertical RMSE 0.1817–0.6194 m | Not local normal-retreat error, no covariance supplied here; reference zeros do not imply perfect terrain |
| Local meteorology | [NASA NSIDC HMA Langtang AWS](https://doi.org/10.5067/TNQV3XTMP9CP), product catalogue inspected | Air temperature, pressure, rain, humidity and wind for 22 Oct 2017–1 Nov 2018 | File-level ingestion not completed; radiation and site representativeness must be checked independently |
| Kyanging AWS | [ICIMOD DOI 10.26066/RDS.22464](https://rds.icimod.org/metadata/dcab6bc7-9dd0-4686-86ce-70664c8589ab), metadata inspected | Potential local atmospheric forcing | Actual schema/access and dates not yet verified |
| Earlier Himalayan meteorology | [NCAR CAMP](https://data.ucar.edu/dataset/camp-himalayas-automated-weather-station-raw-data), DOI 10.26023/5692-F00X-EF08, metadata inspected | Potential independent 30-min meteorological experiment, Jul–Sep 2001 | Does not coincide with these 2013–2018 cliff profiles; no modern cliff validation claim |
| Stream runoff and water temperature | Not supplied in the acquired profile CSV or the AWS variable list | Explicit imposed inputs or genuinely sourced stream measurements | Cannot derive these by relabelling air temperature or station precipitation |
| Obstruction and under-rock geometry | Not classified in the acquired point cloud CSV | Conditional transmission envelope with clearly stated assumed parameters | Hydraulic aperture, porosity and path topology need site measurements or independent resolved checks |
| Collapse | Not labelled in these files | Explicit sensitivity scenarios only after separate solid-ice accounting | No mechanical prediction or collapse validation is established |

## Reproduction

`audit_zenodo_profiles.py` verifies bounded tar headers. `fetch_profiles.py`
fetches only selected members and validates each Content-Range and length.
`profile_download_manifest.json` records local SHA256 values. The full source
archive MD5 is recorded but **not** verified by a partial download.

`summarize_profiles.py` processes every profile using 0.5 m bins, only points
within 0.5 m of its profile and at least 5 points per retained bin. Median
elevation and within-bin IQR are descriptive statistics. The IQR is not a
confidence interval. This processing does not classify a cliff or infer melt.

Dataset attribution: Philip D. A. Kraaijenbrink and Walter W. Immerzeel,
*Unmanned aerial system data of Lirung Glacier and Langtang Glacier for
2013–2018*, v1.2.0 (2025), CC-BY4.0. Associated paper:
[Spatial and temporal variability of the surface mass balance of debris-covered
glacier tongues](https://doi.org/10.1029/2024JF007935).
