# Authorized stake-schema addendum

Explicit user permission to download the public stake CSV was relayed by the parent after the earlier automatic-review refusal. The new download was approved. The2026byte file matches MD5 `973dd111fb0d340235b3792058cfe339`; SHA256 is `700d255c50b7d4fd22cb9563fbb4a372198cf64d7d91717af9992243f31ba382`. The raw file is stored locally, but numerical ablation magnitudes have not been inspected or exported. `STAKE_SCHEMA_MASKED.json` masks every numeric cell; `STAKE_METADATA_ONLY.json` selectively exports depth, elevation, position and period, never the ablation column.

## What the actual schema resolves

The file has three seasonal blocks, not a daily series: September22–October3 2013, November29–December12 2013, April7–April19 2014. Each stake has one scalar under **Ablation in cm**, plus debris thickness, elevation, latitude and longitude. There are seven monsoon rows including dirty ice, seven winter rows, and eight pre-monsoon rows including dirty ice. There are no stake IDs, visit times, daily endpoints, statistical uncertainty, ice-density convention or explicit cm/day operator. Do not reinterpret the scalars as daily rates based on their magnitudes. Published source results discuss rates, but that does not resolve the released CSV's incomplete unit metadata before outcome access.

## Spatial pairing

Monsoon native38cm stake (row8,4093m) is near source site1 profile (d38cm,4093m);42cm stake (row9,4154m) is near source site2 profile (d42cm,4156m). Literal coordinates place them approximately6m and13m from their respective thermistor positions, so these are nearby-site pairings, not verified same-pit measurements. Meter estimates assume the malformed typography in thermal metadata denotes conventional degrees/minutes/seconds. No position uncertainty is supplied. This is usable as an explicitly spatial forcing/geometry assumption, not proof of exact co-location.

Winter/pre-monsoon native stake depths do not include40cm, the stated profile-layer thickness. Several winter depths differ from source2018Table2; retain native-file provenance and do not silently merge rows based on the table. The thermal file's site2 coordinate/elevation also differs between seasons, so do not claim a single fixed pit across the year.

## Temporal pairing

Monsoon profiles cover September22 00:00 through October2 23:00. This could cover a nominal September22 00:00→October3 00:00 interval with a declared hourly left-hold final bin, but visit times are absent. A date-only visit-time uncertainty envelope cannot be fully supported for all possible October3 visit times by the thermal record. Winter profiles startDecember1, missing at least the start of the November29 target period, and have a two-day split. Pre-monsoon profiles startApril8/9, later than the April7 target. Therefore neither later seasonal target has fully native thermal forcing over its declared interval.

## Feasible next test and limits

A prospective *native-schema-sealed external spatial comparison* could reserve monsoon38cm as one local calibration point and42cm as one evaluation point, using only their native thermistor histories. It requires resolving cm versus cm/day, visit-time assumptions and the monsoon site1's duplicate20cm depth label against source30cm, plus predeclared missing-profile interpolation. Compare physical retrieval, buried-temperature dose and constant-rate transfer using the same allowed inputs. There is no released air series to support native air-PDD or a complete surface-energy model. A calibration point and one evaluation point cannot establish consistent superiority, uncertainty calibration, or transfer across seasons.

A fully chronological local-train/local-test benchmark is **not possible from this release alone**, because individual daily ablation measurements are not released. Using the published seasonal scalar as both training and later-test target would leak the same observation. Restricting thermal integration to a well-covered subperiod cannot create the corresponding unobserved stake increment.

**Decision:** keep numerical magnitudes sealed. Obtain operator/visit metadata or daily stake/AWS release before a production melt campaign. This is an actionable source-specific data requirement, not a claim that environmental models require perfect zero-uncertainty observations. The source capacity assumptions can be propagated; the missing target operator and interval cannot be inferred by choosing whichever interpretation gives smaller error.
