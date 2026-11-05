# New direct-stake validation and published-model comparison

The new campaign is Pirámide2022–2023, published by Melo-Velasco et al.2025 (https://doi.org/10.1029/2025JF008360), data https://doi.org/10.5281/zenodo.11581073 . The archive is CC BY4.0; its published MD5 and size pass. The candidate and simple controls were taken from the already fitted Pirámide-excluded folds, not retrained after inspecting new targets. PIRAMIDE_FREEZE.json and its SHA256 precede reading the new MAT bodies. No Pirámide observation, including the older2014–15campaign, contributes to those candidate fits.

## Eligibility and decoding

Three camera-read stakes have co-located measured debris depths0.50,0.45and0.24m. The source stores MATLAB5serialized timetables. A narrow read-only decoder extracts explicit table variable names, row counts and numeric columns and refuses unsupported structures; it does not execute the authors' MATLAB software. Native ablation changes equal0.917times distance changes to1.2e-16m. Despite a README description of 'mm w.e.', the stored ablation column is consistent with metres water equivalent and the authors' calibration code. Primary targets therefore use native camera distance differences directly in metres ice, independently checked against the converted column. This prevents a factor1000unit error.

Only the raw Hour.Pir weather table is used, never its filled/synthetic alternatives. All daily windows require finite03:00camera endpoints and complete measured T/SWIN/SWOUT. No missing weather or target values are filled, no negative increments are clipped and no quality filter uses predictive error. The03:00anchor follows the published night-image schedule. It was operationally specified after inspecting the MAT layout/values but before computing scores, so this is a frozen-coefficient external validation with audited endpoints, **not a pristine publicly preregistered blind experiment**.

MATLAB source datetimes have no timezone. They are treated as an unzoned source calendar, not asserted to be UTC. The first result serialization added a UTC suffix; corrected unzoned labels and all score-equivalent historical outputs are preserved. Numeric scores are identical after correction. The raw weather timestamps are25minutes past the hour. Both start-labelled and end-labelled hourly-mean interpretations are reported, not selected by accuracy. They have slightly different complete-day support. Source-code midday-albedo<=0.3screening is a secondary proxy; it retains all complete days here and does not independently prove no snow.

## First matched comparison

Zero-shift/start-labelled support includes57/46/9daily intervals at sites1/2/3,112stake-days in one shared campaign. End-labelled support includes58/47/9. Equal-site weighting is primary; pooled-day weighting is separately reported.

| Method | Equal-site MAE, mm ice/day | Pooled-day MAE, mm ice/day |
|---|---:|---:|
| Robust candidate |2.865|3.098|
| Previous candidate |3.085|2.988|
| Degree-day/thickness |3.061|2.962|
| Published KP1 2020 fixed effects |3.547|3.434|

The robust candidate beats KP1 on each of the three sites, but only site3improves against both the previous candidate and degree-day/thickness. Its equal-site advantage depends partly on giving nine site3days the same site-level weight as the longer site1/2records. A global superiority claim is inappropriate.

## Published models on identical lag-supported windows

Published equations and tables were read from the original Carenzo2016article XML, Winter-Billington2020Table4and DCG-MIP2026SupplementS21/S22. PDF tables were also visually checked. The2020KP1fixed effects are b0=0.62,b1=-1.46; the2026supplement gives more precise fixed effects0.6215,-1.4572. For an uncalibrated recipient, random effects are zero. KP1 uses the positive daily-mean temperature, not the sum of positive hourly temperatures. No recipient threshold tuning is used.

For DETI, the original published parameter functions are implemented, with constant original albedo0.13. The exact original threshold and fractional-lag execution conventions were not recoverable from author code. All four original/lagged gate and0/1°Cthreshold variants are labelled **equation reconstructions**, not exact author-code reproductions. Lags use nonnegative nearest integer hours. No variant is selected for a headline win. These published parameters were added after seeing candidate outcomes, so this is a benchmark expansion, not part of the original candidate freeze.

All methods on this second table share48/40/9days (97stake-days) after enforcing lagged-forcing coverage.

| Method | Equal-site MAE, mm ice/day | Pooled-day MAE, mm ice/day |
|---|---:|---:|
| Robust candidate |2.864|3.092|
| KP1,2026precise fixed effects |3.522|3.402|
| Degree-day/thickness |3.105|3.031|
| DETI reconstruction, original gate,0°C |2.703|2.877|
| DETI reconstruction, original gate,1°C |2.747|2.907|
| DETI reconstruction, lagged gate,0°C |2.644|2.826|
| DETI reconstruction, lagged gate,1°C |2.645|2.823|

The candidate improves against KP1 but loses to every tested DETI reconstruction. The original novelty/superiority gate remains unmet. This new direct-stake comparison must accompany the earlier Djankuat failure; neither is discarded. No subsequent model change can treat these data as untouched again.

## Verification and uncertainty

The audit independently checks native distance/ablation units,1808predictions,32primary score summaries,recipient exclusion and identical support of the expanded comparator table. The decoder's proprietary-table interpretation is narrowly validated against structural identities and published dates, not against an installed MATLAB runtime. Camera-reading uncertainty and source-clock interpretation are not fully quantified; shared meteorology, serial correlation, few site3days and one campaign limit inference. Differences of fractions of a millimetre per day must not be presented as statistically established gains. Meteorology is observed, not forecast.

## Readiness update

An external direct-ice-ablation campaign is now available with frozen recipient-excluded coefficients. It is not a new independent glacier replicate beyond Pirámide, not a pristine preregistered trial and not evidence of universal superiority. A published KP1comparison is implemented; exact DETI author-code equivalence remains open. Author list/order supplied by Avi: Avi first, Dr. Manish Mehta (Wadia Institute of Himalayan Geology) coauthor. Contributions, funding, conflicts and human approval remain unconfirmed and cannot be certified by software.
