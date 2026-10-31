# Ogier et al. (2021) data collection — retrieval and readiness (no model fitted or scored)

Record: ETH Research Collection, "Drainage of an ice-dammed lake through a supraglacial stream: hydraulics and thermodynamics – Data Collection", DOI 10.3929/ethz-b-000504956, issued 2021-10-26, MIT licence (DSpace metadata via `/server/api/core/items/27f98ac0-…`).
Access log: `fetch.log`. The HTML landing page and most REST calls returned **HTTP 429 (rate limited)** for ~20 min; after back-off the bundle list, bitstream list and the archive download returned 200. The separate README bitstream stayed at 429 (the saved 1560-byte file is the rate-limit page, not data); the archive contains its own README.txt, which was read.
**Retrieved:** `PlaineMorteLake_hydro_thermonynamics.zip`, 337,848,233 bytes, MD5 c446a48f…b15d = the repository's checksum (verified). 454 entries, listing in `zip_listing.json`. Raw archive is untracked (`raw/` is .gitignored); not executed.

| Group (inside zip) | Files | What it is (from names/headers/README only) |
|---|---|---|
| data/CTD (72 csv) | Keller CTD loggers, 1 s sampling, semicolon-CSV with decimal commas, pressure (mbar) + two temperature channels + raw conductivity; first rows show a logger at 34 °C before deployment | **measured, point sensors** (near-bed per the article); 2019-07-04…2019-10-08 |
| data/DC22, SR50 | pressure logger and SR50 sonic ranger series | measured (stage / surface-height proxies); sensor-to-station mapping needs the notebooks |
| data/lake_data | lake level, volume–time series 2012–2019, hypsometry, discharge_2019.csv (lake and glacier discharge, derived), 442 MB DEM GeoTIFF | derived discharge; lake geometry |
| data/Salt_injections | 3 rows | dilution-gauging events (calibrate discharge), few |
| data/Channel_elevation | P1–P5 channel-floor elevation at ≤ 6 dates per station, `Std_deviation error` column = **0.5 m**; stations 464 m apart overall; P4/P5 "levels" file headed "DateTime (hypothesis. Only the time of 30 July was true)" | **geometry is 5 stations**, ±0.5 m (the column name says 'standard deviation' but the article's meaning is not established — not treated as Gaussian) |
| products (50 csv, 19 plots) | hourly/daily discharge (with std), stage1/2/3/5, temp1/2/3/5, Reynolds, Manning, friction factor, **Nusselt** (DB standard, Lunardini, Vincent, Gnielinski, SCR, MR) | **derived** from the above; Nusselt shares stage/discharge/temperature inputs with the incision proxies and is not an independent observation |
| PlaineMorteLake_hydro_thermodynamics.jl (+ .git, 16 notebooks) | processing code | needed to reproduce products and to identify exact definitions |

## What supports a forward test and what does not
* **Forcing:** hourly channel discharge and lake discharge exist (derived, with std files) — usable as *prescribed forcing*, labelled as such. Inlet temperature is not separately defined: `temp*` are near-bed point sensors at P1/P2/P3/P5 — they must not be compared with a model bulk temperature or used as inlet bulk temperature without an observation operator (not built).
* **Geometry:** five stations of floor elevation with ±0.5 m at a handful of dates and a 50 cm DEM of the lake region; no cross-sections or longitudinal profile at forecast resolution. Observed daily incision is centimetres to decimetres, i.e. at or below the stated vertical uncertainty → **a forward test of bed lowering is not supported**. Later elevations cannot be inputs when they are the target.
* **Plausible, narrower test:** reach-integrated cooling (temperature drop between stations) and hourly Q → requires bulk-vs-near-bed operator, sensor timing/calibration checks (not done), and channel geometry assumptions.
* **Not confirmed:** time zone/clock conventions of each logger, calibration files, quality flags, completeness (notebooks not run; products contain NaN for early hours), whether P-stations match the CTD serials (Salt_injections hints 205309/207265 at P2/P1).
Verdict: **inventory complete; data not yet "ready" for a forward bed-lowering prediction.** Heat-export/cooling comparison is the only candidate and needs the operator work above.
