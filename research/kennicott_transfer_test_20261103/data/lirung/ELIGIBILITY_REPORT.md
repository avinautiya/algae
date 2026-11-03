# Lirung public-data eligibility audit

**Decision: not ready for a native paired melt/air benchmark.** Authentic hourly thermistor coverage is verified, but the released record contains no native weather file, no outcome schema was authorized for access, and the source does not yet supply an unambiguous thermistor-to-stake/time-interval binding. No melt predictions, parameter fits or native ablation values were accessed.

## Verified source and safe access

[Zenodo4446387](https://zenodo.org/records/4446387) lists two files: thermal50324bytes, MD5 `65d67be6eaf8e47aa3a492a17e0ac3e7`; ablation2026bytes, MD5 `973dd111fb0d340235b3792058cfe339`. The thermal file downloaded transiently for metadata parsing, matches published size/MD5, SHA256 `624e0cd94a6411d58ff0758aebb3407c1ac41297ffdd5442e0ae58caa1ce4adf`, and uses cp1252 encoding. No raw copy was saved. The audit exports headers, clocks, cadence and missingness only; numerical temperature values are never exported.

Automatic approval review rejected downloading/storing the ablation CSV because its raw bytes contain outcomes. A subsequent thermal attempt that printed unredacted sample temperatures was also rejected before execution. The safe metadata-only thermal audit was approved and completed. No attempt bypassed the rejected ablation action. Its local schema, dates, native measurement operator and uncertainty remain unverified. The2015 published abstract/method-search excerpts exposed aggregate melt context during source review; any future test must say native outcomes sealed, published aggregate context known, rather than fully blind.

## Primary papers

[Chand and Kayastha2018, source article](https://www.climatenepal.org.np/sites/default/files/doc_resources/Chand%20%26%20Kayastha%20et%20al.%2C%202018%20SEAR.pdf), Tables1–2 and methods, was accessible through indexed primary-source passages; a complete PDF fetch timed out. It gives hourly-profile depths and total debris thickness; monsoon site1 is10/20/30cm with38cm total, site2 is5/15/35cm with42cm total. Winter/pre-monsoon profiles have40cm total. Rock density2700kg/m³, heat capacity750J/kg/K and porosity.33 are assumed, not measured. Table2 has monsoon38/42cm stakes near the profile elevations but does not prove pit-level co-location or show40cm winter/pre-monsoon stakes. The source must resolve that binding before scoring; matching elevation alone is insufficient. Thermal conductivity derived from temperatures with assumed capacity is not an independent measured material parameter.

[Chand et al.2015](https://piahs.copernicus.org/articles/368/21/2015/) methods describe daily bamboo-stake measurements and an AWS with five-minute acquisition, hourly analysis. Those weather records are not files in this Zenodo release. Published daily measurement practice does not prove this released stake CSV contains daily rows. The paper also records pre-monsoon snowfall, so snow/phase coverage cannot be assumed away.

## Native thermistor coverage

All below are literal source-local labels with undocumented timezone. See `THERMAL_COVERAGE_ONLY.json` for counts and reproducible source hash.

|Group|Site/season|Clock span|Rows|Depth headers|Issue|
|---|---|---|---:|---|---|
|1|1monsoon|Sep22–Oct2 2013|264|10/20/20cm|Duplicate20cm header conflicts with paper30cm;36/55/13missing|
|2|2monsoon|Sep22–Oct2 2013|264|5/15/35cm|Complete hourly; no surface/ice-interface observation|
|3|1winter first|Dec1–Dec4 2013|96|0/20/40cm|Complete; depth changed after block|
|4|1winter second|Dec7–Dec12 2013|143|0/10/40cm|One7200s gap and one−82800s reversal; do not silently repair|
|5|2winter first|Dec1–Dec4 2013|96|0/20/40cm|Complete; depth changed after block|
|6|2winter second|Dec7–Dec11 2013|120|0/10/40cm|Complete|
|7|1pre-monsoon|Apr8–Apr18 2014|264|0/10/40cm|Complete; source snowfall requires phase handling|
|8|2pre-monsoon|Apr9–Apr18 2014|240|0/10/40cm|Complete; same phase caveat|

The eight blocks are not eight independent glaciers or seasons. The longest block is11days, shorter than the current14dayTRAIN+14dayevaluation contract. Three profile depths also do not support the existing five-thermistor/four-node inverse without a new discretization/observability analysis.

## Pre-outcome protocol proposal and gates

No final melt split can yet be frozen responsibly. A feature-defined *thermal-method development* split is possible for the complete site2monsoon block: September22 00:00→September29 00:00 past-window establishment, then September29 00:00→October2 23:00 later temperatures. This is not yet a paired melt test and does not make the remaining days statistically independent. A7-day adaptive model would have only about4days after initialization; compare a source-adjusted shorter window only if chosen on separate development data, not on Lirung outcomes.

For a new environmental test, freeze the following before native outcome access:

1. Exact stake ID bound to thermistor site, ice component, depth, measurement dates/time uncertainty and density convention; no matching by favorable rates.
2. Endpoint or daily-rate integration operator and an explicit target uncertainty source. If only seasonal average rates exist, score seasonal means without pretending daily targets are available.
3. Weather-source file with actual native coverage; otherwise restrict comparison to thermistor-conditioned heat retrieval and constant-rate/calendar baselines. Do not call a buried-temperature proxy native air PDD.
4. Complete3-sensor forward/inverse model, surface/basal extrapolation assumptions and phase/snow handling. Include fixed-capacity scenarios as assumptions; do not claim independent capacity verification.
5. Freeze model, clock policy, split, all comparator information and retained-family/abstention thresholds before outcomes. If the source binding does not overlap thermal coverage, stop rather than invent forcing or targets.

**Concrete next step:** acquire outcome-schema/date/co-location metadata through an explicitly approved data-only gate, plus the actual AWS release or author metadata. Until then this archive is a genuine Himalayan thermal source and a candidate melt archive, not a ready external validation benchmark. This audit does not improve or reverse any prior failed field comparison.

## Reproduce safe thermal audit

`python3 audit_thermal_metadata.py` downloads only the thermal file, verifies size/MD5, and writes metadata. It never requests the ablation URL. Existing approved network access is required. No external local module or virtual environment is required.
