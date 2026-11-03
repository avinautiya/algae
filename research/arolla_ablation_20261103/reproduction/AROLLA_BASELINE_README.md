# Published DEB_CF Arolla comparator — fetched before evaluation

The three small professional-output members were retrieved with bounded HTTPS byte ranges from [DCG-MIP outputs, record15754456](https://doi.org/10.5281/zenodo.15754456). Each decompressed member matches its central-directory CRC32 and recorded SHA256. The1.22GB archive was NOT downloaded, so no full-archive MD5 claim is made.

`arolla_baseline.load()` verifies every member hash, the13-column matrix and the1686-hour source-prescribed axis. Column1 is the actual6cm site thickness according to the archived output README. The official run-date workbook prescribes forcing indices1030–2715,4July2010 12:00 through12September2010 17:00. The site's [metadata](https://doi.org/10.5281/zenodo.3047649) declares localUTC+2. Output units are m water equivalent per step according to the official experiment protocol. No observed target was loaded or scored by this helper.

`map_window` requires an explicit interval-starting or interval-ending label interpretation. Both mappings are retained in `AROLLA_BASELINE_MAPPING.json`; neither was selected by its score. The example28August–11September midpoint labels are NOT verified visit times. Fractional-hour visits require a separate justified operator. The parent timing audit owns that uncertainty.

## Parameter provenance and scope

The original Arolla metadata gives6cm thickness measured near AWS/stakeH1D. Its roughness row explicitly borrows Miage parameters. The0.96W/m/K conductivity row has no source note; this is not evidence that it was independently measured at Arolla. No Arolla target-tuning of DEB_CF has been identified here, but complete parameter independence is not proven by that absence.

The [paper](https://tc.copernicus.org/articles/20/1895/2026/) describes a snow component in DEB_CF and negligible Arolla snowmelt over the whole intercomparison period. Site metadata lists a30August snow event. The exported matrix does not explicitly identify ice-only versus total melt in its header. This audit has NOT proved the matrix includes snowmelt, nor that it excludes it. Until that component definition is resolved, label the comparator as the published melt output under an unresolved component convention; do not declare strict native pure-ice matching from assumption alone.

Four baseline-helper tests pass. They verify clock interpretation and refusal behavior, not field accuracy. No model output was aggregated against the held-out315mm target before the parent's new freeze.
