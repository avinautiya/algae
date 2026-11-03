# DCG-MIP benchmark red-team and bounded output retrieval

## Verdict

**Conditional GO for a historical professional-baseline reproduction; NO-GO for claiming an independent pure-ice melt test today.** No model was trained and no new accuracy score or favorable result was manufactured. The professional output bytes are accessible; remaining blockers are observation/calibration contracts, not an absence of public files.

## New facts requiring correction to the next-task specification

The primary DCG-MIP paper evaluates final cumulative percentage melt error in its main Table3, not daily MAE alone (§4.3.1). It explicitly says Changri Nup's supplied debris properties were previously optimized against observed melt (§6.1.5). Therefore those supplied values are not independent test parameters. Its DEB_CF sensitivity analysis reports substantial snowmelt at this site (§6.2), so an ice-only target must not be silently substituted for the paper's comparison. The experiment was open and some empirical models used the provided observations for calibration. These are disclosed design features, not evidence that the paper is invalid. Source: https://doi.org/10.5194/tc-20-1895-2026 .

The paper identifies DEB_CF, ROU15 and d2EB as strong energy-balance comparators. Calibrated empirical approaches must remain visible, but their calibration information differs from an independently transferred candidate. An appropriate fair comparison needs matched calibration budgets and an explicit distinction between historical reproduction and independent transfer.

## Actual output archive inspection

The concept DOI15754455 resolves to record15754456. Its archive is1,222,386,989bytes, with published MD5 `4e5e66f44d23fdc0829995ca876afd37`. We did NOT download the entire archive or verify that full-file MD5. Only the small README, final128KiB, central directory1,306,749bytes and three selected members were retrieved by bounded206 HTTP ranges. The directory has15,624entries.

`fetch_small_outputs.py` checks requested Content-Range, byte count, decompressed size and member CRC32, then records SHA256. It retrieves:

* `DEBCF/CN/DEBCF_CN_EXP1_EXP2.csv`:2399rows ×13columns,271,678uncompressedbytes.
* `Dates_Rows/Dates_hourly_CN.csv`:2399naive date labels,16Apr201601:00 through24Jul201623:00.
* `Dates_Rows/Dates_daily_CN.csv`:100labels at12:00,16Apr through24Jul2016.

The output README specifies the FIRST matrix column as actual site thickness (Changri Nup10cm); columns2–13 are the sensitivity-thickness experiments. A sensitivity column must not be selected because it fits better. The archive metadata counts14models/17versions, while the final paper has15models/17approaches; use explicit acronym/version correspondence rather than assuming the metadata counts define the manuscript cohort.

## Why I did not report a baseline score

The selected output CSV lacks a units/header declaration. Its actual-thickness column sums to0.7186453836834981 in its native, unresolved units. This is NOT a claimed melt total in any harmonized unit. The matrix fits the hourly axis by row count, but the hourly interval-boundary/time-zone convention is not explicit in the downloaded README. The daily axis is labeled noon while the observation archive uses24:00 labels. Row-count and calendar overlap alone cannot certify matching interval support or unit conversions. The original AWS metadata and observation audit also retain sensor migration, snow and footprint issues.

The prior loader provides the original daily SR50 lowering and timestamp metadata, not an independently verified pure-ice mass observation. The paper can be reproduced on its declared processed target once its exact conversion/aggregation is reconstructed. That reproduction should then be labeled historical and calibrated as appropriate; it does not prove independent field accuracy of a newly transferred model.

## Implemented gate and tests

`eligibility_v2.py` separates `published_baseline_reproduction` from `independent_pure_ice_melt_test`. Six tests check metric mismatch, optimized properties, string booleans, snow-target mismatch, missing model mapping, and why baseline reproduction is not independence. Passing is metadata eligibility only, not certification of assertions or a performance result.

`CURRENT_GATE.json` deliberately leaves unresolved dates, units, observational support and calibration independence unmet. `SELECTED_OUTPUTS_MANIFEST.json` records the actual retrieved members and byte hashes. `OUTPUT_RECORD.json` and `OUTPUT_README.txt` preserve original public metadata. `ZIP_MEMBERS.json` provides filenames/sizes/CRC/offsets so another investigator can reproduce the small retrieval without a1.22GB download.

## Required next bounded action

Recover the exact experiment input/analysis protocol (DOI10.5281/zenodo.3134763) and supplement, then document the unit, density, snow treatment, period boundaries and observation aggregation used for one Table3 cell. Reproduce that cell without selecting among interpretations by goodness of fit. If an ambiguity remains, keep all predeclared interpretations labeled conditional rather than choosing the published-looking value. Independent transfer then requires parameters and training decisions not informed by the evaluation melt, and a separately frozen site/period contract. Nothing here establishes novelty or coupled cliff validation.
