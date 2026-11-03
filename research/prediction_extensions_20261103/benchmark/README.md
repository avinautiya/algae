# Public benchmark eligibility and verified loader

## Useful new comparison route

The2026 DCG-MIP study compares15 debris-melt models at9sites, including Himalayan sites, and links meteorological/validation archives and submitted model outputs. This supplies an established comparison framework; it does not establish that our candidate is novel or more accurate. Primary paper: https://tc.copernicus.org/articles/20/1895/2026/tc-20-1895-2026.html . Submitted outputs: https://doi.org/10.5281/zenodo.15754455 . Site links include Changri Nup https://doi.org/10.5281/zenodo.3048780 and Lirung https://doi.org/10.5281/zenodo.3050327 . Only the Changri archive was downloaded and inspected in this task; other record alignment is not certified.

## Actual verified data

Downloaded Changri Nup `2_CN-5470_20152017.zip`,4.3MB; published MD5 `4797deeda5c93a1895f08efb42183086` matches. The archive has hourly meteorology, two interface flux/temperature profiles, daily processed surface lowering and a full metadata document. `loader.py` executed successfully:10,813hourly rows and450daily rows. Eight tests pass, including explicit NepalUTC+5:45 conversion, source24:00 end-of-day labels, a scientific-decimal-comma source literal, missingness, ordering and prohibited-feature separation. The loader preserves native NaN as null, never zero-fills, and retains signed flux. It does not train a model or calculate performance.

`CHANGRI_ALIGNED_AUDIT.json` contains aligned input rows, raw-target interpretations, memberSHA256s and open gates. `CHANGRI_ARCHIVE_MANIFEST.json` records source bytes. The metadata text is extracted for reproducibility; no visual-layout or OCR claim is made. The public source record remains the authority.

## Critical observation distinctions

The daily `ABL` is thresholded SR50 surface displacement, not directly measured pure ice melt. Changes below5mm are set tozero. Snow, settling, differing sensor footprints and source clock/averaging conventions require explicit treatment. Independent stakes and dGPS are listed as available upon request, not in this archive. The interface heat-flux plates have an independent instrument but sign/orientation and sensor migration must be audited before scoring. The source declares initial probe depths and a possible later migration; no depth history is invented.

Metadata says `SNOW` and changing sensor heights derive from the same SR50 underlying the target. **Shared source does not automatically imply invalid leakage.** Start-of-interval state available at issue time may be a legitimate conditional input. End-of-day values derived from the future target interval may leak information for a forecast. Our strict independent-source protocol conservatively excludes `SNOW,T_z,RH_z,FF_z,DIST,HEIGHT,ABL,ACC` from model features until an issue-time/provenance audit distinguishes those cases. This changes the input contract and is not an accusation that the published intercomparison is invalid.

## Frozen-before-fit protocol to implement after eligibility is resolved

1. Fix target and units: processed daily displacement, independently measured ice ablation, or signed interface heat flux. Do not rename one as another. Declare conditional hindcast versus future-forcing forecast.
2. Pin site archives and all professional output hashes. Audit time windows, uncertainty, snow screening, footprint and calibration provenance. Do not infer a flux sign by selecting whichever agrees with the model.
3. Assess inspected Swiss/Kennicott only as development or retrospective tests. Their previous negative outcomes stay visible. These sources are one glacier's repeated profiles, not independent glacier replications.
4. Reserve complete sites for external evaluation. Select hyperparameters by leave-one-training-site-out validation only. No random hourly split. If too few independent training sites exist, refuse to rank flexible hybrids as externally validated.
5. Compare: established temperature-index/constant-rate; a reproduced energy-balance/conduction baseline; training-only ridge; and a physics-plus-residual ridge hybrid. Standardization, penalty choice, conductivity and residual fitting use training sites only. Target-derived meteorological state must satisfy the declared issue-time contract.
6. Predeclare primary equal-site MAE and cumulative signed bias, plus all-site win counts and source-scenario sensitivity. Keep model ranking separate by target. Report intervals only if justified by independent groups. No test-set family selection, best-scenario selection or discarded failures.
7. Accessing public targets does not make a test automatically blind. Confirm project-wide inspection history and declare any inspected comparisons retrospective. A new external site may offer a prospective analysis only after its outputs and selection rules are frozen.

## Current go/no-go

**GO:** reproducible source audit and preparation for an established professional component benchmark. **NO-GO:** claiming current public rows validate coupled cliff recovery, a satellite model, pure ice melt, field superiority or novelty. The new archives offer a route to test a debris-melt component; they do not fill the Swiss stream-clearance or cliff-contact inputs.

Run `python3 -m unittest -v test_loader`, then `python3 loader.py changri_nup.zip NEW.json`. Never overwrite the source audit. No performance improvement has been computed or promised.
