# Original public Ogier data recovered, 1 November 2026

## Verified source path

The accepted-stage preprint itself discloses the author-hosted mirror `https://people.ee.ethz.ch/~werderm/pm-data/` in its code/data-availability section: https://tc.copernicus.org/preprints/tc-2021-151/tc-2021-151-manuscript-version2.pdf . The directory is public and lists the ZIP and README; it is a legitimate author-published mirror, not a workaround for the Research Collection access control. No additional requests were sent to the rate-limited Research Collection bitstream.

The mirror index and HTTP HEAD both succeeded. ZIP endpoint:

https://people.ee.ethz.ch/~werderm/pm-data/PlaineMorteLake_hydro_thermonynamics.zip

HEAD reports 337,848,233 bytes, Last-Modified 15 September 2021, ETag `142327a9-5cc08324a4e70`, Accept-Ranges bytes. Three bounded successful HTTP 206 requests recovered the central directory (131,072 bytes), first 8,192 bytes, and the 6,313,555-byte products block. ZIP directory has 454 entries. No archive code was executed. Each of 59 selected CSV/README members was checked against its directory-record CRC32 and uncompressed size; a SHA256 is recorded for every extracted file in SOURCE_MANIFEST.json. This is member verification, not a checksum of the full undownloaded archive.

The DOI record https://doi.org/10.3929/ethz-b-000504956 identifies the public dataset and MIT license. Mirror and repository code attribution belongs to Ogier and collaborators; no new ownership or relicensing is asserted here. The mirror archive is from 15 September 2021 while the deposited record is dated 26 October 2021. Their byte identity has NOT been established. Identify this package as author-mirror version, not checksum-equivalent deposited ZIP.

## Recovered measurements

- Original `data/Channel_elevation/P{1,2,3,4,5}_channel_elevation.csv` direct floor-elevation records with original uncertainty columns, plus P4/P5 water-level-cut elevations and station coordinates.
- Corrected measured temperature products `products/temp{1,2,3,5}_daily.csv` and `_hour.csv`; no graph digitization is needed for these values. No calibration offsets should be added again.
- Corrected measured stage products (1/2/3/5), discharge products, and observation-derived hydraulic products. These are processed observations with author selection/calibration, not raw unprocessed sensor measurements.
- Products also contain fitted/derived Nusselt estimates and author-computed standard-correlation predictions. Keep these separately labelled; they are not new independent observations or a valid independent benchmark if fitted to the same incision targets.

The fast complete selected set is `/private/tmp/ogier-author-original-data-20261101/` (59 members, 44,491,946 uncompressed bytes). The additive workspace package includes members below 1 MB; the several multi-MB `*_corr_raw.csv` high-frequency processed series remain in the private working copy, with hashes retained in the manifest. These `_corr_raw` filenames still denote corrected products, not original logger files.

## Newly resolved and remaining issues

P5 original floor CSV has explicit header `DateTime (UTC)` and per-observation ±0.1 m standard-deviation-error column during July. P5 floor elevations are 2731.74 m at 16 July 14:00, 2728.10 m at 25 July 10:28, and 2726.24 m at 30 July 14:31. These are more direct targets than daily top-water-surface cuts; retaining all recorded observations includes the slightly rising 30→31 July bed reading instead of silently filtering it.

The original UTC header must be reconciled with article plots labelled CEST and the notebook's naive timestamps before precise time integration. The same-file convention may be consistent, but that must be demonstrated, not guessed. Hourly/daily averaging support, sensor calibration uncertainty, reference-frame motion and thermal versus other geometric change remain measurement-model issues. Public data recovery does not by itself establish consistent model superiority.

## Reproduction

`recover_from_mirror.py --output DIRECTORY` makes the three same bounded public requests, checks archive version/size and per-member integrity, and safely extracts only selected CSV/README paths. It refuses incompatible source metadata instead of silently accepting a changed archive. The script has been syntax checked and its extraction function was rerun offline against the saved HTTP response blocks: all 59 members reproduce the exact original SHA256 and pass CRC32/size verification. The network recovery wrapper has not been re-run, to avoid redundant requests. Do not repeatedly run it against a server reporting rate limits.
