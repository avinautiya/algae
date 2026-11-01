# Original logger header timezone check

A small original Keller logger export was acquired from the author-disclosed mirror with a single bounded HTTP206 range, bytes2110104–2121860 (11,757 bytes). Original member: `PlaineMorteLake_hydro_thermonynamics/data/CTD/2019-07-11/205144-10mH2O_11_07_2019-11_45_00.csv`. Uncompressed size62,488 bytes, ZIP CRC32 1742085422, SHA256 `f3c9e52b1105df1a342654072bea882828fde25b9fa85cb1b91bea2e05bb8fc9`. Both CRC32 and size verified before reading.

The actual original header identifies serial205144, device20514410mH2O, sourceDX5 filename, and `Record Start Time;11/07/2019 11:45`. The table columns are `No;;Date;Time;P1;T;ConRaw`, with first record11/07/2019 at11:45:00. There is **no UTC, GMT, CEST, DST, offset or timezone field in this header**. The archive filename likewise contains only a naive clock time.

This does not resolve the source-floor CSV explicitUTC header versus articleCEST labels. A logger clock may have been set to local orUTC time without recording that fact in the CSV. It would be unjustified to infer a timezone from the filename or the author's workstation export path. Preserve source-naive timestamps, and retain explicitclock-offset scenarios for comparisons requiring subdaily alignment. No source timestamps or benchmark results were modified by this inspection.

The complete small original sampleCSV is retained as `timezone_header_sample/205144-10mH2O_11_07_2019-11_45_00.csv`; it is an uncalibrated logger export and must not replace calibrated corrected temperature products.
