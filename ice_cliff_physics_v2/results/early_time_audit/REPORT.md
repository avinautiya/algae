# Post-study early-time numerical audit

This audit does not alter frozen production outputs. It evaluates cumulative ice heat at t=0.1 tau against the independent classical Robin series (512 terms), for all12 frozen boundary cases. The same FV physical inputs are used without tuning.

|Cells|Steps/tau|Relative excess ice heat range|
|---|---|---|
|64|4000|1.321% to 1.875%|
|64|16000|0.436% to 0.665%|
|128|4000|1.217% to 1.679%|
|128|16000|0.331% to 0.469%|

Both time and mesh errors matter. Increasing steps/tau4000→16000 reduces temporal error substantially; doubling cells64→128 reduces the remaining spatial component. Fine128/16000 remains0.331–0.469% high, so even these runs are not exact early-time data. The original64/400 production early-time error must remain reported at its larger measured level; a closed energy ledger is not evidence of transient accuracy.

The implicit Euler surface-flux integration is first order in time; the cell-centred spatial method approaches second order for this smooth linear case. Positive excess heat here is consistent with diffusion of the initially sharp thermal transient by backward Euler. This diagnosis is numerical, not a new physical mechanism or environmental validation.
