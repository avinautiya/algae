# Arolla native ice-ablation comparison, frozen before predictions

This is a retrospective temporal holdout on independently measured ICE ablation, not a blind experiment. Source metadata and target values were audited before writing the code. The second interval never fits parameters or selects models. This does not validate the stream/cliff solver or assert a new physical law.

## Target and data

Source DOI10.5281/zenodo.3047649, source workbook MELT[mmw.e.], ice density900kg/m³. Four stakes H1D/G1D/F1D/E1D have3 dated visits. Fit each stake's first interval ONLY, predict its second interval. H1D is the PRIMARY location because the source identifies it as closest to AWS and provides its measured6cm±10%debris thickness. The other3 are secondary: their thicknesses are unknown and weather is measured at the AWS, not individually. They are not independent weather events or glaciers.

The held-out period is28August–11September2010. Initial calibration dates are28July forH1D/G1D and30July forF1D/E1D. All visit hours are unknown. Primary nominal hours are12:00local at each visit, declared before predictions. Propagate all27 joint combinations of0/12/24h for the3visits. The shared middle visit uses the SAME hour in calibration and evaluation. Time zone is source-declared UTC+02.

Hourly values are recorded means. Nominal convention: label is bin END. Also report CENTER and START conventions, all3 without outcome selection. Exact piecewise-constant hourly integrals are calculated over each proposed physical interval. Positive hourly temperature defines degree-day input; it does not assert exact subhourly positive-degree exposure. Full support is required. No gaps or missing values are replaced by zero.

SNOW is missing after4July. It is never silently set to zero. The models below do not explicitly model snow; snow-related structural errors remain in their evaluation. Do not add an event mask after seeing model errors. Do not use supplied conductivity or stake-specific test quantities for fitting. Measured weather during the evaluation interval is prescribed input, so this is a weather-conditioned ablation hindcast, not a future-weather forecast.

## Frozen models

All3 comparators are established one-driver empirical baselines with ONE nonnegative local coefficient fitted from the FIRST stake interval only:

1. Constant-rate persistence: melt per observed calibration day.
2. Positive-degree-day (PDD): sum of max(hourlyT,0) integrated inKdays.
3. Net-shortwave-dose: sum of max(SWIN−SWOUT,0) integrated inMJ/m².

PDD is the predeclared PRIMARY comparator of interest. The other2 are retained strong alternatives, not removed if better. Neither PDD nor radiation is claimed to be a novel melt model or complete energy balance. Conductive properties, heat storage and snow dynamics are not resolved. Two independent heat coefficients cannot be calibrated from the single local ablation interval; that model is deliberately not fit. Different calibration periods across stakes do not justify sharing unknown local coefficients when their debris properties differ.

For each weather-clock/visit-time combination refit coefficients ONLY using its calibration interval, then predict the second interval. Report all81 scenarios per stake, all4stakes and all3models:972 rows. This is a sensitivity set, NOT972 independent field experiments. No scenario is chosen by error. The nominal results and full primary timing ranges are both shown.

## Metrics and accuracy statement

Report primaryH1D signed/absolute error inmmw.e. and relative error; all4 secondary stake results and equal-stakeMAE; and every timing scenario. A predeclared descriptive robust-win flag requires PDD to have strictly smaller absolute error than BOTH alternatives atH1D in EVERY declared timing scenario. It is not a p-value or scientific-validation certificate. No arbitrary post-hoc accuracy threshold is introduced.

The source says uncertainty0.1m without specifying a standard deviation, interval/reading association or a probability distribution. Report TWO conditional bounded interpretations, not confidence intervals. With f=test/train dose ratio and δ=0.1m×900kg/m³=90mmw.e., per-reading independent bounds give residual uncertainty2(f+1)δ, because the middle reading is shared. If0.1m instead applies per ablation interval, the bound is(f+1)δ. Neither is asserted to be the authoritative statistical meaning. Compatibility under these assumptions is not proof of high accuracy.

## Freeze and reporting

Publish source/data hashes and this protocol BEFORE computing predictions. Preserve complete records, source coefficients, all timings and failures. If source changes after freeze, create a new version and do not overwrite results. If PDD loses, retain that verdict. Even a favorable result would establish only a small retrospective weather-conditioned, sub-debris ice-ablation evaluation at one glacier—not exceptional generalization, independent cliff validation, originality or guaranteed publication/competition outcomes.
