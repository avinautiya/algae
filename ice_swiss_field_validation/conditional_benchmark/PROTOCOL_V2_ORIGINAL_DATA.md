# Original-data retrospective comparison freeze

Supersedes PROTOCOL.md input limitations after recovery of the original archive. Preserve both protocols and their hashes. No thermal score has yet been computed at this freeze. Targets and the source paper have already been inspected: this is a retrospective temporal evaluation, not a blind validation or independent preregistration.

## Inputs and timestamp conventions

Use the CRC/size/SHA-verified original files under `/private/tmp/ogier-author-original-data-20261101`, pinned by SOURCE_MANIFEST.json. Original direct floor-elevation CSVs now supply the epoch errors, unlike the previously recovered comment arrays. Base interpretation uses their explicit UTC header. The paper uses CEST for narrative times and processed temperature files lack a zone: report manual-epoch shifts of −2, 0, +2 hours relative to temperature products as separate engineering alignment scenarios, not statistical bounds.

Temperature forcing uses native `temp3_corr_raw.csv` and `temp5_corr_raw.csv`, plus the final author export corrections **once** (+0.035 and +0.040°C respectively). Integrate by elapsed time. The daily products use sample arithmetic means despite variable logging frequency, so they are not exact heat-budget integrals. Do not substitute drawn daily curves or sample means as the primary forcing when native data exist. Retain missing stretches and reject unsupported intervals consistently across methods. The native integration module records its declared reconstruction and gap rule. Sensor shifts −0.05/0/+0.05°C are engineering scenarios; temperatures below melting are clamped to zero only with counts and the convention reported. The shift does not bound unknown floor-to-bulk bias.

Use author-processed `Reynold_daily.csv` as imposed calendar-day Reynolds numbers, constant within each day. It is a **derived hydraulic forcing** (including reconstructed discharge over a gap), not independently observed daily velocity. No author Nu or fitted heat-transfer products enter predictions. Diameter scenarios are 1.08/1.30/1.52 m, not fitted to the test targets. Re is retained at its authored value across these diameter scenarios: this varies the transfer length at fixed imposed nondimensional forcing rather than claiming a dynamically consistent geometry perturbation. Report that limitation explicitly.

Fixed water/ice properties: Pr=13.5, k=0.56 W/(m K), rho_i=917 kg/m³, latent heat=334000 J/kg. These are declared 0°C engineering constants, not fitted values. The test estimates local vertical thermal lowering at melting-point ice, without atmosphere, mechanical erosion, moving reference correction, pressure-melting correction, or changing cliff contact. No claim of full coupled-system field validation follows.

## Split, models and endpoints

Use all adjacent manual P3 and P5 bed observations. For each station the split epoch is its last manual visit on or before 25 July 2019. Fit on intervals ending no later than that epoch; evaluate every later interval beginning at or after it that has complete common forcing support. This includes negative/small observed bed changes; never delete them because a melt model cannot reproduce them. All models use the identical accepted intervals. Report unsupported intervals and reasons before scores.

All six models are fixed now:

1. One pooled nonnegative constant melt rate, fitted on training interval lowerings.
2. One pooled nonnegative constant heat-transfer coefficient h, fitted on training temperature integrals.
3. Unfitted Dittus–Boelter with Pr exponent0.4, as used in the published comparison.
4. Unfitted DB with fluid-cooling exponent0.3.
5. Unfitted smooth-friction Gnielinski.
6. One pooled nonnegative multiplicative coefficient on the DB0.4 heat-transfer law, fitted on the same training intervals.

Fits minimise unweighted squared errors in interval lowering metres using one scalar and closed-form nonnegative least squares. Candidate6 and baselines1/2 each have one fitted parameter. No station-specific fitting, test retuning, family selection or learned final temperature correction. Refit training-only scalars separately under each engineering forcing scenario; preserve every scenario, never choose its best test score.

Report nominal interval predictions, per-station MAE/RMSE/bias in metres, pooled scores, observed/predicted totals and every exclusion. Compare candidate6 with **both** fitted simple baselines, not merely with unfitted DB. Aggregate and component intervals share observations and cannot count as independent replications. Paired bounds use original epoch uncertainty conservatively, with clear distinction between standard deviations and deterministic engineering bounds. Do not label ±2 summed epoch errors as a confidence interval.

Success would mean improved conditional temporal prediction at both measured stations against both fair fitted baselines, surviving reported engineering scenarios. A failure at one station must remain visible. Even success is one Swiss stream event and cannot establish general environmental consistency or professional glacier-model superiority.

Sources: [Ogier et al. methods and equations](https://tc.copernicus.org/articles/15/5133/2021/); immutable author repository commit9116f98dc7fef42be82dc4de50316f4a676ca8e9; original archive DOI10.3929/ethz-b-000504956. Standard correlation implementation provenance is in PROTOCOL.md and standard_closures.py.
