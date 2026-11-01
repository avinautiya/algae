# Conditional Swiss channel thermal-component comparison

Written before thermal scores, after observation-target inspection. This is an analysis freeze, not independent preregistration or a blind discovery dataset. No thermal-score run has been performed in this package.

## Observation target and support

Use **direct manual channel-bottom elevation changes**, not daily water-surface cuts. Source: author-comment arrays in `Channel_geometry.ipynb` at commit `9116f98dc7fef42be82dc4de50316f4a676ca8e9`, independently extracted in `../../ogier_author_repo/channel_elevation_authored_comments.csv`. Original standard-deviation columns remain unavailable. Do not relabel authored comments as downloaded original measurement CSVs.

| Station | Start (authored local clock) | End | Bed lowering |
|---|---|---|---:|
| P5 | 16 July 2019 14:00 | 30 July 2019 14:31 | 5.50 m |
| P3 | 16 July 2019 13:17 | 30 July 2019 13:20 | 4.32 m |

Primary endpoints are these two aggregates. Secondary disjoint interval endpoints are the intermediate manual visit on 25 July at 10:28 for P5 (3.64 then 1.86 m), and 24 July 14:38 plus 25 July 09:48 for P3 (2.05, 0.73, 1.54 m). Do not count aggregates and their components as independent replications. The two stations share one stream and weather event.

Paper methods place measurements relative to drilled stakes/crossbeams, with GPS elevation; this is not an independently corrected moving-ice reference. Typical crossbeam vertical accuracy is 0.1 m; P3/P5 temperature accuracy is 0.05°C. The thermistors were near the bed, so recorded temperature need not equal bulk-water temperature. See [Ogier et al. methods](https://tc.copernicus.org/articles/15/5133/2021/#section3). Preserve these systematic issues rather than treating them as Gaussian errors. Original epoch uncertainties are missing; use ±0.2 m endpoint-sum sensitivity on changes only as a documented typical-error scenario, not a measured interval or confidence bound.

The P5 waterline imprint aggregate is 5.65 m, a **different** target, with an observed approximately zero endpoint-depth correction and a 31-minute endpoint discrepancy. Report as an auxiliary observation check only. Its difference from 5.50 m direct bed lowering is 0.15 m; do not pick whichever target yields a better score.

## Gate before a comparison can run

1. Recover measured corrected daily P3/P5 temperature vertices, including date/noon convention and source pixel uncertainty. Do not add the author's correction offsets again or invert incision to obtain temperature.
2. Derive daily averaging windows from helper source. A noon label alone does not establish the averaging window; inspect the actual grouping. Confirm common clock convention. Keep naive source dates until justified rather than assigning UTC silently.
3. Exact integration to manual visits must disclose constant-within-daily-bin assumptions and endpoint fractions. No extrapolation or unobserved gap filling.
4. Velocity records are instantaneous **segment** tracer velocities, not a daily local series. Freeze linear interpolation as an imposed hydraulic scenario, not an observation. Required support 15–31 July exists; no extrapolation. Width and hydraulic diameter are insufficiently observed for a uniquely measured local time series. Therefore use the same imposed diameter scenarios (1.08, 1.30, 1.52 m) for all methods; these are declared study geometry scenarios, not confidence limits. Document their paper provenance separately.
5. Floor-to-bulk temperature bias is unknown. The ±0.05°C sensor scenarios do not bound that bias. Any score remains conditional on treating measured near-bed temperatures as bulk forcing.

## Unfitted comparators and evaluation

`standard_closures.py` implements DB as used in the author's notebook (A=0.023, Re exponent 0.8, Pr exponent 0.4), the conventional cooling exponent 0.3, and smooth-friction Gnielinski. Both DB variants must be retained; neither may be selected after scoring. Formula and heat-to-lowering conversion follow [Ogier equations 2, 7, 9–10](https://tc.copernicus.org/articles/15/5133/2021/#section3). These are existing models, not claimed new methodology.

Smooth Gnielinski calculates its Darcy factor from Re; a supplied measured-friction variant is a separate optional scenario. Never set a convenient friction value or fit Nu using the same incision targets and call the result independent. Constants must be explicit inputs, with source provenance in the executable campaign spec. No fitted physical constants are permitted.

Output each method's predicted and observed lowering, signed and absolute errors in metres, relative error, and all geometry/temperature/hydraulic scenario results. Include both stations and every accepted interval. Do not average percent errors around near-zero targets or report nominal sensor coverage as calibrated uncertainty. Aggregate error and interval error answer different questions.

Numerical Re/Pr gates do **not** establish application validity: these correlations were developed for tube flow, whereas this dataset is a rough free-surface ice channel. The [INL MOOSE implementation](https://mooseframework.inl.gov/docs/doxygen/modules/ADWallHeatTransferCoefficientGnielinskiMaterial_8C_source.html) provides a reproducible smooth-friction expression; this package conservatively restricts Re to 3000–5e6 and Pr to 0.5–2000. [ORNL's implementation report](https://info.ornl.gov/sites/publications/files/Pub69409.pdf) distinguishes DB fluid-heating/cooling exponents. The original Gnielinski paper has not been read in full here.

This can test **conditional local thermal-closure adequacy** against a field event. It cannot establish new coupled hydraulics skill, cliff survival, independent environmental transfer, or consistent superiority over professional glacier models. A new proposed model must be specified before scoring and evaluated alongside all comparators on identical inputs and endpoints. A successful two-station comparison would remain one-event evidence, not proof of cross-environment consistency.
