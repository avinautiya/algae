# Labelled profile artifact after the external parameter freeze

Nine relative profile values were recovered at the explicitly labelled sections A, B and C for 25, 26 and 29 August 2005. They are samples of the drawn measured longitudinal profiles, not original raw survey files or fitted model curves. The native embedded 709×367 JPEG was extracted from publisher PDF page 94 without resampling. Its SHA256 is enforced by the executable extractor. The selection proof was visually inspected; black circles lie on the red, purple and cyan profiles at the selected section lines.

The published raster clips the leading altitude digits: vertical labels appear as `40,20,00,80,60`. To avoid inventing absolute altitudes, the extracted elevations use the central `00` tick as a zero reference. The spacing establishes metres, and differences between survey dates do not require the missing additive altitude offset. No absolute elevation was assigned.

Source pixel coordinates, nearby colored-line extent and conservative graphical halfwidths are retained per row. These bounds are approximately 0.82–0.96 m, reflecting the coarse original raster. They are **not** survey-instrument errors or confidence intervals. Bed-change bounds are conservative sums of endpoint graphical bounds; shared errors make differences dependent. No quantitative high-accuracy claim follows from this resolution.

The horizontal tick calibration places the plotted A/B/C section lines roughly 3/17/33 m downstream of the plotted lake-breach origin. These graph-derived positions have about ±0.53 m selection bounds and do not prove that section A coincides with a critical hydraulic control. In particular, the lake temperature cannot automatically be used as local temperature at all sections.

`FIG5_DATE_ONLY_BED_CHANGES.csv` supplies geometric lowering for three overlapping date pairs at each section. Survey hours remain unknown. Consequently it retains calendar-date duration envelopes rather than inventing noon timestamps or precise rates. The adjacent-date pair permits a very broad elapsed-time envelope; it must not be used as an exact 24-hour rate. These nine differences share endpoints and satisfy algebraic dependencies; they are not nine independent glacier events.

The external parameter freeze was reported before target extraction. No prediction was scored, no model coefficient was changed and no observation operator was selected by agreement with these values. The unlabelled Figure 6b series is still not assigned to A/B/C. Exact field scoring remains gated on section binding, time support and a matched point-bed operator.

Reproduce:

```sh
python extract_labelled_profiles.py page4_original-000.jpg --out new_extraction
```

The source DOI and PDF hash are in `OBSERVATION_MATCH_GATE.json`; extraction anchors and the parameter-freeze hash are in `FIG5_EXTRACTION_MANIFEST.json`. This bounded artifact makes the usable observations explicit without converting unresolved input/target matching into a claim of field superiority.
