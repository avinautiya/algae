# Verified melt targets, not another heat-flux proxy

Downloaded and checksum-verified five public DCG-MIP site archives from their original Zenodo records. Arolla is the strongest immediate component-test candidate. No prediction or accuracy score was computed by this audit.

## Primary: Arolla 2010 stake-inferred ice ablation

Source DOI10.5281/zenodo.3047649; archive `1_ARO_2010.zip`. Its original workbook contains raw stake readings in centimetres, cumulative surface-height differences, and a separate **MELT [mm w.e.]** block. Source density is900kg/m³. The loader reads the provided water-equivalent targets directly and verifies their identity with the source schema, rather than converting a heat-flux sensor fit into alleged melt.

The primary stake H1D is explicitly the closest to the AWS. A measured debris thickness0.06m has a stated generic10% uncertainty and is linked to that stake/nearby AWS. Meteorological data cover22May–12September2010, fully enclosing every possible hour of the stake dates28July,28August and11September. The first interval provides846mm water equivalent; the second provides315mm. These are a train/evaluation design opportunity, not evidence that any proposed predictor is accurate. The two intervals share the same stake and glacier; only one held-out weather interval exists for the primary design.

Three additional stakes provide two intervals each. They are farther from the AWS, and their site-specific debris thicknesses are not supplied. Their local training calibrations may support secondary comparisons but must not be relabelled as independent glacier replications or independently measured material parameters. G1D's intervals are666/315mmw.e.; F1D's909/315; E1D's630/90.

## Timing, snow and uncertainty must remain explicit

The source declares local timeUTC+02:00. Stake visits have day precision only; joint visit-hour ranges0–24h can be propagated because the forcing covers the complete envelope. Do not invent visit times or call a chosen noon scenario a native timestamp. Hourly means' precise label interpretation remains a one-hour sensitivity, not permission to adjust phase by test fit.

Source metadata state only `Uncertainty =0.1m` for stakes. It does not specify Gaussian sigma, per-reading versus interval interpretation or a confidence level. One may propagate alternative explicitly labelled engineering bounds; one cannot invent95% intervals or claim that fitting inside such a range proves statistical calibration.

The native SNOW variable is unknown/NaN after4July. Missing snow must NOT be changed silently to zero. Source metadata describe transient events on24July,30July,6August and30August (melting within24h), plus15August snow persisting3–4days. A predeclared metadata-calendar mask is a conditional closure, not a measured snow record. The stake target is ice ablation; predicted snowmelt must not be added to it.

Conductivity0.96W/m/K is supplied but its DOCX table row has no provenance note. The explicit Miage borrowing note belongs to the roughness row; it must not be silently reassigned to conductivity by flattened text extraction. Neither independent local measurement nor a particular borrowed source for conductivity is established here. The candidate test should avoid using it or resolve its provenance separately. Excluding target-fitted material inputs is essential.

## Other inspected sites

* **Miage3050557:** eight real ablation stakes, three visits and observed weather; no snowfall reported. However, the supplied per-stake conductivity was explicitly derived from the SAME ablation measurements and must not enter an independent test. Forcing begins13June; stake101 begins12June. Forcing ends8September00h; seven final visits are8September with unknown hours. Stake114 has all visit-day envelopes inside the forcing (16June/28July/7September). Other cases need explicit boundary reconstruction or restricted support, not an unnoticed missing day.
* **Pirámide3056072:** one stake200m from AWS, literature conductivity rather than same-target inversion. Measurements4December2014/10February2015/8April2015 fall inside forcing; increments1.08m and0.50m. Local-clock convention and site-offset assumptions remain. Initial15April measurement is outside the forcing and is unnecessary for those two interior intervals.
* **Suldenferner3056524:** UTC and independently derived later-year conductivity are useful, but the public target is SR50-derived displacement, not a separate stake series. Snow and height-derived forcing need the same issue-time and target-meaning audit as Changri Nup.
* **Tapado3362402:** public processed ICE ABL and snow/quality flags exist, with sensor reinstallation and smoothed SR50 processing. This is a conditional processed-displacement target, not the preferred independent stake target for the current request.

## Claim ceiling

Arolla permits an actual observational test of **sub-debris point ice-ablation prediction** under declared visit-time and snow assumptions. It does not validate the stream/cliff model, vertical-face retreat, a whole glacier, Himalayan transfer, historical originality or superiority over professional models. Prediction/calibration/test rules must be frozen before computing metrics. Do not select a winning model or site after seeing its error.

## Execute

```sh
python3 -m unittest -v test_arolla_loader
python3 arolla_loader.py 1_ARO_2010.zip AROLLA_RECORD.json NEW_PARSED.json
```

Five tests pass, including real archive targets/density and complete visit-day support. Output is refused if it already exists. Original source members, archive checksums and target precision remain in the parsed record.
