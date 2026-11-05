# Observation operators and uncertainty

This document separates verified archive facts, declared engineering assumptions and unresolved quantities. It does not assign probability distributions where sources provide only generic errors or no error estimates. Source metadata were re-extracted from the original archive DOCX files into this directory. Original values and source hashes remain preserved.

## Training and development observations

| Cohort | Target and geometry | Verified limits |
|---|---|---|
| Arolla 2010, H1D | Stake reading differences 0.94 and 0.35 m ice; measured debris approximately 0.06 m near AWS/H1D. | Metadata gives stake uncertainty 0.1 m without a fully specified error distribution, and generic thickness ±10% plus spatial variability. Metadata documents snowfall, including about 15 cm on 15 August with residual snow for several days. These intervals must not be labelled continuously snow-free. Adjacent intervals share a reading endpoint. |
| Miage 2005 | Native stake lowering and per-interval measured debris thickness. | Thickness ±0.02 m. Stakes are within 700 m horizontally and less than 80 m vertically of AWS. Several intervals share weather and visit endpoints. Conductivity estimates use these ablation observations; they cannot serve as independently measured validation features. A stake-error distribution was not identified in the extracted metadata. |
| Pirámide 2014–15 | Native PI3 stake reading differences 1.08 and 0.50 m ice, approximately 0.18 m debris. | Thickness ±0.02 m; stake approximately 200 m horizontally from AWS. One earlier reading lacks matching forcing and is excluded. No supported stake-error probability distribution was found. |
| Tasman 2017–18 | One manual stake interval, 2.0 m ice, approximately 0.30 m debris. | Thickness ±0.05 m, measured October 2017 and January 2018. AWS/stake footprints differ. Snow data are missing. The stake conversion ratio 0.917 differs from the draw-wire product's 0.9. No unmeasured later-season thickness is inferred. |
| Kennicott 2023 | Three native net-loss observations 1035, 873 and 688 kg/m² divided by the declared 900 kg/m³ ice density to give 1.15, 0.97 and 0.7644 m ice. | Weather is shared. P12 geometry uses a nearby profile rather than perfectly co-located stake geometry. Net mass loss is not independently measured gross latent melting. Density and component-process assumptions remain conditional. Source2020 forcing has118 unsupported hours and is refused. |

The common development cohort has15 intervals from five glacier campaigns. It does not have15 independent glacier replications. Within each site, weight stakes equally and then intervals within stake; weight glaciers equally. Report interval duration as well as total-loss error. Weather is observed, so scores evaluate weather-conditioned ablation hindcasts, not skill of future meteorological forecasts.

## Distinct processed products

CN ABL is positive lowering; Sulden ABL is negative lowering. The corrected exporter preserves those conventions explicitly. Earlier daily results with the wrong CN sign are superseded. Tapado supplies a filtered author-processed water-equivalent product, with a24-hour filter and source flags; its debris depth is a lower bound in metadata. These products must not be pooled with direct ice stakes without explicit observation mappings. Sonic lowering may include debris settling, snow and other surface changes. Source classifications derived from the same sensor are not independent truth.

## Tasman2018 comparison

The earlier frozen test targets minus the signed native Melt_mmwe sum,1203.3mm water equivalent. Instrument identities verify Melt_mm equals the Length_mm increment and Melt_mmwe=0.9Melt_mm. The PDD control's seasonal error1.05% coexists with maximum cumulative error156.8mm. The draw-wire anchoring, snow, reset/error procedures and 2018 geometry are not fully documented. Do not call this direct gross-melt validation or1% accuracy at every time scale. It was opened before the latest robust-model revision and is no longer an untouched test for that revision.

## New Djankuat archive audit

Official record: https://doi.org/10.5281/zenodo.3049871 . The API returned open files on this pass; the archive's MD5 and size were verified. This updates the older local 'restricted' audit without using any alternate access route.

The DOCX specifies sonic-ranger calculated lowering and provides no stakes/dGPS/dDSM validation. It gives h=0.65±0.02m, while the debris CSV gives0.61m. This unresolved discrepancy is preserved, not tuned to outcomes. Both depths lie beyond the maximum training thickness of approximately0.30m. Snow is missing; some meteorological observations were replaced from a nearby debris-free AWS. Air-temperature manufacturer error is±0.2°C at0°C; SWIN error is15W/m². Those specifications do not establish total field forcing uncertainty.

HEIGHT accounting verifies deltaHEIGHT=ABL+ACC in native row order, establishing negative ABL as lowering. The sign is not selected by predictive score. All complete target days with finite matched temperature/SWIN/SWOUT are retained, including zero lowering. No target-magnitude exclusions or clipping are applied. The weather and target HOUR24 placement is inconsistent in row order;0and+1hour complete-window mappings are both reported, not selected by score. This treatment does not fully resolve the physical clock convention.

There are81 eligible2007days and56 eligible2008days under the zero-shift convention. The shifted2007support has80days. Separate clock tables therefore have slightly different support and must not be interpreted as pure paired clock sensitivities. No native snow-free claim is possible. These are correlated days at one glacier in two seasons, not137 independent glacier trials.

## Uncertainty reporting policy

1. **Measurement uncertainty:** report provider statements and unsupported quantities separately. A generic±value is not automatically a Gaussian standard deviation.
2. **Input sensitivity:** recipient-only±20% thickness changes robust development MAE from0.282m to0.343/0.328m. Label these engineering scenarios, not confidence intervals. At Djankuat, h0.61/0.65gives robust daily MAE18.704/18.958mm in2007and17.263/17.526mm in2008; neither resolves failure.
3. **Sampling uncertainty:** five glacier campaigns and shared endpoints prevent treating intervals as independent replicates. No universal-performance confidence claim is supported.
4. **Selection uncertainty:** many model families were tried on inspected development outcomes. The reported9.7%gain is model-development evidence, not a preselected confirmatory effect or multiplicity-adjusted significance result.
5. **Structural uncertainty:** snow, changing debris properties, turbulent energy exchange and gross/net melt separation remain unmodelled. Observation screening cannot make those processes disappear.

Missing measurement definitions cannot be 'fixed' by assigning invented precision. A true independent ice-melt validation still needs native stakes or another justified ice-melt operator, measured applicable geometry, co-supported forcing and independent error information.
