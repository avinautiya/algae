# Independent discovery decision — do not substitute a promise for evidence

## Verdict

**No-go for claiming that the present model already has a distinct, consequential environmental discovery or field superiority.** The numerical repairs are useful. The fixed-forcing ideal recovery exclusion is a defensible mathematical result *within its stated model*, but the kinematics and finite-heat laws are established. Its agreement with an algebraically identical comparator cannot establish an improved predictor. No review can certify Nature significance or STS placement.

**Conditional go:** reproduce a professional-model benchmark on the exact published observation support, then test a single independently constrained source of prediction failure. A meaningful new result must survive the best appropriate comparator, measurement error and mesh uncertainty simultaneously.

## Compact decision matrix

| Candidate | Existing distinction / closest work | Present quantitative evidence | Current decision / what would change it |
|---|---|---|---|
| Fixed-Q debris lock-out and pulse rescue | Petersen2024 and Kneib2022 already associate stream disconnection with cliff abandonment; earlier audited Strickland2025 models evolving topography/debris movement. | The declared common-h/common-T ideal model cannot repair61mm in21days. Some small partial repair exists, but no sustained recovery. | Stop the rescue headline. New physical inputs must be independently justified; do not increase h/T to obtain a win. |
| Route-specific upstream cooling | Isenko2005 and other stream-cooling literature already model downstream heat depletion. SaDS already couples melt sources/routing. | Executed local ranking reversal has only0.0438mm/day drift in the supplied example. | No impact/priority headline. Require independently measured route inputs and a substantial target effect beyond uncertainty; compare to a standard cooling model with the same information. |
| More detailed heterogeneous debris conduction | Established transient conduction and energy-balance models, including professional DCG-MIP participants. | Earlier reserved temperature test improves on homogeneous diffusion but loses to ridge at every profile; no demonstrated ablation improvement. | No claim that complexity alone produces a better melt model. Physical parameters need independent constraints and a distinct predicted quantity. |
| Historical DCG-MIP reproduction | Direct comparison against submitted professional outputs, not an invented fractional-cover proxy. | Actual small output members are retrieved. Their exact units, interval mapping and snow/observation support are still gated. | Go only for reproduction after protocol recovery; this is not a new discovery or independent transfer test. |
| Physics-constrained uncertainty/measurement-aware prediction | Inverse problems, partial identification, uncertainty quantification and experimental design are established fields; priority in this glacier application has not been established. | The exact equivalence fixture below proves why the current sparse calibration cannot uniquely determine forecasts after a forcing change. | Conditional methodology candidate, NOT a new physical law. Must outperform equally informed uncertainty baselines at equal forecast coverage and show consequential observational/decision benefit. |

## Primary evidence and what it does not establish

The2026 DCG-MIP paper directly documents spatial mismatch between forcing and validation footprints, uncertain debris properties, and calibration/metric dependence. It explains that Changri Nup's supplied properties had been optimized to observed melt, and that cumulative and daily metrics can rank calibrations differently. Those are disclosed limitations of an intercomparison, not evidence that its authors committed an error. Source: https://doi.org/10.5194/tc-20-1895-2026 , sections6.1–6.4.

Other primary neighbors retained from the existing source audit:

* Petersen etal2024: https://esurf.copernicus.org/articles/12/727/2024/ — stream/cliff maintenance and abandonment are known.
* Kneib etal2022: https://tc.copernicus.org/articles/16/4701/2022/ — sub-seasonal cliff processes, including stream disconnection.
* Hill & Dow2021 SaDS: https://doi.org/10.1029/2021JF006309 — coupled supraglacial drainage/melt already exists.
* Strickland, Covington & Gulley2025 preprint: https://doi.org/10.22541/essoar.175647655.50202412/v1 — a closer evolving-topography/debris comparator identified in the pinned audit. Preprint status and any subsequent versions require explicit checking before a priority claim.

This is a bounded audit. It does not certify that every2026 paper, supplement or code release has been excluded as prior art.

## A different, coherent question if the mechanism path remains blocked

**When does apparently accurate melt calibration fail to identify a forecast under changed forcing, and can independent thermal observations restore a useful prediction bound?**

This produces melt-power/volume forecasts, uncertainty bounds and a precise minimum-observation requirement, rather than a list of correlating features. It does not promise superior average point prediction. If the original objective demands that result, this alternative needs a separate decision rather than quietly redefining success.

For finite water heat, `P=F T [1-exp(-h A/F)]`, with known F and A, one observed P defines a continuum of pairs `(h,T)`. `forecast_ambiguity.py` constructs two such pairs that fit training power exactly but predict different power at a different F. The inverse relation for fixed known T is `h=-(F/A)log(1-P/(F T))` when `0<P<F T`. In the film limit the observation predominantly identifies the product hT, not the separate factors. Additional network weights or fitted coefficients do not supply missing information. This is elementary established inverse-model reasoning, not an original-theorem claim.

The code's synthetic parameters are not estimates of Swiss field uncertainty. A real forecast interval must use independently justified bounds, correct observation operators and omitted-physics allowances. An assumed bound cannot become an environmental guarantee.

## Single fair falsification for a candidate prediction framework

**Question:** does an explicitly constrained identifiable-parameter forecast improve *out-of-calibration* melt predictions and uncertainty on the same observations, compared to the strongest equally informed physical and statistical baselines?

1. **Data gate:** reproduce one published baseline score first. Resolve units, interval boundaries, snow versus ice target, sensor movement, physical support and all target-derived inputs. Label supplied target-optimized parameters as calibrated. Choose a genuinely unused site/time partition after the metadata audit and disclose that audit.
2. **Numerical gate:** use a fixed physical domain; demonstrate independent space, time and source/geometry convergence. A short bare-floor test cannot validate covered-cliff lifetime. Predicted method differences must exceed combined numerical and measurement uncertainty.
3. **Identifiability gate:** derive the observable parameter combinations on training data only. An extra parameter needs an independent observation or a defensible prior. Test observable-equivalent parameter sets under withheld forcing, not merely optimizer covariance at one fitted point.
4. **Comparator gate:** include the best applicable submitted physical model, an unfitted physically justified correlation, a training-only statistical baseline, and the simpler model. Match support, available forcing and calibration budget. Do not compare an independently transferred candidate to a target-calibrated model without labeling the information difference.
5. **Freeze:** prerecord one primary melt metric on100% eligible support, a proper probabilistic/interval score, and secondary time-resolved diagnostics. Use a training/validation stage for every hyperparameter and uncertainty bound. Keep final targets out of model/threshold selection. Cluster evaluation by site/event, not individual sensor timesteps.
6. **Abstention rule:** if the new model withholds predictions, every baseline gets the same allowed retained fractions, e.g100%,80%,60%. Report error and coverage jointly. Never advertise lower error after simply excluding difficult cases.
7. **Kill:** if its forecasts equal an existing model with the same inputs, improvements disappear on matched support, gains are smaller than uncertainty, or the strongest baseline still wins, the claimed advantage fails. Keep that result. Do not rename the same failed method or retune on the test.

## What would count as a contribution

A transferable, quantitatively validated reduction of prediction error or false confidence under a clearly identified physical/observational regime would be useful. A precisely demonstrated limitation plus a minimal independent measurement that resolves it could also be useful, but is a different claim. Neither is historically original, high-impact or award-worthy by declaration. The work must first pass the stated test and a closer methods-literature check.

## Delivered computation

`forecast_ambiguity.py` and four passing tests check exact training equivalence, changed-forcing forecast divergence, finite-supply/film limits and inversion when independent temperature is supplied. No Claude code, frozen target or production job was changed.
