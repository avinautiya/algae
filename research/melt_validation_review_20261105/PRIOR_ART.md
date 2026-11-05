# Focused prior-art review

Review date: 5 November 2026. This is a targeted review of the current ablation-index model, not an exhaustive review of all glacier science. Searches covered debris enhanced temperature-index models, debris thickness, empirical model transfer, robust regression, and model intercomparison. Primary sources below were checked in accessible full text, author records or source metadata; access depth is distinguished. Absence from these searches does not prove historical originality.

## Closest work and claim boundary

| Primary work | What is already established | Consequence for this manuscript | Access |
|---|---|---|---|
| [Carenzo et al. 2016, An enhanced temperature index model for debris-covered glaciers accounting for thickness effect](https://doi.org/10.1016/j.advwatres.2016.05.001) | Combines temperature, radiation and thickness-dependent empirical melt responses, with temporal transfer and an Arolla test. | Neither those inputs nor accounting for debris is new. Our positive interval-dose regression differs in calibration and aggregation, not in discovering debris physics. A faithful DETI implementation remains an important missing comparator. | Author repository record and accessible indexed methods; attempted full-text PDF retrieval was denied. No claim of a complete equation-level reproduction. |
| [Winter-Billington et al. 2020, Evaluating the transferability of empirical models of debris-covered glacier melt](https://doi.org/10.1017/jog.2020.57) | Evaluates pooled and hierarchical empirical models with glacier/year structure, transfer and uncertainty. | Cross-glacier empirical transfer and accounting for grouped observations are established. Equal-glacier weighting is not a substitute for their hierarchical uncertainty model. | Full publisher HTML: data preparation, model fitting and evaluation sections inspected. |
| [Pellicciotti et al. 2026, DCG-MIP](https://doi.org/10.5194/tc-20-1895-2026) | Compares multiple debris-melt formulations on a common observational basis and assesses debris-property sensitivity. | Our small benchmark must not be described as the first multi-glacier comparison. Published scores cannot be mixed with our different targets/support to claim superiority. Actual matched implementations are required. | Full publisher HTML; PDF retained locally. |
| [Reid and Brock 2010, An energy-balance model for debris-covered glaciers including heat conduction through the debris layer](https://doi.org/10.3189/002214310794457218) | Resolves surface-energy and debris-conduction processes that empirical doses compress. | Our model is not a full energy-balance solver; it cannot claim new thermal physics. A radiometric Fourier comparator is not this full model. | Bibliography and discussion in checked primary intercomparison; original full equations not re-audited in this pass. |
| [Nicholson and Benn 2006, Calculating ice melt beneath a debris layer using meteorological data](https://doi.org/10.3189/172756506781828584) | Uses meteorology and a debris heat-transfer formulation to estimate sub-debris melt. | Simple debris-mediated melt prediction is longstanding. | Publisher abstract/search record; full derivation not re-audited in this pass. |
| [Reid et al. 2012, Including debris cover effects in a distributed model of glacier ablation](https://doi.org/10.1029/2012JD017795) | Applies debris-aware ablation modelling spatially and reports field parameter constraints. | Local point predictions here do not demonstrate a new distributed glacier model. | Publisher methods excerpts and archive metadata; not a full software reproduction. |

The candidate minimizes a nonnegative, penalized absolute-error objective. Absolute-error regression and linear programming are established statistical tools. Replacing squared error with absolute error is an implementation/methodological choice, not a new theorem. A publishable incremental claim would require a reproducible, independently tested advantage or a consequential explanation of failure regimes.

## Exact current method

For interval i with debris thickness h_i, define degree-hour dose D_i = sum max(T,0) delta-t and radiation dose S_i = sum max(SWIN-SWOUT,0) delta-t. The fitted interval prediction is a positive combination of D_i/h_i and S_i/h_i. Feature scales are donor-only weighted RMS values and the target scale is the donor-only weighted mean. The loss is weighted absolute residual plus a nonnegative-coefficient L1 penalty. Optional rate normalization and penalty are selected within donor-glacier cross-validation. The outer recipient supplies neither fitting targets nor normalization statistics.

This form has no snow state, surface heat storage, turbulent energy flux, cold-ice state or explicit debris-moisture dynamics. It does not calculate calving, cliff retreat, flood risk, algae optics or future weather. Those earlier project directions must not be included as implemented achievements in this paper.

## Review decision

Historical novelty of the physical model is **not established**. The five-glacier development gain is real under its declared operator but is incremental and sensitive to geometry. The new frozen Djankuat surface-lowering test fails against simple controls. At present, a claim of a newly superior transferable melt model is unsupported, and Nature-level scientific significance has not been demonstrated. A transparent calibration/transfer audit is the defensible manuscript direction; its journal fit requires mentor assessment.

Next literature task before a model-superiority submission: obtain and reproduce DETI and the published mixed-effects model under the same observation support, units and no-recipient-calibration rules. Do not substitute an informal approximation while retaining a published model's name.


## Follow-up retrieval and comparisons

The follow-up retrieved the original Carenzo2016full-text XML through Europe PMC and visually checked the Winter-Billington2020and DCG-MIP2026parameter tables. This supersedes the earlier access-depth limitation. A recipient-excluded new2022-2023Piramide direct-stake comparison now exists. KP1fixed-effects equations are implemented; four DETI equation reconstructions retain explicitly unresolved threshold/lag execution conventions. The candidate loses to all four reconstructions. See followup/FOLLOWUP_REPORT.md; original historical comparisons are unchanged.
