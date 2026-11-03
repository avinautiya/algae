# Three-sensor weak-energy algorithm and observable target

## Implementation status

`three_sensor.py` and `predict_thermal.py` provide an executable three-sensor thermal reconstruction with integrated consistent-P1 energy balances, a two-variable passive NNLS fit, rank/excitation refusal, signed basal-gap storage accounting, exact source-edge aggregation, an explicit output-scale quotient, and a past-only parameter-update API. Eleven independent synthetic tests pass. No Lirung ablation magnitude or uninspected outcome was read. No field performance result is asserted.

The mathematics is established finite-element heat balance and constrained least squares. The implementable candidate contribution is a source-audited short-campaign thermal-to-basal-energy pipeline with explicit observability/phase abstention. Its historical originality and predictive improvement remain unestablished.

## Source geometry and campaign support

The prior source audit is `round4/external_lirung/ELIGIBILITY_REPORT.md`, based on [Chand and Kayastha2018](https://www.climatenepal.org.np/sites/default/files/doc_resources/Chand%20%26%20Kayastha%20et%20al.%2C%202018%20SEAR.pdf) Tables1–2 and [Zenodo4446387](https://doi.org/10.5281/zenodo.4446387). Complete monsoon site2 has sensors5/15/35cm and total debris42cm, leaving an unmeasured7cm basal gap. Site1's duplicate20cm CSV header conflicts with paper30cm; this implementation does not repair it. Some winter blocks change sensor depths or reverse clock order. Exact native block mapping is a separate source gate; matching elevations is not proof of co-location with an ablation stake.

The longest block is11days, so the CLI accepts actual positive source durations up to11days and never invents14-day training/evaluation endpoints. It requires the exact source-native timestamp labels, source SHA256 and three ordered depths. It refuses nonmonotone clocks, missing values, unresolved excessive gaps and unbound aggregation endpoints. The source's rock2700kg/m³, specific heat750J/kg/K and porosity0.33 imply a dry matrix contribution1.35675MJ/m³/K; these are assumptions, not independent field measurements. No heat-capacity value is silently estimated from a test melt target.

## Integrated inverse equations

Let depths z1<z2<z3, h1=z2−z1, h2=z3−z2. Assume piecewise-constant diffusivities alpha1,alpha2 and a common constant volumetric heat capacity C. For the middle P1 node, divide the heat equation by C and integrate over an observed interval[a,b]:

    h1/6 DeltaT1 +(h1+h2)/3 DeltaT2 +h2/6 DeltaT3
      = alpha1/h1 integral(T1−T2)dt
        +alpha2/h2 integral(T3−T2)dt.

Every term has units K*m. The consistent mass terms include neighboring endpoint changes. Temperature integrals use declared piecewise-linear reconstruction. Windows are assembled at actual observed edges; partial terminal windows are retained with duration normalization. Default6h balance windows are a declared engineering configuration, not an outcome-selected optimum.

Each window supplies one row for two unknown diffusivities. NNLS searches the unconstrained interior candidate, both nonnegative faces and the origin. This is solved directly without a numerical-optimizer dependency. Multiple dynamic windows can make the two columns independent. A steady profile generally constrains only the diffusivity ratio and fails the rank gate, even when a numerical least-squares solution exists. A flat profile carries no diffusivity information. The configurable default minimum singular-value ratio0.001 is a declared practical safeguard, not a universal identifiability threshold or confidence level. Values on zero-conductance faces fail the positive-material point-estimate gate.

With an externally justified per-sensor error bound delta, a necessary residual outer bound is

    |r| <= delta*(h1+h2)/duration +2delta*(alpha1/h1+alpha2/h2).

The code reports consistency counts only. Correlated errors in the response and regressors remain; this is not a posterior covariance or proof that one coherent sensor-error history fits all windows.

## Basal energy target and unmeasured gap

For a declared ice depth zi>=z3, gap g=zi−z3, and a declared constant basal temperature Tb, approximate the temperature field linearly between measured sensors and from sensor3 to the ice. The sensible-temperature moment below the midpoint of the deepest measured edge is

    H/C = h2*T2/8 +3h2*T3/8 +g*(T3+Tb)/2.

The normalized integrated heat reaching the ice is

    Qice/C = alpha2/h2 integral(T2−T3)dt −Delta(H/C).

This is top-of-control-volume conductive heat minus its stored sensible-energy change. It does not assume the deepest thermistor is at0C or at the ice. It also does not independently determine the actual basal gap, cold ice, wet-pore latent energy, advection, lateral heat flux or evolving geometry. Unknown geometry/phase is reported rather than invented. The incoming deepest-edge flux is an effective P1 gradient approximation, not a resolved microscopic basal-flux measurement.

Absolute heat needs externally specified C. Simultaneous k→s*k,C→s*C leaves observed temperatures and diffusivity unchanged while scaling Qice by s. The temperature record alone therefore cannot identify absolute heat capacity, conductivity or ablation scale.

## Scale quotient and phase gates

If an independently bound, positive TRAIN melt target exists, a calibrated empirical prediction is

    Yhat_eval = Y_train * (Qice/C)_eval / (Qice/C)_train.

The common absolute capacity scale cancels exactly. This is not identification of physical C. An aggregate TRAIN target must not also be used as the evaluation target. If the archive supplies only one campaign-mean ablation measurement for the entire11-day thermal block, there is no within-block independent calibration-and-validation split unless independent subperiod outcomes are actually provided.

The absolute physical conversion Qice/L is exposed only behind explicit independent temperate-interface and complete-interface-budget declarations. Negative heat steps cause refusal: they require an actual cold-ice/refreezing state, not max(Q,0) clipping. A positive full-period total does not prove a valid melt-only trajectory. The normalized signed thermal retrieval can remain useful when melt conversion is unsupported, but it must retain that label.

## Causality

`reconstruct` fits diffusivities using all temperature covariates of the supplied block. It is a conditional thermal reconstruction/hindcast, not a future-weather forecast. It does not use camera targets. `fit_past` uses only temperatures with source labels<=asof and a declared history window; even future NaN temperatures do not affect it. A72h past-window configuration can be tested, but its warm-up leaves the first3days without a causal fitted estimate. Do not fill those days with hindsight parameters while claiming a full-campaign causal prediction. An external fixed-material prior could initialize them only as a separately declared model.

## Information-fair comparisons

1. An externally fixed conductivity times the same signed measured gradient is the conventional fixed-material flux comparator. Its C/k provenance must be independent of evaluation melt.
2. The same deepest-edge gradient with one local TRAIN melt gain tests whether inference/storage adds information beyond the gradient. Retain signed energy and identical clocks/footprints; any positive clipping is a separate model.
3. A constant-rate comparator uses the same one TRAIN target and exact duration ratio. Do not calibrate its rate using the evaluated seasonal total.
4. If native air/weather is absent, do not relabel buried temperature as air PDD or claim comparison to a faithfully reproduced full professional energy-balance model.
5. All calibrated methods receive the same local TRAIN target and exact native evaluation operator. The candidate adds two temperature-derived diffusivity parameters; that extra complexity and its rank/uncertainty diagnostics must be reported.

Full campaign-integrated thermal retrieval is executable now. Paired native melt validation remains conditional on the independent source-binding and phase gates. No source table or frozen previous result was changed.

## Reproduce

`python3 -B -m unittest discover -s algorithm -p test_three_sensor.py -v`

The CLI takes a JSON block with `source_sha256`, `source_block_id`, `depths_m`, `source_local_labels` and `temperature_C` (one3-value row per exact label):

`python3 -B algorithm/predict_thermal.py --input BLOCK.json --output NEW_RESULT.json --ice-depth-m 0.42 --basal-temperature-C 0 --balance-window-hours 6`

Output files are exclusively created, never overwritten. The generic block must be produced by the separately audited native-source decoder; this package does not guess the contradictory depth header.
