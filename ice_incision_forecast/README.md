# Sharp forecast bounds for glacier channel incision

This contribution contains a proof manuscript, an executable conditional forecasting tool, independent numerical checks and six reproducible method figures. It is separate from the active CFD milestone and changes no existing solver or case.

## Proved results

`THEOREMS.md` establishes, under an explicit reduced thermal model:

1. An exact nonlinear map including meltwater mass addition, and an observational-equivalence result for outlet records.
2. Sharp regional heat/incision bounds with complete attaining spatial profiles.
3. A common pair of spatial profiles that attains the finite-horizon bounds across all permitted inlet histories and common conductance multipliers.
4. The optimal worst-case scalar forecast and sharp timing bounds for a declared incision threshold.

The general scalar midpoint-minimax principle and conservation laws are established mathematics. The project-specific sharp allocation formulas and forcing-uniform construction are documented with proofs and checks. Their mathematical priority has not been certified; `PRIOR_ART.md` records what was actually examined. These results do not establish a new free-surface instability or environmental validation.

## Run the forecast

The forecasting core and tests use Python 3.9+ standard library only, with no CFD installation or project-package dependencies.

```sh
python3 -m unittest discover -s ice_incision_forecast/tests -v
python3 -m ice_incision_forecast demo --output /tmp/ice_case_01.json
python3 -m ice_incision_forecast predict --case /tmp/ice_case_01.json --output /tmp/ice_forecast_01
python3 -m ice_incision_forecast.validation --output /tmp/ice_vertex_check_01.json --cases 500
```

Choose new output paths. Existing results are preserved. Outputs include cellwise and regional lowering ranges, a feasible nominal bed profile, scalar regional minimax forecasts, threshold timing, outlet temperature/discharge, mass/energy ledgers, explicit witnesses and hashed provenance. Future outlet temperatures or future bed measurements are not accepted inputs.

This is a forward **scenario forecast** driven by prescribed inlet conditions and a conductance prior, not a forecast of meteorological forcing. It can be run with an interval on total conductance; an exact endpoint constraint is not required. The default demonstration uses an exact total to expose the remaining spatial ambiguity.

## Figures

Figure generation requires NumPy and Matplotlib; these are optional for the forecasting core. Versions used are recorded by the figure environment/source manifest.

```sh
python3 -m ice_incision_forecast.figures --results /tmp/ice_forecast_01 --output /tmp/ice_figures_01
python3 -m ice_incision_forecast.plot_published_cfd --output /tmp/ice_cfd_figure_01
```

The six-figure script expects the supplied exact-total, no-prefix demonstration with an interior target region. General case forecasting does not have this plotting restriction.

1. Equal outlet predictions with different regional incision.
2. Sharp marginal spatial envelope and feasible-profile illustrations.
3. Forward regional forecast and declared incision threshold.
4. Ideal interior thermal information and remaining local ambiguity.
5. Independent ODE/extreme-value checks and meltwater-addition effects.
6. Local versus reach-scale information requirements.

The second command produces a seventh figure by replotting already committed OpenFOAM benchmark metrics from the pinned `9a52524` source. It does not rerun the solver. A numerical reporting discrepancy in the earlier steady-sensitivity wording is documented with the actual stored differences.

Every figure is saved as PNG and vector PDF with its source data, captions and provenance. The first six are synthetic conditional-model demonstrations, not field observations or resolved-CFD validation. The seventh uses prior repository verification metrics, not new CFD runs. A marginal upper-envelope curve is not a single jointly attainable spatial forecast.

## Input contract

The demo JSON is the complete schema. Each cell specifies length, width, initial bed/depth, and lower/upper/nominal conductance in W/K. The prior specifies a total-conductance band and optional **exact** interior prefix conductances at cell boundaries. Real uncertain point-temperature data must not be silently converted into exact bulk-temperature/prefix constraints.

Forcing intervals specify positive duration and inlet Q, nonnegative temperature excess up to 4 K, and a nonnegative common conductance multiplier. Properties are constant and explicit. Named contiguous regions can include a positive `critical_lowering_m`; it is an incision threshold, not an automatically inferred hydraulic-failure criterion.

The model assumes a fixed temperate ice bed, no extra heat sources or thermal storage, prescribed inlet discharge and a fixed conductance shape. All meltwater joins the flow at the melting point. It does not resolve turbulent/free-surface motion, pressure controls, lateral wall erosion, cold ice or evolving hydraulic geometry. The code flags user-declared geometry-change and transit-time screening limits; those thresholds are not universal applicability theorems.

## Scientific use

The theorem makes a specific distinction: matching downstream thermal observations can leave local incision unidentifiable in this model. The tool computes the best possible conditional range and scalar worst-case forecast under the stated information. That is useful for evaluating what a model/data combination can support, even before asserting a new dynamical mechanism.

To improve a real glacier forecast, the conductance prior must be physically justified and the reduced-model assumptions checked against the verified CFD/observations. The final field still needs independent accuracy and prior-art assessment. See `HANDOFF.md` for later integration without interrupting the active solver work.
