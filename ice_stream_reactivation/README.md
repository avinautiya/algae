# Glacier-stream reactivation research package

Start with [WORKING_RESULT.md](WORKING_RESULT.md). The leading candidate is a thermally evolving reactivation boundary after finite debris disturbance, including the case where all debris is exported yet ordinary flow fails to return.

This independent package contains mathematics, code, exploratory simulations, an independent numerical reference and figures. It does not modify Claude's project. All forecasts are conditional on a reduced two-route model; historical priority and resolved environmental transfer remain research questions.

## Reproduce

Python 3.9+, NumPy, SciPy, Matplotlib. From `compute/`:

```sh
python3 -m unittest -v test_model test_certificate test_fixed_stock test_topology
python3 spatial_check.py
python3 spatial_reference.py
python3 spatial_environment.py
python3 spatial_figures.py
python3 fixed_stock_forecast.py
python3 topology.py
```

The 25 component tests check hydraulic limits, stock and energy ledgers, melt accommodation, critical-depth force balance, patchy heat transfer and conditional exclusion tests. Separate `test_alternative_budget.py` at the package root tests a sediment-budget necessary condition.

For the forecasting CLI:

```sh
python3 forecast.py example_forecast.json
python3 forecast.py example_forecast.json --simulate-days 10 --dt 225
```

The certificate can exclude reactivation within its stated assumptions. A failure to exclude is not a prediction that reconnection succeeds. The simulation uses illustrative default parameters unless replaced explicitly in the input file.

## Main records and figures

- `compute/records/spatial_coverage_check.json`: sparse/dense and uniform/footprint comparisons.
- `compute/records/spatial_pulse_check.json`: every flood in the 12-case spatial-coverage timing scan.
- `compute/records/spatial_independent_reference.json`: adaptive-reference comparisons.
- `compute/records/accepted_spatial_histories.json`: higher temporal resolution histories used in the main figure.
- `compute/records/spatial_environment_sensitivity.json`: bounded external-heat and heat-transfer sensitivity.
- `compute/figures/spatial_coverage_and_reactivation.png` and `.pdf`: primary current figure.

Legacy uniform-layer comparisons and invalid early hydraulic runs remain in the separate local research archive. This branch packages the current footprint model and its explicit uniform-layer comparison in spatial_coverage_check.

## Failures and corrections

1. Initial hydraulic code used head as the grain-scale flow depth. Corrected to critical depth `2H/3`, then reran. Invalid old runs are segregated.
2. A zero-conductance branch produced an infinitesimal nonzero flux through a floating-point floor. Corrected to exact zero and verified the insulating limit.
3. The uniform insulating layer exaggerated sparse-rock shielding. Replaced by explicitly computed grain-footprint coverage and retained the counterexample.
4. The model still represents hydraulic obstruction through a mean crest. Spatial bypass between grains is the principal unresolved physical test.
5. Independent adaptive integration produced tiny negative numerical grain stocks, no larger than roundoff/error tolerance; these values are reported in the reference records, not claimed as physical negative volumes.
6. Exploratory JSON uses null for absent event times; the command-line tool refuses to emit nonstandard NaN.

All campaigns are exploratory. Do not label their case selection as preregistered or quote local melt contrasts as a glacier-model accuracy improvement.
