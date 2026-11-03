# Independent post-score audit — frozen representative-point plate diagnostic

**Verdict: the declared primary superiority rule fails.** The retained diffusion model is worse than the strongest declared comparator at both sites. Both numerical gates pass. This is a valid negative result for the declared conditional raw-signal prediction problem, not evidence of improved glacier melt prediction.

The pre-execution commit is `2a0828dd2a4c8c42562a277eec8dfb2b991cf035`, prefix `research/independent_evidence_20261103/extension/`. I fetched the remote freeze and three source files at that commit, verified their Git blob hashes against local bytes, then checked local source/data SHA256 values against the frozen manifest. The audited owner source is unchanged. I did not change a parameter, refit a model, select a different test window, change sign, or modify owner results.

## Independently recomputed results

| Site | Retained hours | Diffusion MAE W/m² | Causal ridge MAE W/m² | Instantaneous gradient MAE W/m² | Training median MAE W/m² |
|---|---:|---:|---:|---:|---:|
|1|336|1.2599108201|1.0076641038|4.9608209499|10.1835372619|
|2|336|1.0729556760|0.6854185275|2.8161612217|8.9121706548|
|Equal-site mean|672 total|**1.1664332480**|**0.8465413156**|3.8884910858|9.5478539583|

Diffusion's primary equal-site MAE is **37.79% larger** than ridge. Its wins over the median and instantaneous-gradient comparators are secondary results; they do not reverse the frozen primary failure. The frozen rule requires diffusion better at both sites, at least5% lower aggregate MAE, successful numerical gates and a non-endpoint alpha. Only the latter two requirements pass.

## Exact checks performed

* Verified all672 site/timestamp pairs are unique and match the complete finite-plate source support on the frozen Nepal-time interval. Each site has336 rows. No difficult records were dropped.
* Matched every saved observed raw plate value exactly to the source-aligned data. No target conversion to latent melt was used.
* Independently recomputed every saved model/site n, MAE, RMSE and signed bias, including the predeclared opposite-orientation sensitivity. Maximum discrepancy is3.55e-15.
* Replayed predictions with **saved fixed parameters**, constructing endpoint features, training median and instantaneous gradients independently. Replayed diffusion uses the frozen solver with the same saved alpha/coefficient; it is not a new discretization or physical validation. Maximum prediction discrepancy is7.55e-15. No fitting function was called.
* Checked that saved alpha is the minimum of the retained61training objectives and is not a search endpoint. Both sites use that same alpha; their plate amplitudes are locally fitted on their own training records as declared.
* Recomputed64versus128mode differences using the same fixed coefficient: site1maximum0.0230623W/m², site2maximum0.0185384W/m², both below the frozen0.5W/m² thresholds. This does not establish physical model adequacy or unique parameter identification.
* Recomputed the frozen decision from saved predictions and gate values. It is `false`, matching the owner's result.

The executable audit is `audit_plate_scores.py`; `POSTSCORE_AUDIT.json` records checks and metrics. The script assumes the restored source-aligned archive path used during this session and can be adapted to a new input location without changing scientific parameters.

## Interpretation limits

This experiment predicts archived signed plate signal under observed endpoint temperatures using a representative-timestamp approximation. Native averaging interval, physical orientation/placement and sensor migration remain unresolved. There is no demonstrated latent melt, native interval-flux, weather forecast, unseen-glacier transfer or professional-model superiority. The two nearby sites and their correlated hours are not independent glacier replications.

The positive numerical result is that source/data integrity, fixed support, arithmetic and code replay are consistent. The scientific comparison remains negative. Preserve both; do not tune the evaluated window or change the primary comparator to manufacture an improvement.
