# Corrected signed error assessment — not a closed causal budget

Proposed additive correction, revision `corrected_evidence_v1`, based on pinned commit `05e7312e3bbf248acf9f646267a5b8ea06a8b32e`. Historical documents and results remain unchanged. See `CLAIM_REGISTRY.json`, `SOURCE_MANIFEST.json` and `check_claims.py`.

## Observed summaries

The six nominal signed errors are all positive. Their mean bias and mean absolute error are both **0.584872726098m**. The value **0.630937225384m** is a linear-imputer marginal over27nuisance scenarios, not the nominal case. The engineering scenario range+0.52…+0.82m is not a confidence interval. Preserve original units, calibration convention and scenario aggregation in every figure.

## One sign convention

Use `delta_error = (prediction_after - observation_after) - (prediction_before - observation_before)` for one specified counterfactual. Positive worsens the current positive error; negative reduces it. State what changes and what is held fixed. Do not mix an omission's current bias with the effect of adding the omitted process. This table records separate sensitivities/hypotheses; its rows cannot be summed into a causal attribution.

| Item | Explicit counterfactual | Sign of delta_error | Evidence and qualification |
|---|---|---|---|
| Calibration influence | Refit on different training subsets only | varies | Leave-one-out scale3.6–9.6 vs4.4 is demonstrated instability. No demonstrated attribution of held-out bias; leverage count1.339480is not statistical ESS. |
| Imputer | Change imputer, keep the published per-scenario fitting convention | varies | Reported marginal bias0.60–0.66m. Marginals are not paired nominal perturbations or independent causes. |
| Clock | Apply reported±2h clock scenario and its fitting convention | varies | Marginal bias0.56–0.75m; largest reported single-factor shift about0.19m. True clock correction not established. |
| Temperature offset | Apply±0.05°C scenario | varies | Reported marginal bias0.61–0.65m. Does not demonstrate actual sensor bias. |
| Width | Change1.5–2.5m, recalibrating each scenario as originally done | near zero | ≤0.005m reported shift is confounded by re-fitting; not an independent geometry constraint. |
| Reconstructed support | Replace missing forcing under declared rules | unknown causal sign | Test reconstructed share≤4.1%is descriptive. A small duration share alone does not bound its error contribution without forcing/response bounds. |
| Pressure depth | Increase/decrease assumed depth with closure held fixed | regime dependent | ±0.1m gives arithmetic heat-transfer changes about13–18%. Independent water-file station/time binding and floor interpolation are hypotheses; mismatch under that mapping does not identify which source is wrong. |
| Clearance/bulk temperature | Replace proxy by independently measured clearance or bulk relation | unknown here | Physical hypothesis; source values and whole-pipeline sensitivity are unavailable. |
| Cold ice | Add a physically consistent cold-ice energy sink, same latent target | negative if all else fixed | Less heat available for melt. Unknown ice temperature; current model does not implement the full process. |
| Viscous dissipation | Add positive dissipative heat to prediction, same observation | positive if all else fixed | Would add melt, contrary to a correction that reduces overprediction. Magnitude not computed; do not label its added-process effect negative. |
| Extra observed surface lowering | Increase observation by an independent non-modeled lowering component, same prediction | negative | Observation correction, not a direct thermal prediction effect. Magnitude and even source attribution unmeasured. |
| Discharge mismatch/open-channel closure | Replace assumed relation with an independently justified one | unknown here | Hypotheses; no test-outcome fitting permitted. |
| Survey uncertainty | Change target within its stated source uncertainty | both signs | Measurement uncertainty is not a demonstrated random realization or explanation of a particular residual. |

## Nonclosure

No individual tested sensitivity explains the whole reported bias within its stated range. This does **not** rule out interactions, a joint explanation, omitted inputs or structural errors. The budget remains unclosed; no total attributable bias or residual causal allocation is reported. Calibration leverage is a demonstrated diagnostic, not a proved causal mechanism. No coefficient is tuned to make the budget close.
