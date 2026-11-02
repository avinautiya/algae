"""Render corrected proposed release documents from the versioned registry."""
import json
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def documents(r):
 f=r['facts'];v=f['mechanism_versions']
 intro=f"Proposed additive correction, revision `{r['revision']}`, based on pinned commit `{r['audited_commit']}`. Historical documents and results remain unchanged. See `CLAIM_REGISTRY.json`, `SOURCE_MANIFEST.json` and `check_claims.py`.\n\n"
 table='# Corrected claim–evidence table\n\n'+intro+'Evidence classes and statuses distinguish recomputed arithmetic, reported numerical verification, retrospective comparison, conditional demonstration and unsupported independent validation. Reported results have not been re-run by this correction. **No claim establishes field superiority or independent coupled-cliff validation.**\n\n| ID | Claim | Evidence class | Status | Qualification |\n|---|---|---|---|---|\n'
 for c in r['claims']:table+='| '+' | '.join(str(c[k]).replace('|','/') for k in ['id','claim','evidence_class','status','qualification'])+' |\n'
 table+='\n## Source/version rules\n\nUse the raw-trace rescore for legacy ROI/recovery/history statements. Use solver versions/addendum for operator and boundary corrections, not as a completed full-v4 matched campaign. M3_FINAL is historical development evidence only where superseded. Different uncertainty calculations, endpoints and solver versions must not be pooled.\n\n'
 table+=f"The source-reported removal classes are3maintained/1genuine reconnection/2failed, with2pre-loss histories and only1solid-strip loss comparison. The contact-persistence boundary is{v['contact_persistence_boundary_days'][0]}–{v['contact_persistence_boundary_days'][1]}days under the1-day rule; the strength threshold is about{v['strength_defined_boundary_days_approx']}days. A day10removal rewets for2.6days and loses contact again at day12.4. None is permanent recovery.\n"
 memo='# Corrected decision memo\n\n'+intro+'''## Decision

Do not draft a new-environmental-mechanism or superior-field-predictor headline from the present package. A conditional mechanism study is a candidate, not a publication-ready conclusion. The failed field verdict remains visible. A future suitable public dataset is not ruled out merely because the present package fails eight metadata gates.

## Evidence that must govern the decision

- The candidate fails the declared all-baseline consistency target. It systematically overpredicts the six inspected test intervals. Recalibration against those test outcomes would not repair independent evidence.
- Calibration leverage is demonstrated instability: one interval dominates and leave-one-out scales vary widely. It is not demonstrated causal attribution of test bias. A leverage-equivalent count is not an independent-sample effective size.
- The historical0.02-day strip-loss agreement is one shared pre-removal trajectory under common thermal closures. The old0.37-day numerical GCI is a development loss estimate, not a statistical interval or a bound on every recovery forecast.
- The latest rescore separates maintained contact, genuine reconnection, temporary reconnection with reloss, and recovery strength. Its contact and strength boundaries answer different questions.
- Current solver source isv4, but complete authenticated v4 matched-ROI scoring is absent. Legacy locked-strip ROI errors are−54%and−67%(−0.0020/−0.0034m³), well above the5e−6m³ storage precision. Do not claim all melt predictions are accurate or dismiss these errors as rounding.
- Model-to-model agreement shares closures and assumptions. It does not independently validate those closures, porous clasts, vertical cliff melt, lifetime or Himalayan transfer.

## Candidate question, not an established discovery

'''+r['candidate_question']+'''

Peak runoff sufficient to rewet a route and integrated differential incision sufficient to retain access after the pulse are separate conditions. Contact at the return instant is not persistence across a subsequent observation window. Define that window and the heat/contact-strength endpoint before running a new comparison.

The broad facts that streams maintain cliffs, abandonment permits reburial, and debris limits circulation already have precedents. The independent novelty review found an evolving-topography/debris-transport comparator by Strickland, Covington and Gulley2025; a finite-head/finite-water-heat distinction is a comparison opportunity, not proof of priority. See the companion novelty review; the source documents here do not establish historical originality.

## Go/no-go sequence

1. Correct the release statements and metadata gate. Finish or explicitly keep blocked the v4 matched-ROI scorer. Preserve old results rather than overwriting them.
2. Check phase-generated water locally per route, or integrate a separately versioned phase source and verify full coupled balances and bed-update convergence. Fluid-only conservation is insufficient.
3. Compare against an evolving-topography baseline with matched initial head–discharge relation, partition and thermal power. A fractional-cover barrier with a different initial hydraulic rating is a weak comparator.
4. Freeze a bounded history/geometry challenge only after a falsifiable distinction and physically possible recovery condition are derived. Retain all outcomes. If the stronger baseline explains the results equally well, reject the new-method claim.
5. Independent observations remain necessary for environmental accuracy. Prepare data requests if authorized, but do not send messages or assert that missing measurements already exist.

**No-go today:** field superiority, unconditional irreversibility, all-melt accuracy, a new mathematical theorem, or Nature/award readiness. **Conditional go:** a narrowly falsifiable quantitative mechanism after the repairs and prior-art distinction. No new campaign was launched to write this memo.
'''
 budget='# Corrected signed error assessment — not a closed causal budget\n\n'+intro+f'''## Observed summaries

The six nominal signed errors are all positive. Their mean bias and mean absolute error are both **{f['nominal_pooled_bias_m']:.12f}m**. The value **{f['linear_imputer_marginal_bias_over_27_scenarios_m']:.12f}m** is a linear-imputer marginal over27nuisance scenarios, not the nominal case. The engineering scenario range+0.52…+0.82m is not a confidence interval. Preserve original units, calibration convention and scenario aggregation in every figure.

## One sign convention

Use `delta_error = (prediction_after - observation_after) - (prediction_before - observation_before)` for one specified counterfactual. Positive worsens the current positive error; negative reduces it. State what changes and what is held fixed. Do not mix an omission's current bias with the effect of adding the omitted process. This table records separate sensitivities/hypotheses; its rows cannot be summed into a causal attribution.

| Item | Explicit counterfactual | Sign of delta_error | Evidence and qualification |
|---|---|---|---|
| Calibration influence | Refit on different training subsets only | varies | Leave-one-out scale3.6–9.6 vs4.4 is demonstrated instability. No demonstrated attribution of held-out bias; leverage count{f['leverage_concentration_equivalent_count']:.6f}is not statistical ESS. |
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
'''
 return {'CLAIM_EVIDENCE_TABLE.md':table,'DECISION_MEMO.md':memo,'ERROR_BUDGET.md':budget}

def main():
 r=json.loads((ROOT/'CLAIM_REGISTRY.json').read_text());out=ROOT/'proposed_release';out.mkdir(exist_ok=True)
 for name,s in documents(r).items():(out/name).write_text(s)
if __name__=='__main__':main()
