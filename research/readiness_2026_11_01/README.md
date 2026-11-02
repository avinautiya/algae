# Independent pre-campaign readiness review

Audited repository: `avinautiya/algae`, branch `claude/integration-m4`, pinned commit `05e7312e3bbf248acf9f646267a5b8ea06a8b32e`. Three independent agents reviewed prior art, mathematics/physics, and numerical/environmental evidence. The parent added two isolated repair modules. No Claude source, job, checkpoint or historical result was changed. No production campaign was launched.

**Decision: NO-GO for a new-mechanism, field-superiority or Nature-readiness headline. Conditional GO for developing a precisely defined finite-head recovery boundary, after the concrete defects below are addressed.** This review is bounded; it does not certify exhaustive literature coverage or a competition outcome.

## Concrete findings

1. **The readiness validator is not fail-closed.** The committed checker accepts seven deliberately invalid metadata examples. It treats strings, nonfinite/negative uncertainties and empty containers as evidence; its documented numerics group is absent. `strict_contract/` supplies a typed version2 checker, blank template and12 tests. Complete syntax is explicitly NOT a certificate of scientific validity. Its output never authorizes an environmental-accuracy claim. This is a schema migration, not a silent drop-in accepting old truthy fields.
2. **The release uses superseded recovery evidence.** The latest rescore must supersede the earlier M3 final report. Repeated loss times share one pre-removal strip history. Maintained contact, genuine reconnection, eventual reloss and strength are distinct endpoints. The latest contact boundary is10.0–10.5days; the strength boundary is about9.2days. A0.02-day warm-start fixed-medium-bed grid shift is not full coupled uncertainty. Old0.37-day loss GCI cannot be assigned to the recovery boundary.
3. **The reported nominal bias is mislabeled.** Six nominal errors average0.5848727261m, not0.63m. The latter is the linear-imputer marginal across nuisance scenarios. Calibration concentration is not a statistical effective sample size or a demonstrated cause of test bias. Individual sensitivity shifts do not exclude a joint explanation. Error-budget signs mix changes to prediction, changes to observation, and omission bias; use one declared `delta(prediction-observation)` convention.
4. **Phase mass is missing from the current SWE bed operator.** It lowers the bed without adding generated meltwater to depth. Its fluid-only ledger does not close the combined ice–water system. The isolated `phase_source/` operator adds the density-consistent source while preserving zero-born sensible heat and horizontal momentum under stationary temperate ice. Eight tests pass and a separate agent checked the derivation. It is NOT integrated or a full moving-boundary momentum/energy solver. A source-step test is not verification of the full SWE solver.
5. **Fixed-discharge lock-out is not environmental irreversibility.** Within the ideal two-route rating, toe flow is zero when `g+d >= (Q/Kb)^(2/3)`. A larger inflow can restore contact. Our isolated example doublesQ and reopens the route, but it recloses whenQ returns. Durable reset requires sufficient differential incision, not merely temporary wetting. See `math/` for equations, assumptions, exact bounds and executed checks.
6. **General mechanism claims overlap existing work.** Stream abandonment/reburial and debris-limited circulation have direct precedents. A newly located closer process model couples channel incision, debris movement and internal drainage. The distinctive candidate is finite-head/finite-water-heat recovery under evolving topography, not the broad statement that debris arrangement affects flow. See `novelty/NOVELTY_RED_TEAM.md` for primary sources, equations and code evidence.

## The sharper, falsifiable research question

**Can a finite-head, heat-limited model predict whether a runoff pulse merely rewets an abandoned ice route or restores sustained heat access after the pulse?**

Two quantities must be separated: peak discharge sufficient for initial wetting, and cumulative differential incision sufficient for post-pulse access. For gap `g=b_toe-b_bypass`, the exact reduced identity is `g_dot=E_bypass-E_toe`. Returning to baselineQ0 permits positive toe flow only if the post-pulse gap falls below `Hc(Q0)`. Thus the required integrated repair is greater than `g_initial-Hc(Q0)` under the stated reduced assumptions. This elementary deduction is not asserted as historically new mathematics. Equal-area symmetric bare closures can prohibit durable reset; geometry and thermal footprints must be specified before predicting it.

The repair integral at the pulse end certifies access only at the return instant, not persistence afterward. For a declared subsequent interval, positive contact requires `g(t)<Hc(Q0)` throughout that interval; continued differential incision can cause reloss. This is exactly why temporary reconnection, duration and strength must not be conflated.

The novelty audit finds a concrete comparison opportunity in Strickland, Covington & Gulley2025, but does not prove priority. Their preprint/code includes evolving topography and debris transport; its stream-power incision lacks the candidate's explicit finite-head and finite-water-heat states. Batley & Smith2025 is an important closed-channel methods neighbor, not necessarily the closest environmental competitor.

## Minimum conditions before asking Claude for more production results

1. Correct the release using the pinned latest rescore, consistent endpoints and signed error summaries. Finish or explicitly gate the missing v4 matched-footprint scorer; do not transplant legacy accuracy claims.
2. Replace the truthiness validator with a typed metadata gate. Review evidence independently. Bound generated meltwater locally per route, or integrate a new versioned phase source and verify the complete coupled ledgers, positivity, time/space and bed-update convergence. No silent replacement of frozen outputs.
3. Define a stronger comparator with evolving topography and equal initial head–discharge/thermal support. An area-averaged barrier whose initial hydraulic rating differs is an unfair stand-in for current models. Remove unsupported porous-clast claims from the main claim set unless independently checked.
4. Fix one primary endpoint and a finite observation window before testing: e.g. sustained post-pulse heat delivery. Separately report first wetting, delivered heat, contact strength and subsequent reloss. A conditional0.02-day internal agreement cannot be the success criterion.
5. Only then freeze a small independent forcing/geometry challenge that could falsify the proposed boundary against a strong baseline. Predict all cases before execution and retain failures. If the stronger baseline explains the outcomes equally well, reject the novelty claim. Synthetic success still does not establish field accuracy.

## Reproduce the bounded checks

```sh
cd strict_contract
python3 -m unittest -v test_contract
python3 validity_contract_v2.py --template blank.json
python3 validity_contract_v2.py blank.json  # expected exit2
```

```sh
cd phase_source
python3 -m unittest -v test_phase_source
```

The agent audit folders retain their source manifests, independent reproductions and explicit limitations. No metadata can manufacture missing environmental observations, and no review can guarantee journal acceptance or award selection.
