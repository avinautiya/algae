# Continuous calibrated-goal glacier heat inference

This round develops an executable melt/heat-inference algorithm and audits new field data before using native melt outcomes. It follows the mixed adaptive-graph development results; it does not relabel them as independent validation.

## Algorithmic contribution under investigation

A thermal model can fit temperatures and one observed melt interval while leaving later melt poorly constrained. The desired quantity is the calibrated heat ratio, not necessarily unique material parameters. For an explicitly justified linear branch, the new module optimizes that ratio over a continuous feasible parameter polytope using Charnes–Cooper linear-fractional programming. It returns numeric lower/upper melt estimates, extremal witness models, a conditional minimax point and an abstention flag. Coefficient-box uncertainty is also supported, without claiming the physical errors are independent.

The improvement over the earlier finite-family checker is continuous optimization over the declared polytope, rather than a scenario grid that can miss feasible cases. Independent vertex enumeration on100random2Dproblems agrees within4.55e−13mmwater equivalent. This is an executable numerical result, not evidence of environmental prediction accuracy.

Charnes–Cooper optimization, inverse heat conduction, goal-oriented inference and minimax intervals are established mathematics. The research novelty question is the specific source-error-constrained calibration-transfer algorithm for glacier melt and its demonstrable environmental value. A Nature-level discovery, superior field accuracy and STS placement are not established by these derivations or numerical tests.

## Run the prediction tool

```
python continuous_goal/predict.py continuous_goal/SYNTHETIC_EXAMPLE.json --output NEW_OUTPUT.json
```

The example is explicitly synthetic. The output is conditional inference, not automatically an advance forecast. The CLI refuses absent unit/phase/feasible-set declarations; declarations themselves are not scientific validation. A nonlinear unknown phase history or an invalid physical feasible set must not be hidden behind a numerical interval.

`continuous_goal/DERIVATION.md` gives the transformation, conditions, uncertain-coefficient extension and limits. `review/` contains independent numerical and primary-source reviews. `algorithm/` generalizes weak-energy inference to the three-sensor thermal source rather than inventing additional sensors.

## Field-data status

The user explicitly authorized Lirung ablation-file download after automatic review had rejected the earlier attempt. Published size/MD5 match. Schema/site/date metadata were inspected while native numerical melt magnitudes were withheld from fitting and scoring. Published paper abstracts have exposed aggregate melt context; no claim of fully blind prospective data collection is made.

The native file contains seasonal totals, not daily targets. Source interval, unit and stake/profile matching are insufficient for the previously proposed chronological test. Forcing coverage and site differences cannot be silently filled or scored as complete paired observations. All eligibility failures and conditional alternatives are retained. A separate metadata audit is investigating longer Kennicott observations before any outcome release.

No newly inspected native melt target is used to select a candidate, exclude an unfavorable test or claim superiority. Previous published branches/results remain preserved. New authorized metadata are additive; old release snapshots are restored to their published hashes.
