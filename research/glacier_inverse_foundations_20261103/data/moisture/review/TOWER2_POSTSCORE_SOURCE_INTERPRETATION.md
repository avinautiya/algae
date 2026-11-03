# Tower2 released-source verification and interpretation

The training audit followed normalized-prediction publication `1fc6a5780523ba7b1032e08a0ec818ef3b38b4a1`; evaluation followed final-prediction publication `66cc20de0b0f9317f69cbe8d851bfd74ec489487`. The role-gated source audit validates the respective receipts and reads only their fixed endpoints. Frozen model code and observation operators were unchanged.

| Interval | Raw distance increment×917, mmwater equivalent | Supplied ablation increment×1000 |
|---|---:|---:|
| Dec2–Dec16TRAIN | 79.359886028146 | 79.359886028146 |
| Dec16–Dec30 | 54.256457338122 | 54.256457338122 |
| Dec30–Jan13 | 48.880581238058 | 48.880581238058 |
| Jan13–Jan27 | 66.179102025554 | 66.179102025554 |
| Jan27–Feb10 | 75.930078267978 | 75.930078267978 |

All endpoints are at03:00in the source's naive local clock. Unit differences are at floating-point precision. The supplied ablation field is consistent with metres water equivalent, despite the README's millimetres description. The predeclared raw exposed-distance conversion is correct and unchanged; the failed performance cannot be explained by a1000-fold source-unit mistake.

No complete interior camera series was inspected in this audit. Positive endpoint differences do not establish absence of interior reference resets or validate manual corrections. Camera precision remains unquantified. No numerical confidence interval or guaranteed field-error bound follows. Tower3and the excluded Tower1February24endpoint remain unopened.

## What the held-out comparison supports

The published results give exponential-moisture candidate periodMAE13.353386mmwe versus7.265305for airPDD; the air model has smaller nominal error in every one of the four periods. All candidate predictions overestimate the observed increments. The predeclared improvement test fails.

This is a real held-out stake comparison with one local training interval and fixed development parameters. It is evidence against this **specific transferred exponential moisture correction**, not against all moisture physics. It does not establish atmospheric energy-balance superiority, Himalayan transfer, glacier-scale skill or Nature-level novelty.

The invariant model exactly cancels a constant additive sensor offset, but actual performance shows that algebraic invariance did not produce better transfer here. Soil calibration drift, multiplicative calibration, sensing-depth differences and other structural effects remain possible explanations, not demonstrated causal diagnoses. In particular Tower2SM3is nominally6cmaboveice, versusTower1SM3at the interface; raw sensor trends are different spatial observables. Those limitations were documented before release and must not become a post-hoc rule for excluding unfavorable periods.

A retrospective probe change, differentbeta, negativebeta, gradient-moisture blend, period deletion, or calibration adjustment made using these evaluation errors would turn Tower2into development evidence. It would require a different independent evaluation before claiming restored superiority. Preserve all nine models, previous failed studies and source flags.

## Reproduce

From the round3root with the parser environment available:

```
../glacier-round2-20261103/data/venv312/bin/python data/moisture/review/postscore_tower2_source_audit.py train --output data/moisture/review/TOWER2_TRAIN_SOURCE_AUDIT_REPRO.json
../glacier-round2-20261103/data/venv312/bin/python data/moisture/review/postscore_tower2_source_audit.py evaluation --output data/moisture/review/TOWER2_EVALUATION_SOURCE_AUDIT_REPRO.json
```

The wrapper locates the released Pirámide package through the existing sibling/environment locator. Receipts, source hashes, unit checks and unchanged-operator flags are retained in both auditJSONfiles. Output files must be new; existing evidence is never overwritten.
