# Post-release independent result audit

`audit_actual_results.py` must run **only after the parent announces that evaluation outcomes are released**. It never reads camera files. It reads the resulting observation totals from the scored result artifact, then independently checks arithmetic and sources. The explicit `--evaluation-released` flag gates every artifact read. No actual calibration/prediction/result artifacts were opened during preparation.

Example after release, with actual paths substituted:

```
python3 audit_actual_results.py --evaluation-released \
  --normalized NORMALIZED.json --predictions PREDICTIONS.json \
  --seal PREDICTION_SEAL.json --results RESULTS.json \
  --spec ../analysis/ANALYSIS_SPEC.json --decoded-root ../data/decoded \
  --out NEW_INDEPENDENT_AUDIT.json
```

The checker verifies the complete15-model family, every thermal dose-ratio prediction, train-only gain arithmetic, normalized-heat/capacity units, fixed primary identity, same native917 target operator, independent PDD feature integrals and source hashes, observed/predicted totals, errors and nominal comparison labels. Alternative geometry/hold models remain secondary; they cannot replace the frozen primary after seeing its result.

For each primary-baseline pair it derives the conditional absolute-bound precision needed to identify the nominal ranking. If the true cumulative readings differ from their observed readings by e0,e1,e2, each bounded byδ, training/evaluation increments shift as e1−e0 and e2−e1. With ratios fa,fb, the rank-determining affine quantity is `(fa+fb)X−2Y`, with shared-reading radius `[2(fa+fb)+4]δ`. The reported threshold includes the positive-increment gate. It is converted to metres ice using the fixed917 convention.

That threshold is a **post-score diagnostic**, not an observed camera precision, newly selected model criterion, statistical confidence interval or retrospective pass rule. Source precision is still unknown. The script also distinguishes the already-issued fixed-forecast error, for which training input is not recalibrated when measurement errors vary.

Five synthetic tests pass, including all15-model recomputation, deliberately corrupted errors, an independent sign-crossing integral, release refusal before file access and shared-reading-corner verification of the rank threshold. Actual field success/failure awaits the authorized audit.
