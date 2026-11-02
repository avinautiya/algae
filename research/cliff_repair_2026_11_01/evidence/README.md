# Corrected evidence release, revision1

This is an additive proposed repair, not another audit-only memo. It contains corrected replacements for three release documents, a versioned claim registry, eleven authenticated source snapshots, and an executable regression gate. No Claude worktree, old document, solver output or freeze was changed.

Source repository: `avinautiya/algae`; pinned integration revision: `05e7312e3bbf248acf9f646267a5b8ea06a8b32e`. `SOURCE_MANIFEST.json` records each original path, exact Git blob hash and source URL. The original documents remain under `source/`; the corrections are under `proposed_release/`.

## Use

```sh
python3 build_docs.py
python3 -m unittest -v test_claims
python3 check_claims.py --target evidence-consistency --out EVIDENCE_CHECK.json
python3 check_claims.py --out V4_GATE_RESULT.json
```

The first check exits0 when the corrected claims match the pinned evidence. The second deliberately exits2: a complete authenticated v4 matched-ROI campaign/scorer is absent. This is the expected fail-closed result, not a broken test. Twelve regression tests pass. No production simulation or recalibration is launched.

## Actual repairs

- Nominal bias/MAE is recomputed from six signed errors as0.584872726098m;0.630937225384m is labeled as the linear-imputer marginal instead.
- Latest raw-trace rescore governs recovery classification/history and legacy ROI claims. Contact persistence and recovery strength remain separate endpoints; day10rewetting followed by reloss is not permanent recovery.
- The0.02-day strip-loss agreement is one pre-loss history, not three independent validations. Old loss GCI is not a statistical confidence interval or a universal recovery uncertainty.
- The table explicitly keeps the failed superiority verdict, unavailable strict-native score, unverified porous closure, unavailable v4 matched scoring and unsupported field/lifetime claims.
- The error assessment uses one explicit change-in-prediction-minus-observation convention. Its independent sensitivities/hypotheses are not summed or forced to close. Calibration leverage is instability, not established test-bias causation. Package-specific data gaps are not a proof that future public data cannot exist.
- Corrected documents are generated from the registry and tested byte-for-byte against it, preventing drift between the machine record and the proposed prose.

## Limits and migration

The checker recomputes actual nominal arithmetic from pinned JSON. Mechanism figures are checked against authenticated report text, not rerun from the full raw traces or a new full v4 campaign. These distinctions are in every output. It does not validate source truth, historical novelty, physical applicability or journal/award readiness. A future source or solver revision needs a separate registry, evidence manifest and scientific review; do not toggle a gate to true.

The candidate question is temporary rewetting versus sustained post-pulse heat access under finite head and heat-limited evolving geometry, tested against a hydraulically matched evolving-topography baseline. This is a falsifiable comparison opportunity, not a guaranteed novel result. Read the companion mathematical and prior-art reviews before freezing any new campaign.
