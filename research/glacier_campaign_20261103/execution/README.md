# Bounded execution preflight

Audited source: repair commit `f1f85ee496dd39617fedf0f58202b943b3af3168`. The old private directory does not survive in this environment; source/schema were inspected directly from that pinned GitHub commit. The parent subsequently restored its verified release under the campaign inputs folder.

```sh
python3 -m unittest -v test_preflight
python3 preflight.py --repo /absolute/path/to/ISOLATED-checkout --out /absolute/path/OUTSIDE-checkout/NEW-preflight --mode physical_bare
```

For a deliberately conditional quasi-steady debris scenario use `--mode conditional_quasi --covered`. This does not make it physical-time verified. `--mode physical_bare --covered` refuses immediately.

The runner verifies every repair-manifest hash and the two installed phase-solver copies, copies the packet into a NEW output directory, then runs its66 pinned repair tests and an explicit solid/porous physical-source atomic-refusal probe. It never edits a frozen result or Claude checkout. Its logs and regenerated test report belong to the copy. Source is rechecked after execution. Existing output directories refuse.

Resource controls: one thread per BLAS/Numba/OpenMP process, a600-second wall limit per stage by default, and process-group termination on timeout. Two stages mean at most2×the configured timeout; this is not a RAM quota or a production scheduler. A free-space gate defaults256MiB; actual production runs need a separately measured budget. No cleanup of shared data is attempted. The checker records a GitHEAD when available rather than falsely asserting that an unversioned snapshot is the pinned checkout.

**Execution blocks**

1. Hash mismatch, fewer/more than66 expected pinned tests, failed source-refusal probe, insufficient disk or stale output: stop and create a versioned investigation; never edit a freeze to silence a gate.
2. Native physical-time mode does NOT support covered meltwater, porous storage or the uniform-incision counterfactual. Do not bypass refusal to complete a campaign.
3. A quasi-steady covered simulation must retain its drainage/time assumptions and separate ledgers. It does not obtain the physical-time certificate by passing software tests.
4. New long-run event-time and morphology-step convergence are still required. Existing two-second smoke results are not evidence of21-day prediction accuracy.
5. The current field-superiority result is failed. No preflight or metadata result changes it, and no output is labelled a novel environmental discovery.

Five independent preflight regression tests pass. The complete66-test real-source run is not claimed here unless an actual `PREFLIGHT_RESULT.json` records it.
