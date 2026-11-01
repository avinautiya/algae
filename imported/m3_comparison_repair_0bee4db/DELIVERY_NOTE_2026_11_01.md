# Delivery correction and current-solver interpretation

This repair package previously existed only in the Codex local workspace and
downloadable archive. That path was not shared with Claude's remote container.
It is now published additively as `m3_comparison_repair/` on a separate branch.
Do not search for it in the active branch until it has been fetched/imported.

All 41 tests passed again from a clean extraction on 1 November 2026. They are
component/scoring checks; they launch no CFD campaign and do not validate the
current production solver automatically.

## Same clamp defect, not a competing patch

The package identified BOTH mutation of non-ice tail cells and heat-implied
melt without corresponding geometric ice lowering. V4's unclamped operator
addresses those same defects. Do not revert it to the pinned snapshot or
blindly replace it with this package's independent diagnostic operator.
Retain its fixing commit and its actual phase/geometry checks.

`snapshot/` is intentionally pinned to legacy commit
`f3602e35e1595bcdbee4208a354a34944f1ae974` for reproduction. It is not production
source. Historical classes and diagnostic numbers in the package are not v4
results.

## Run from this folder

```sh
python3 verify.py
python3 -m melt.repair_eval /ACTUAL/ORIGINAL_MEDIUM_RESULTS snapshot/frozen/production_spec.json --output /SEPARATE/matched_melt_original.json
python3 classification/recovery.py /ACTUAL/ORIGINAL_RUN.json --output /SEPARATE/classified_original_run.json
```

The scorer validates the legacy schema/config/initial profiles. An incompatible
new export must be rejected or handled with a separately versioned adapter;
never weaken its checks merely to get a score. Preserve original raw results.

`classification/recompute.py` regenerates REPORT-DERIVED legacy accounting and
frozen-prediction classifications. It is NOT a raw-v4 classifier or an independent
recomputation of current shared histories. For current runs, use `classify_run`
and `history_key` with actual causal input/provenance records. Keep the raw
one-day sustained-contact label distinct from newer R-based partial/full
recovery-strength labels; they answer different questions.

## Finest-grid result wording

Warm-started fine hydraulics over a transferred medium-grid bed history tests
the post-removal hydraulic crossing on that history. It does not establish full
coupled morphodynamic convergence. Preserve the reported 0.02-day shift as a
conditional sensitivity result, and keep unresolved near-boundary fractions
explicit. Do not rerun expensive fine campaigns merely to import this folder.
