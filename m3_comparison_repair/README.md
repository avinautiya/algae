# M3 comparison repair

Independent additive audit of pinned commit
`f3602e35e1595bcdbee4208a354a34944f1ae974`, 1 November 2026.
Claude's source, jobs, checkpoints and frozen results are untouched.

Read `CLAUDE_HANDOFF.md`, then run `python3 verify.py` (41 tests).

- `melt/`: common central 1 m x 0.5 m observation support, saved-profile
  reconstruction, scope-matched prediction scoring and safe input validation.
- `classification/`: maintained contact versus actual recovery versus failed
  recovery, persistence/censoring checks and shared pre-loss history accounting.
- `numerics/`: conditional GCI reproduction and a gated isolated finest runner.
- `bed_audit/`: zero-time bed mutation reproduction, conservative ice-only
  operator and a bounded paired coarse development check with saved results.
- `snapshot/`: ten exact pinned upstream blobs; provenance records their Git IDs.

The report-derived nonporous removal cases are two maintained, one genuine
reconnection and one failed recovery. E04/E06/E07 share one pre-loss history.
The old GCI is 0.3672 days, a conditional numerical estimate, never a CI.

The pinned global bed clamp raises non-ice tail cells at zero elapsed time and
can count latent melt without actual lowering. The paired coarse 3.5-day check
shifts loss by 0.0237 day and removes 0.0166 m3 of counted-minus-geometric melt.
It does not quantify the correction at medium/fine resolution or over 21 days.
The fine runner refuses that defective pinned operator; repair/refreeze first.

Raw M3 CFD JSONs are absent from the complete public pinned tree. No corrected
21-day melt ranking is claimed. Claude can run the scorer on its existing files:

    python3 -m melt.repair_eval /ACTUAL/results/eval snapshot/frozen/production_spec.json --output /SEPARATE/matched_melt.json

This needs no CFD rerun. Shared thermal closure, physical approximation and
model accuracy remain distinct from observation-support correction.

See component reports for exact assumptions, tests, raw-data limitations and
runtime estimates. Dependencies: NumPy, SciPy, Numba; Matplotlib for figures.
