# Bed-operator / solver versions (reactivation experiment)
| version | defect | effect | status |
|---|---|---|---|
| v1/v2 | `lower_bed` clamped EVERY cell to pool_floor: the overfall tail floor was lifted to −1.0 m after the first bed update | hydraulics never settled (storage ≤3 % Q); early-removal late decline artefact | superseded |
| v3 | clamp restricted to ice cells; **still discarded melt** for ice cells at pool_floor (ramp bottom): geometric bed loss was 18–21 % smaller than the heat-implied volume after 21 d (coarse strip none: 1.925 vs 2.446 m³; bare: 2.360 vs 2.882) | ledger gap in the bed operator; reach/footprint melt unaffected (reach cells never reach the floor) | superseded by v4 |
| v4 | no clamp; `ledger_gap` recorded each update; unit test `TestBedOperatorLedger` (identity to 1e-14, also for cells sitting at the old clamp level) | coarse 21 d: geometric vs heat volume rel. diff 1.8e-9 / 8e-10; f_toe differs from v3 by ≤9e-4 (no removal) and ≤4.8e-3 (removal at 8 d, post-removal); loss day +0.011 d (2.7212 → 2.7322) | current |
The attachment `m3_comparison_repair` (CLAUDE_HANDOFF.md, verify.py, matched-melt scorer) was NOT available in this environment (searched the filesystem and all remote branches); the repair above was derived independently from the heat-vs-geometry volume check.

## v4 finest-grid boundary repeat (fixed bed, warm-started dx=0.0125; `results/v4/fine_warm_day*.json`)
| removal day | dx=0.025 | dx=0.0125 |
|---|---|---|
| 9.5 | 0.0241 | 0.0291 |
| 10.0 | 0.00677 | 0.00933 |
| 10.5 | <1e-6 | <1e-6 |
Crossing of f=1e-3: 10.43 d (0.025), 10.45 d (0.0125). Medium v4 strip-none: geometric/heat volume gap 1.4e-9, loss day 2.4828 d (v3 2.4775).
Not done (inputs unavailable): matched-melt scorer, history-group recomputation and `verify.py` from the missing `m3_comparison_repair` folder.
