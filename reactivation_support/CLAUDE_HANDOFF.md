# Handoff while the main milestone continues

This commit adds only `reactivation_support/`. Import that folder or cherry-pick this single support commit; do not merge or reset the active work to the support branch. All source audits are pinned to c70c553 and may already be superseded by your current edits.

Read `DIAGNOSTIC_FINDINGS.md` first. The main experimental lead is entrance treatment/reference mismatch: in a matched frictionless fixed-bed test, a graded inlet changes the strip's reference discrepancy from +0.10475 to +0.01152 at dx=.05, then -0.00384 at dx=.025. Manning removal, face-flux extraction and changing only the outlet drop do not explain the original large discrepancy.

Use the smooth geometry as an applicable verification test, not a way to force the production experiment to pass. A fixed graded approach can become an artificial sill during incision; preserve a physically consistent moving-boundary treatment and justify the production entrance model.

Ready helpers:
- `hydraulic_reference.py`: ODE-free flat-channel GVF reference (~416× measured routine speedup; restricted geometry).
- `event_schedule.py`: exact removal and endpoint times.
- `campaign_analysis.py`: actual removal timestamps, sustained ordinary-flow reactivation, censoring and conditional refinement diagnostics.
- `connectivity.py`: winding-path geometric crest rather than straight columns; not a discharge solver.
- `phase_update.py`: verified local mass/energy exchange contract. Read the multirate caution before integrating its adapter; it is not a complete global-ledger fix.
- `plot_campaign.py`: figures from completed run files with explicit evidence labels.

The package's 46 helper tests pass. Pinned snapshots and development fixtures do not establish final-data eligibility. Keep your existing verification gate and freeze procedure.
