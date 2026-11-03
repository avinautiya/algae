# Actual field campaign handoff

Read `FINAL_FIELD_REPORT.md` and `review/FINAL_SCIENTIFIC_EVIDENCE_VERDICT.md` before running anything. The project does **not** meet all-baseline field superiority or novelty goals. Preserve that result. This additive branch does not modify the active Claude cliff/chemistry worktrees or prior frozen results.

## Executed

- Decoded/checksummed native Pirámide2022–23 thermal/moisture/weather MCOS tables.
- Matched camera distance to ice ablation and verified water-equivalent units independently.
- Preserved V1 physical no-go and V2 unscorable native endpoint.
- Frozen V3 eligibility amendment: one native camera endpoint moved three hours earlier by a deterministic missingness-only rule, before its value was read; all models/gains unchanged.
- Scored39 locked predictions at one initial period, then all39 unchanged predictors at three later source-eligible periods, publishing predictions before target release.
- Mesh/time/deep-boundary checks, source audits, source-bounded scenarios, exact energy accounting, conservative shared-reading error arithmetic and reproducible figures.

## Verdict to preserve

Three later-period MAE in mm water equivalent: forward candidate14.220367; derivative OLS14.349710; constant rate18.173057; buried-temperature index14.192518; air-temperature index9.870827. Including the earlier result, candidate11.754572 versus air7.906105. The nominal primary loses to the strongest baseline. No parameter scenario supplies established all-baseline superiority. Camera precision remains unknown; no confidence intervals or independent-glacier claims are justified.

The raw primary needs a training gain near0.235, so do not call the empirically corrected output an independently validated physical energy budget. Raw solver conservation is verified separately. The derivative OLS baseline is not source NTLS or the full professional energy-balance model.

## Where to start

`model/frozen_piramide_future_batch_final/melt_comparison/` holds predictions, seal, primary results and secondary four-period summary. Its parent freeze, source inputs and parameters are all pinned. The role receipts bind each camera release to a published prediction lock. See `data/camera/` for original source MATs and unit audits. Python dependencies are pinned in `data/requirements_decode.lock`; the local venv is intentionally not included. Run `python3 verify.py` after installing dependencies.

Figures are reproducible with `make_field_figures.py` and `make_future_figures.py`. No model is fitted by either plotting script. Rendered PNG/PDF figures are in `figures/`.

## Further work requires a justified new method

Do not use the four now-inspected periods as fresh validation of variants tuned to improve their errors. Do not rename a secondary capacity/geometry scenario as primary or remove the winning air baseline. Do not manufacture source camera precision, fill targets, clip negative melt or reinterpret source units to improve scores.

`contact_research/` is a separate mathematical and numerical foundation, not a validated new field model. Contact-resistance physics has prior art. The training-only steady check finds that a constant resistance large enough to explain the full raw amplitude mismatch would demand a near-base temperature jump inconsistent with source temperatures under those assumptions. Read `NEXT_CONTACT_NOVELTY_GATE.md` and `FOUNDATIONS_AND_GATE.md`; the recommendation is no production contact campaign without a stronger identification argument.

Tower2 metadata/features are audited in `contact_research/data/`. **Its camera outcome values remain unopened.** They can support an independently declared stake test if a new hypothesis survives prior-art, geometry, state/parameter identifiability and source-support checks. Sharing a glacier/weather record does not make it an independent climate replication. The unsupported February10–24 Tower1 camera outcome also remains unopened.

No Nature publication, STS placement, universal glacier-model superiority or new mathematical theorem has been established by this work.
