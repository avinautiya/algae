# Release README — ice_integration (milestone M4 + follow-up)

**Status: verified computational integration and a retrospective, non-independent field comparison. The integrated model did NOT improve on every baseline consistently. Nothing here is independent validation.** All frozen results, the failed-superiority verdict, the original demo runs and superseded versions are preserved. No recalibration against evaluation intervals and no new production campaign were run for this package.

## Candidate repair branch — versioned changes

The new `validity_contract` is a typed metadata gate; exit0 means metadata complete, never environmental validation. The original untyped checker and historical Swiss result remain archived in `research/cliff_repair_2026_11_01/contract/legacy/` and Git history. See `VALIDITY_CONTRACT_V2.md` and the explicit migration report.

New `ice_cliff_survival.reactivation.swe2d_v5.Sim` provides a phase-source repair and physical-time split step. It is a separate, opt-in solver version; the older `swe2d.py`, frozen runs and demo commands below are not silently rewritten. Native physical mode currently supports bare connected temperate ice only. Covered sources, porous storage and uniform-incision counterfactuals refuse in that mode. Quasi-steady covered runs remain explicitly conditional and do not gain a physical-time conservation certificate. This phase-source v5 is unrelated to the earlier Swiss point-model V5 name.

The corrected decision memo, claim table and error assessment supersede the earlier release summaries. No failed superiority verdict has been reversed. The full-v4 matched-melt gate is still incomplete; the new fixed-ROI exporter is for future versioned records, not invented legacy data. Numerical source/ledger tests do not certify field accuracy, novelty or award/journal outcomes.

Reproduce repairs with `python3 research/cliff_repair_2026_11_01/verify_repairs.py`. The three short physical-time examples under that packet are verification, not a glacier campaign.

## Historical commands (verified before this repair)
| Command | Result |
|---|---|
| `python -m pytest ice_integration` | 15 pass |
| `python -m ice_integration.score_integrated OUT_DIR` | Reproduces Codex V5 predictions to 4e-16; recomputes per-station/pooled errors; prints the fixed-rule verdict (False) |
| `python -m ice_integration.demo_coupled OUT.json DX N [strict\|reconstructed]` | Coupled plan-view demonstration with imposed Q; resumable, atomic checkpoints (0.05 m: ~5 min for 56×6 h; 0.025 m: ~40 min; RSS ≈ 230 MB) |
| `python -m ice_integration.diagnosis.station_diagnosis OUT` / `clearance_check OUT` | Station table, training-only identifiability, independent-depth check |
| `python -m ice_integration.diagnosis.toe_convergence settled\|history\|windows …` | Settled-state and identical-window convergence, wet/dry sensitivity (private solver copy per threshold) |
| `python -m ice_integration.release.error_budget OUTDIR` | Signed budget numbers + `calibration_influence.png` |
| `python -m ice_integration.release.validity_contract META.json` | Contract check; exit 2 if any item unmet |
| `python -m ice_integration.rochemelon.check_readiness INPUT.json` | Readiness gate; exit 2 while gates are open (template fails all) |
| `python -m ice_integration.make_figures` | Three evaluation/demo figures (labelled as non-CI) |

## Layout
`forcing.py` (support accounting), `hydraulics.py` (pressure head vs depth, derived section quantities), `finite_heat.py` (depletion, ledgers), `coupling.py` (coupled driver, restart), `point_model.py` (V5 equivalent), `score_integrated.py`, `demo_coupled.py`, `diagnosis/`, `rochemelon/`, `release/`. Imports (unmodified): `ice_swiss_field_validation/`, `ice_external_forecast/` (pin 67f9b6ee…). Specification committed before scoring: `INTEGRATED_EVALUATION_SPEC.md`.

## Documents
`CLAIM_EVIDENCE_TABLE.md` · `ERROR_BUDGET.md` · `MODEL_VALIDITY_CONTRACT.md` · `DECISION_MEMO.md` · `../../INTEGRATION_REPORT.md` · `../../FOLLOWUP_DIAGNOSIS_REPORT.md`.

## Scientific claims that remain UNSUPPORTED
1. That the integrated model out-predicts constant-rate, constant-h or standard correlation baselines (it does not, consistently).
2. Any strict-native scored result (no interval has full native support).
3. A P5-specific parameter, or any parameter beyond one weakly identified scale.
4. That pressure-derived stage equals full depth; any sensor-clearance value.
5. Environmental prediction of incision, retreat, melt volume or lifetime at any real site.
6. Natural boulder/porous behaviour (drag closure unverified).
7. Field validation by the coupled demonstration (imposed Q, assumed strip, no target).
8. Confidence intervals: all ranges are engineering sensitivities or leave-one-out descriptions.
9. Rochemelon results (gates unresolved; no inputs invented).
10. Originality of any mechanism claim (see `DECISION_MEMO.md`).

## Known limitations of this release
Hydraulics in the original demo runs were capped at 60 s per interval (9/14 early intervals not fully settled; effect ≤4.4e-4 in f_toe). Contact timing at dx 0.05 is ≈5 h late and ≈9–13 % high versus 0.025. Engineering scenarios recalibrate per scenario, so the width effect is absorbed and is not an independent control.

