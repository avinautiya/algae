# Release README — ice_integration (milestone M4 + follow-up)

**Status: verified computational integration and a retrospective, non-independent field comparison. The integrated model did NOT improve on every baseline consistently. Nothing here is independent validation.** All frozen results, the failed-superiority verdict, the original demo runs and superseded versions are preserved. No recalibration against evaluation intervals and no new production campaign were run for this package.

## What runs (verified in this checkout)
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
