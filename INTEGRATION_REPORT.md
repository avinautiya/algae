# Integration report (M4)

Branch `claude/integration-m4` (base dc4a696). Imported unmodified: `ice_swiss_field_validation/`, `ice_external_forecast/` (pin 67f9b6ee…). No solver, frozen spec or historical result was altered. Spec: `INTEGRATED_EVALUATION_SPEC.md` (committed before scoring).

## Direct answer
**No.** The integrated model (candidate = scale-fitted Dittus-Boelter 0.4 with finite-heat accounting) does not improve the declared target against every baseline consistently.
Nominal test MAE (m, 6 test intervals; P3 n=4, P5 n=2), from `ice_integration/results/EVALUATION_RECORD.json`:

| baseline | P3 | P5 | pooled |
|---|---|---|---|
| constant rate (fitted) | 0.658 | 0.280 | 0.532 |
| constant h (fitted) | 0.330 | 1.188 | 0.616 |
| DB0.4 (unfitted) | 0.577 | 0.580 | 0.578 |
| DB0.3 (unfitted) | 0.642 | 0.662 | 0.648 |
| Gnielinski smooth (unfitted) | 0.454 | 0.420 | 0.443 |
| **candidate** | **0.374** | **1.007** | **0.585** |

Candidate wins in P3 against 4 of 5 baselines (not constant-h: 0/81 scenarios) and loses at P5 against all baselines except constant-h. Pooled it beats only constant-h in nominal. Over 81 engineering scenarios (clock ±2 h, T ±0.05 °C, width 1.5/2/2.5 m × 3 imputers) pooled wins vs constant-h 81, DB0.3 44, DB0.4 28, Gnielinski 6, constant-rate 4.
**Failing regime:** P5 (2 test intervals, fitted h inflated) and P3 versus constant-h. **Evidence needed for the next justified change:** independent section geometry / sensor clearance, a per-station calibration constraint at P5, and an independent event not yet inspected. Existing Swiss cases are inspected retrospective data, not a blind test.

## Implemented / verified
- `ice_integration/forcing.py`: offsets applied once; time-weighted integration; observed / reconstructed / unsupported seconds kept separate; strict-native view. Strict-native scoring: no complete manual interval has full native support, so **0 intervals** could be scored strictly.
- `hydraulics.py`: pressure head ≠ depth (clearance required); area, D_h, velocity, Re recomputed per section; basis labels; refusals (no silent P3→P5).
- `finite_heat.py`: expm1-stable depletion, incoming heat allocated once, bed/wall split by wetted perimeter, fresh-ice cost L + c_i(−T_ice); ledger residuals; restart state + config fingerprint.
- `coupling.py`: coupled driver + atomic checkpoints (refuses provenance/dx/mode mismatch).
- Review of `predict_external.py` vs existing coupling: same physics in plug-flow limit (≤1 % vs coupled solver); the coupled solver's h_eff·T sink is the same flux that lowers the floor. Imposed Q is not updated by meltwater.
- Tests: 15 integration tests pass; Codex's 39 Swiss + 13 external tests were run at import. Equivalence with Codex V5 predictions: max rel. diff 4.2e-16 (n=5994).
- Not claimed: moving geometry, cold-ice PDE, full atmospheric exchange.

## Coupled demonstration (P5 native forcing 2019-07-20→08-03, 56×6 h, Q=0.07 imposed)
| run | wall | RSS | melt toe / bypass (m³) |
|---|---|---|---|
| dx 0.05 strict | 301 s | 227 MB | 0.483 / 0.703 |
| dx 0.05 reconstructed | 310 s | 227 MB | 0.508 / 0.742 |
| dx 0.025 strict | 2374 s | 228 MB | 0.481 / 0.703 |
| dx 0.025 reconstructed | 2413 s | 228 MB | 0.505 / 0.740 |

- Water and energy residuals ≤ 1e-12 in all runs; support accounting identical across resolutions (observed 1,167,235 s; reconstructed 42,365 s or counted unsupported in strict mode).
- Melt volume converged (coarse vs fine ≤0.6 %). **Flow-contact (toe fraction) is not converged**: mean 0.0141 vs 0.0124, max interval difference 0.0086. Reconstruction adds ≈5 % melt. Checkpoint 0.42 MB. Cost scales ≈ dx⁻³ (≈8× per halving); full 56-interval case at dx 0.025 ≈ 40 min. No parameter search was run.
- This is a demonstration with a strip deposit, not field validation: no coupled cliff-model field validation exists.

## Rochemelon gate: unresolved (not scored)
`ice_integration/rochemelon/`: `INPUT_TEMPLATE.json` and `check_readiness.py` (exit 2 while any gate open; template currently fails all). Open: measured A/B/C ↔ prediction section binding (≈3/17/33 m), survey timing/uncertainty, inlet-to-target thermal path, independent geometry or qualified closure (no Rochemelon-fitted coefficients), alignment flags. No inputs were invented.

## Assumed / unsupported
Assumed: tube correlations for rough open channels (Re gate 1e4–5e6), naive source clock, gaps ≤24 h bridge. Unsupported: strict-native scoring, confidence intervals (ranges are engineering sensitivities), independent validation, any claim of coupled-model skill.

## Reproduce
`python -m pytest ice_integration`; `python -m ice_integration.score_integrated`; `python -m ice_integration.demo_coupled OUT.json DX N [strict|reconstructed]`; `python -m ice_integration.make_figures` (figures in `ice_integration/results/figures/`).
