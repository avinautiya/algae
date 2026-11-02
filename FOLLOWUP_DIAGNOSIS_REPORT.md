# Follow-up: figure audit, station diagnosis, identifiability, geometry sources, toe-contact convergence, go/no-go

**The failed accuracy verdict stands unchanged** (see `INTEGRATION_REPORT.md`): the integrated candidate does not beat every baseline consistently. Nothing below re-scores, re-selects or re-tunes anything. Existing demo runs, scores and figures' data are retained; only figure labelling was edited.

## 1. Figure audit (all three figures opened and inspected)
| Figure | Defects found | Fix (labels only; data unchanged) |
|---|---|---|
| `nominal_mae.png` | Axis said MAE but no x-unit label, no n, no statement that there is no interval; P5 panel on a different x-scale unmarked | x-label "mean absolute error of incision (m)"; panel titles carry n (4/2/6); footnote: single nominal scenario, no confidence intervals, retrospective inspected data |
| `scenario_wins.png` | Bars could be read as probabilities/intervals; zero-win bars are invisible so P5 looked "missing"; no statement of what a scenario is | title states these are COUNTS over engineering scenarios, not probabilities or confidence intervals, scenarios are not independent draws, zero-height bar = 0 wins; axis "scenarios won (of 81)" |
| `demo_dynamics.png` | Overlapping lines hid runs; x-axis unlabelled; title implied settled hydraulics; no "not validation" statement | strict solid / reconstructed dashed; x-label with interval length and naive-clock caveat; footnote: not field validation and hydraulics were capped at 60 s per interval |
Scenario ranges (clock ±2 h, T ±0.05 °C, width 1.5/2/2.5 m, three gap imputers) are engineering sensitivity envelopes. They are not confidence intervals; the only uncertainty with a statistical meaning in the data is the source's per-epoch standard-deviation column (0.1–0.3 m), reported as `obs_sigma_m` (sum of two epochs, a scenario value, not a CI).

## 2. Station-level error diagnosis
Source: `ice_integration/results/STATION_DIAGNOSIS.json` (`python -m ice_integration.diagnosis.station_diagnosis`). Nominal scenario, frozen calibration (scale 4.396 recomputed from training rows only; matches the frozen value exactly).

| id | role | days | obs (m) | pred (m) | err (m) | temp reconstructed share | engineering pred range (m) | obs inside range ± σ? |
|---|---|---|---|---|---|---|---|---|
| P3:2 | train | 5.1 | 2.19 | 0.98 | −1.21 | 1.9 % | 0.76–1.13 | no |
| P3:3 | train | 8.1 | 2.05 | 2.56 | +0.51 | 12.8 % | 2.53–2.57 | no |
| P3:4 | train | 0.8 | 0.73 | 0.40 | −0.33 | 1.9 % | 0.37–0.47 | no |
| P3:5 | test | 5.1 | 1.54 | 2.14 | +0.60 | 3.4 % | 2.01–2.35 | no |
| P3:6 | test | 0.8 | 0.07 | 0.22 | +0.15 | 4.1 % | 0.18–0.24 | yes |
| P3:7 | test | 8.1 | 0.66 | 1.13 | +0.47 | 0.7 % | 0.60–1.46 | yes |
| P3:8 | test | 6.0 | 1.16 | 1.43 | +0.27 | 1.6 % | 1.16–1.60 | yes |
| P5:0 | train | 0.8 | 0.13 | 0.01 | −0.12 | 1.9 % | 0.00–0.02 | yes |
| P5:3 | test | 5.2 | 1.86 | 3.49 | +1.63 | 4.1 % | 2.98–4.39 | no |
| P5:4 | test | 0.8 | −0.01 | 0.37 | +0.38 | 2.6 % | 0.32–0.47 | no |
(P3:0/P3:1 are tiny training intervals, omitted for space; both inside range.)
Source of uncertainty per station (descriptive; test rows are inspected data, so this explains errors, it does not select anything):
- **Observed inputs.** Temperature is ≥87 % natively observed in every interval (reconstructed share ≤13 %, ≤4.1 % on every test interval), so gap bridging cannot explain a 0.4–1.6 m error. A floor-point sensor is not a bulk-water thermometer (labelled, not quantified). Survey error is 0.2–0.5 m summed; it explains the P3:6, P5:0 and small-interval scatter only.
- **Reconstructed/proxy inputs.** Discharge and P3 depth are author reconstructions (lake-balance rating, bridging 13–24 July pressure gap); P5 depth is a pressure proxy; Q is assumed conserved P5→P3; width is an assumption (2 ± 0.5 m). Clock ±2 h, T ±0.05 °C, width 1.5–2.5 m, imputers move predictions by 10–40 %, which is smaller than the error at P3:2, P3:3, P3:4, P3:5, P5:3 and P5:4 (observation outside the range ± σ, error sign identical in all 81 scenarios in all of these except P3:7/P3:0). Perturbing the labelled inputs therefore does not remove these errors.
- **Calibration assumptions.** One scale for both stations. Training residuals alone have opposite signs at one station (P3:2 −1.21 m, P3:3 +0.51 m), so no single multiplicative scale is consistent even within the training data. P3 test errors are all positive (+0.15…+0.60 m) and P5 test errors are positive (+0.38, +1.63 m): a systematic over-prediction at test time, i.e. the scale fitted mostly to P3:3 transfers badly.
- **Structural approximations.** Pipe-flow Nusselt scaling (Re 6e5–9.5e5, outside the validated pipe range and for a rough open channel), temperate ice at a point, imposed Q without meltwater feedback, floor-point temperature as bulk driver, no atmospheric exchange, daily-mean hydraulics inside the nonlinear closure. The sign pattern (over-prediction at P5 where the constant-h baseline fails worst) is consistent with the heat-transfer scale being fixed by P3 conditions, but this is a hypothesis; it is not tested by any fit here.

## 3. Parameter identifiability (training data only)
- One calibration scalar, six training intervals. Share of Σb² (basis²): **P3:3 85 %, P3:2 12 %, P3:4 2 %**, the other three ≈ 0.003 % in total. Effective sample size ≈ 1.3 intervals.
- Leave-one-out scale: **3.62 – 9.58** versus 4.40 (−18 % to +118 %). Propagated survey sigma on the scale alone: ±0.32 (7 %), so data-noise is not the dominant uncertainty; interval choice is.
- Implied per-interval scale (obs/basis): P3:0 5.6, P3:1 −20.1, P3:2 9.9, P3:3 3.5, P3:4 7.9 (P3:1 negative because obs < 0). P5:0 52.6 ± ~30 (observation 0.13 m vs σ 0.2 m, 0.76 days).
- Conclusion: the scale is **weakly identified for P3 and not identifiable for P5** from training data. P5 has one 0.76-day training interval with a signal-to-noise ratio < 1. No P5-specific parameter was fitted to P5 evaluation outcomes, and this report does not recommend one; evaluation outcomes must not be used for that.
- Eight exclusions in the nominal scenario (P3:9, P3:10, P5:1, P5:2, P5:5, P5:6, P5:7, P5:8) lose most of the long, high-signal P5 intervals because depth/discharge is missing or non-positive; this is why P5 is constrained so poorly. Exclusions are availability, not outcome-based.

## 4. Do independent geometry and sensor-clearance records exist in the original public sources?
Checked: author-hosted ZIP mirror (CRC-verified members), `source/` notebook, `original_mirror/data/Channel_elevation/*.csv`, products, OBSERVATION_OPERATOR, HEAD_DROP audit.
| Item | Exists? | Detail |
|---|---|---|
| Surveyed channel-floor elevation | **Yes** | P1–P5, 6–12 dates each, ±0.1 m (±0.2/0.3 m for late epochs), UTC header (unreconciled with CEST plots) |
| Station along-channel distances | Yes | P5→P4 151.3 m, →P3 253.7 m, →P2 312.6 m |
| Water-surface (edge) elevation | Partly | `P5_levels_elevation.csv`: 15 daily values 16–30 July, ±0.05 m; **header says "P4 water level cuts"** and times are "hypothesis, only 30 July true" |
| Sensor elevation / clearance above the bed (z_s − z_b) | **No** | Only qualitative notebook notes ("sensor located at the bottom", cells 79/91 "not on bottom"/"in the air"); no numeric clearance record |
| Cross-section width/shape | **No (measured)** | Notebook uses assumed width with 0.5 m uncertainty; no surveyed cross-sections |
| Direct discharge at P5/P3 | Partial | Tracer-derived velocities; channel discharge is an author lake-balance reconstruction |
Descriptive cross-check (`ice_integration/results/CLEARANCE_CHECK.json`, no model, no fit): independent depth η − z_b (water-level cut minus time-interpolated floor; ±0.15 m) versus author pressure-derived stage on 13 days. The pressure proxy is **0.02–0.58 m shallower** during 17–27 July and **0.10–0.17 m deeper** 28–30 July: differences of 2–4σ. This shows the stage proxy does not represent full depth consistently but **does not determine clearance**: the level file's station label is ambiguous (P4 vs P5), the floor is interpolated across 5–9-day gaps while the floor drops ~0.4 m/day, and the level cut is at the channel edge. **Conclusion: independent geometry exists only coarsely and at sparse times; sensor clearance is not available in the public sources; only a numeric correction from original logger metadata or the authors could supply it.** No clearance was inferred or assumed.

## 5. Toe-contact convergence (existing runs retained; new runs in `ice_integration/results/convergence/`)
Definition: f_toe = Q_toe / Q at the mid-section. Scripts: `ice_integration/diagnosis/toe_convergence.py`.
- **Problem in the earlier statement.** The earlier 12 % difference averaged f_toe over all 56 intervals, 42 of which are zero: it mixes contact and non-contact periods. In addition the demo capped hydraulics at 60 s per interval; 9 of the first 14 settle only after 60–116 s, i.e. the existing runs were not fully settled in those intervals (water/energy residuals were still ≤1e-12). The effect is small (existing vs fully settled coarse run: max |Δf_toe| = 4.4e-4) but is recorded here, and the old numbers are not deleted.
- **Identical bed, fully settled state (k=0, T=0.312 °C, Q=0.07 m³/s, tol 2e-4).**
| dx (m) | f_toe | abs diff to next finer | rel diff to finest |
|---|---|---|---|
| 0.05 | 0.12135 | 0.00807 | 8.7 % |
| 0.025 | 0.11328 | 0.00164 | 1.5 % |
| 0.0125 | 0.11165 | – | – |
Ratio 4.9 ⇒ observed order ≈ 2.3; Richardson-extrapolated limit ≈ 0.1112 (indicative only, 3 grids).
- **Identical windows, fully settled coupled history (hyd cap 3000 s; dx 0.05 vs 0.025, first 14 intervals, reconstructed forcing).** Mean f_toe, k0–5: 0.0949 vs 0.0867 (abs 0.0082, rel 9.4 %); k0–13: 0.0566 vs 0.0499 (abs 0.0067, rel 13.4 %). Per-interval relative difference grows from 7 % (k=0) to 33 % (k=10) and is undefined/large when f_toe → 0 (k≥12), so relative differences there are not meaningful. Contact loss time (f_toe < 1e-3): 82.2 h (coarse) vs 77.1 h (fine), Δ = 5.2 h (6.7 %). Settled cumulative melt: toe 0.1468 vs 0.1447 m³ (1.5 %), bypass 0.1928 vs 0.1920 m³ (0.4 %). Existing 56-interval melt totals agree to ≤0.6 %.
- **Wet/dry threshold sensitivity (HDRY 1e-5, 1e-4, 1e-3 m).** Settled f_toe at dx 0.05 and 0.025 is identical to 12 digits; the 14-interval coupled history at dx 0.05 is bit-identical (max |Δ| = 0, melt identical). Each value was run on its own private copy of the solver (the constant is baked into cached numba kernels, so an in-place change would be stale). In this strip geometry the threshold is therefore not controlling; the test does not cover geometries with extensive thin-film flow.
- **Verdict.** Toe-contact fraction is converging (observed order ≈ 2) but **coarse-grid 0.05 m values are ≈ 9–13 % high and the contact loss time ≈ 5 h late**; medium/coarse results are sufficient for melt volume (≤1.5 %) but not for contact timing at better than ~10 %. Wet/dry threshold is not a source of error here.

## 6. Go / no-go: what could support a genuinely independent test
| Candidate observation | Public status (verified or not) | Closes which gap | Decision |
|---|---|---|---|
| Rochemelon 2005 graph extraction (floor profiles A–C, 3 survey dates, lake T, lake-level discharge) | Figures only; no numeric files found; heat-transfer coefficients in paper were tuned to this event | Distinct event, would be untouched by current Swiss development | **NO-GO now.** Gate inputs unresolved (section binding, survey hours, inlet-to-target path, geometry); GO only after audited digitization with uncertainty and with no use of the paper's tuned coefficients; templates in `ice_integration/rochemelon/` |
| Plaine Morte original logger metadata for sensor elevation/clearance, hourly width/stage | Not in public archive; would need authors | Clearance, cross-section | **NO-GO** until obtained; no inference allowed |
| Plaine Morte 2019 later events beyond the exposed intervals | Same season already inspected | Not independent | NO-GO as a "new blind" test; usable only as extra development data |
| Pohle et al. 2022 Rhonegletscher moulin deposit | Citation verified, archive listing not verified (HTTP 429), pressurised vertical channels | Different physics | **NO-GO** for this open-channel model; separate applicability audit first |
| Karlstrom 2014 Llewellyn | Temperature/discharge but no channel-bed incision series | Missing target | NO-GO |
| Isenko & Mavlyudov 2002 Fisht | Four-day mean incision, widths, slope; no paired T/Q series | Plausibility only | NO-GO as test; usable as order-of-magnitude check |
| Plaine Morte water-level cut file corrected for station/time | Header says P4, time hypothesis | Independent depth at ≤3 dates | **CONDITIONAL**: GO only if the station and times are confirmed from the notebook/authors; would validate depth proxy, not heat transfer |
| New field survey with floor, water surface, clearance, width at sensor stations plus bulk-water temperature and discharge | Not public | All | The only unconditional GO (future data) |
Nothing here is field validation; the demonstration runs are not field validation; no Rochemelon input was invented.
