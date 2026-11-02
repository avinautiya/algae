# Signed error budget (overprediction)

Sign convention: **+ = increases predicted incision (worsens the observed over-prediction)**; − = decreases it. Observed pooled bias of the frozen candidate on the six test intervals: **+0.63 m (nominal), +0.52…+0.82 m over all 81 scenarios, positive in 100 % of them.** The budget does **not** close: no item below is shown to account for the bias, and none is tuned to do so. Numbers: `ERROR_BUDGET_NUMBERS.json`.

Evidence grades: **D** demonstrated from data or computation in this package; **C** computed sensitivity with a hypothetical input change (sign/magnitude of the change is real arithmetic, whether the input is actually wrong is not shown); **H** hypothesis with no quantification here.

| Item | Sign on prediction | Magnitude (if any) | Grade | Basis |
|---|---|---|---|---|
| Calibration scale dominated by P3:3 (85 % weight); LOO range 3.6–9.6 vs 4.40 | + at test (P3:3 over-predicted +0.51 m at train, P3:2 under-predicted −1.21 m) | scale uncertainty from interval choice −18 %/+118 % | **D** | `calibration_influence.png`; the fitted scale is not representative of the other training intervals |
| Gap imputer (linear/left/right) | ± | mean pooled bias 0.60–0.66 m; ≤0.06 m shift | **D** | scenario table: cannot close bias |
| Clock shift ±2 h | ± | mean bias 0.56 (+2 h) … 0.75 m (−2 h): ≤0.19 m | **D** | cannot close bias |
| Temperature offset ±0.05 °C | ± | 0.61–0.65 m: ≤0.05 m | **D** | cannot close bias |
| Channel width 1.5–2.5 m | ≈0 after per-scenario recalibration | ≤0.005 m | **D** | width absorbed by the refitted scale; not an independent control |
| Temperature observed vs reconstructed | ± | reconstructed share ≤4.1 % on every test interval | **D** (small) | cannot explain 0.15–1.6 m errors |
| Pressure stage shallower than independent depth (17–27 Jul, P5: −0.02…−0.58 m) | + | heat-transfer coefficient +18 % per −0.1 m, −13 % per +0.1 m, −30 % per +0.3 m at Re≈9e5 | **C** | `depth_proxy_sensitivity…`; independent depth is itself indicative (edge cut, station label ambiguous, interpolated floor); sign flips on 28–30 Jul (proxy deeper by 0.10–0.17 m) so net effect on P5:3 unknown |
| Sensor above bed (clearance unknown) | + if sensor high ⇒ depth under-estimated | not quantifiable | **H** | notebook flags "not on bottom"; no numeric clearance |
| Floor-point temperature warmer than bulk/water-ice boundary layer | + | not quantifiable | **H** | floor sensor is not a bulk thermometer |
| Pipe-flow Nusselt scaling applied to rough open channel at Re 6e5–9.5e5 | + or − | not quantifiable | **H** | outside validated range; baselines Gnielinski/DB0.3 (lower Nu) do better at P5 |
| Heat consumed by cold ice (ice below 0 °C) | − (reduces melt) | latent cost + c_i(−T_ice); unknown T_ice | **H** | would reduce, not worsen, overprediction; ice temperature not measured |
| Assumed Q at P5 equals Q at P3 (melt/lateral inflow ignored) | + at P5 if true Q is lower upstream | not quantifiable | **H** | P5 over-prediction larger than P3 is consistent but not tested |
| Daily-mean depth inside nonlinear closure | ± | small (h ∝ depth^−0.2) | **H** | Jensen effect, not computed |
| Floor survey ±0.1–0.3 m per epoch, interval sums 0.2–0.5 m | ± random | up to ±0.5 m per interval | **D** | explains scatter of short intervals only |
| Lowering measured as floor elevation: surrounding-ice surface ablation, channel migration | − (observed lowering larger than melt-only) | unknown | **H** | omitted process that would increase observed lowering, i.e. act opposite to the bias |
| Viscous dissipation (omitted in point predictor) | − | small relative to 1 °C·h forcing; not computed | **H** | would add melt (opposite sign) |
| Observation-location mismatch (sensors ~2 m downstream of station) | ± | not quantifiable | **H** | source note |

Reading: demonstrated sources (imputer, clock, T, width, reconstruction share, survey σ) are all too small to remove a +0.5…+0.8 m bias; the demonstrated calibration-leverage problem explains why a scale fitted mainly to one interval transfers poorly; the remaining candidates (clearance/depth proxy, floor-vs-bulk temperature, Nusselt closure, Q conservation) are hypotheses that **cannot be tested with the present public data**. Residual unexplained bias is reported as unexplained.
