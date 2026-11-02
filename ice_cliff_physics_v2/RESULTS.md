# Results of the cliff-physics upgrade

## What was actually implemented

The reduced prediction CLI now couples finite stream heat to bare ice and separate wet/dry debris columns, optionally resolving normal debris diffusion and height-dependent face retreat. Material transfers conserve enthalpy. True profile volume is distinct from the basal routing proxy. Water/cover boundaries are inserted losslessly. Explicit local transfer validity gates, route-specific meltwater checks and zero-inflow drainage contracts prevent unsupported calculations from looking successful.

The production core passes 131 isolated tests. Its original cumulative totals reproduce the previous implementation exactly. Tests establish numerical properties, not environmental accuracy.

## Pre-execution controlled comparison

Code, cases and source hashes were published before execution at `b9d6f39d09f7ead306ea87c852da8f14c53e9e89`. All 30 coupled cases and 12 fixed-boundary comparisons completed; no failed or refused case was dropped. The maximum coupled time-series normalized heat-ledger residual is 2.26e-15; maximum geometry–phase mismatch is 4.30e-15 m³. `results/ANALYSIS.json` and `reviews/FINAL_RESULT_AUDIT.md` independently reproduce the contrasts.

On an identical fixed-elevation covered footprint, the single-temperature approximation overpredicts the refined model's ten-diffusion-time melt by about 1.08%. FV16 versus FV32 differs by about 0.0023%. Disabling water–debris heat exchange changes covered-footprint melt by −49% or +78% relative to the coupled model. The sign depends on whether the finite stream heats or cools the debris. The disabled-exchange calculation is an intervention, not an implementation of a professional glacier model. Equal prescribed net heating above and below water is a controlled absorbed source, not identical weather. These are not field-accuracy percentages.

Height-resolved bare-face calculations produce approximately 5.06 cm and 10.68 cm maximum geometric undercut after 48 h under the supplied 0.5 K and 1 K water boundaries. Height refinement from 32 to 64 initial bands agrees after adaptive class alignment. This is a geometric prediction conditional on adiabatic horizontal ledges; no failure, collapse, cliff-lifetime or natural-boulder claim follows. Whole evolving-floor outputs do not establish fixed-footprint floor accuracy.

## Independent analytic reference exposes an early-time numerical limit

At one diffusion time, the single-temperature closure's cumulative ice-heat error spans −30.27% to +13.17%; the refined FV comparison is within 0.0541% of the classical Robin solution. At ten diffusion times FV error is below 0.00019%.

**At 0.1 diffusion times, however, FV64 with 400 steps per diffusion time still has 11.78–16.41% relative error.** Temporal refinement dominates the remaining spatial difference. Conservation alone does not certify early-time accuracy. A separate 48-run post-study audit (all 12 cases, 64/128 cells, 4,000/16,000 steps per diffusion time) reduces the early error to 0.331–0.469% at the finest tested settings. Both time and space matter. The complete original results and the audit are retained; no frozen result is replaced by a better-looking setting. Relative errors at vanishing heat/cancellation scales are not evidence of physical significance.

## Public field tests: real gains and real losses

The separate Kennicott temperature test uses measured shallowest/deepest boundary probes and predicts 27 interior series across 11 retained profiles. It trains on 13 days after a one-day initialization and scores the next 14 days, with no interior-state reset. It is a conditional thermal hindcast, not a forecast of future weather or water contact.

| Model | Equal-profile mean absolute temperature error (°C) |
|---|---:|
| Instantaneous quasi-steady interpolation | 1.1933 |
| Homogeneous diffusion | 0.7552 |
| Causal lag-filter baseline | 0.7534 |
| Training-only boundary ridge regression | 0.2863 |

Diffusion improves on quasi-steady interpolation at all 11 profiles, roughly 36.7% in aggregate. It does not outperform the strongest baseline. Profiles from one glacier are not independent glacier replications. A constructed zero-temperature boundary in the earlier 2020 profile is excluded from the primary measured-endpoint temperature test.

The separate ablation transfer benchmark retains all 1,728 edge/model/metadata combinations. For the 2020-source to 2023-target transfer, distributed conduction gives approximately 21.1–23.0% mean absolute percentage error depending on predeclared thickness choice; source-rate transfer gives 17.5%. Across all 12 directed transfers, conduction gives about 29.3–29.6%, versus 24.0% for source-rate transfer. Quasi-steady conduction is essentially tied with diffusion on seasonal totals. The 64→128-mode ablation difference is at most 0.00334 kg m⁻².

2023 thermal records miss the first two days of the stake interval. Full-period estimates therefore include explicitly conditional leading-day reconstructions and date/depth scenarios, not complete observed forcing. Conductivity is calibrated on each source origin's ablation only; temperatures identify diffusivity but not conductivity and heat capacity separately. These tests do **not** establish improved field ablation, coupled cliff-melt accuracy, Himalayan transfer or professional-model superiority.

## Prospective vertical-heterogeneity test

A separate protocol, verified two-layer solver and reserved scoring window were published at `4549f6c92b8ea0bca32b58d70f5d28f3007f994c` before evaluating days28–42. Training remains the same early13-day interval. The fixed midpoint/common-capacity two-layer model searches a predeclared9×9 diffusivity grid; all81training objectives are retained. Nine profiles,23interior series and7,751 reserved predictions satisfy the availability criteria.

Two-layer mean absolute temperature error is0.5217°C versus0.7477°C for homogeneous diffusion,1.1823°C quasi-steady,0.7664°C lag-filter and0.2314°C ridge. The physical-model improvement is30.2% in aggregate; it wins6/9homogeneous comparisons,9/9quasi-steady and0/9ridge. **The frozen all-baseline superiority gate fails.** This is a bounded physical-model improvement, not universal superiority. The40→80-cell difference is at most0.00515°C and an independent raw-value/metric audit passes. Neither conductivity, melt flux nor the unique physical cause of the observed profile curvature is established by temperature fitting. The two-layer model is an effective hypothesis, not measured layer composition or a new heat-conduction law.

## Scientific contribution and remaining decision

The deliverable is a more defensible, executable and reproducible coupled model, an audited analytic/numerical benchmark, and an actual public-field comparison that diagnoses omitted physics. The general heat-storage, Robin-conduction and wet thermal-resistance mechanisms are established physics. Historical originality of the precise application remains a literature question; a groundbreaking natural discovery has not been demonstrated by these runs.

The strong field baseline's advantage is scientifically useful: it shows a homogeneous conductive column cannot be assumed adequate merely because its solver is accurate. The separately frozen vertical-heterogeneity test improves on homogeneous conduction but still loses to ridge; tuning against either now-inspected test is not an independent confirmation. Coupled cliff field accuracy still requires observations matching water forcing/contact and the same moving ice footprint. More parameter combinations cannot supply missing observations.
