# Fixed-parameter numerical refinement

Computed from the published frozen normalized-heat predictions, without camera
outcomes, gain fitting, diffusivity refitting or changes to the study criteria.
Executable: `model/refinement.py`. Record:
`model/frozen_piramide_v1/NUMERICAL_REFINEMENT.json`. Two synthetic ratio/unit tests pass.

The frozen primary forward-fit diffusivity is`2.365417846e-6 m²/s`; its original
evaluation/training heat ratio is`0.7271392208`. The source-geometry primary
derivative-OLS comparator uses`6.654687215e-7 m²/s`, ratio`0.7558948594`.

| Method | Cells | Step s | Evaluation/training heat ratio | Change versus original ratio | Negative basal-flux substeps |
|---|---:|---:|---:|---:|---:|
| Forward fit |32|150|0.7271051521|−0.004685%|637|
| Forward fit |64|150|0.7270978099|−0.005695%|637|
| Forward fit |64|75|0.7270806142|−0.008060%|1292|
| Derivative OLS |32|150|0.7559094197|+0.001926%|0|
| Derivative OLS |64|150|0.7558933590|−0.0001985%|0|
| Derivative OLS |64|75|0.7559007462|+0.0007788%|0|

The ratio changes are also the predicted-ablation multipliers when applying the
same one-TRAIN-target gain-calibration rule at each resolution. If the original
numerical gain itself is held fixed, the finest-grid prediction multipliers are
`0.9999998813` for forward fit and`0.9999746705` for derivative OLS. Both
interpretations are recorded; no camera target is required to derive them.
Relative normalized-energy ledger residuals are at most`3.83e-13`. These finite
mesh/time differences are sensitivity evidence, not rigorous error bounds or
statistical uncertainty intervals.

**Physical gate:** the primary forward-fit prediction fails the frozen
no-refreezing melt-conversion condition. Its original32cell/300s record already
has314negative basal-flux substeps and50negative half-hour heat intervals.
Their negative integral magnitude is`0.588874 K m`, about0.523% of the net28-day
heat, with a minimum interval`−0.0256604 K m`. Refinement preserves the negative
flux; it is not explained away as roundoff. Stable numerical heat ratios do not
establish that a signed heat withdrawal equals observed ice ablation without
the missing interface-water/refreezing physics. The existing campaign criterion
must remain intact; selecting the comparator or deleting the cooling intervals
would change the declared primary experiment.
