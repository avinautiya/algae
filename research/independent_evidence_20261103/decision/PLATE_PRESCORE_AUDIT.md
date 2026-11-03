# Independent pre-score audit — conditional plate-signal experiment

Scope: source review of `extension/plate_experiment.py`, `dirichlet_slab.py` and `PROTOCOL.md`, plus the seven synthetic/unit tests. I did **not** execute held-out scoring or inspect held-out target values, fitted model results or prediction errors.

## Must-fix found and resolved before freeze

The original `prepare` accepted time spacings shorter than one hour, while `lagged` interpreted row offsets as1/3/6/etc hours. This could change the declared temporal features on irregular data. The owner added `validate_hourly` requiring every consecutive difference to be3600s (absolute tolerance1e-5s), and a regression test. That resolves the source contract for these row-based lags.

## Checks passed by source inspection

* Training and evaluation masks use explicit Nepal UTC+05:45 boundaries and half-open intervals. Each model uses the same finite-plate target timestamps for evaluation; `scores` refuses missing/nonfinite predictions rather than dropping difficult target rows.
* Site1alpha uses only training interior-temperature targets. The initial profile uses the first source measurement before training and is never reset at the split. Held-out interior temperatures are not supplied to the evolving model or ridge features.
* `lagged` internally selects only the first/last endpoint temperatures, even when a three-probe array is passed. The suspected held-out interior-feature leak is **not present**.
* Every flux-amplitude fit uses only finite training plate readings. Sign is inferred on training only and retained; no test-driven orientation reversal or clipping is present. The synthetic held-out-flux perturbation test leaves fits/predictions unchanged.
* Ridge features are causal in the representative-timestamp convention: current observed endpoints and previous endpoint rows. Its penalty is selected within training and refitted on all training only. The training median and instantaneous-gradient baseline have their own training fits.
* Site2 uses site1alpha, but fits its **own** local amplitude on site2training plate readings. It is an alpha-transfer/local-calibration time-holdout test, not an unseen-site or zero-shot parameter-transfer test.
* The source metadata declares plate readings in W/m², and the loader preserves those values. No latent-heat or melt conversion appears in this pipeline.
* Default execution refuses native interval-flux interpretation. The explicit `representative_point_diagnostic` option permits only the documented conditional raw-signal comparison. This does not solve unknown averaging convention, plate orientation/position or sensor migration.
* The numerical check compares64/128modes using the same fitted coefficient. Its tolerance depends on training plate RMS, not the held-out error. A failed gate remains visible; the code must not be changed after scoring just to pass it.

## Qualifications to retain in the final report

1. Site1diffusion receives an auxiliary training interior-temperature series to estimate alpha; the baseline regressions use endpoint temperatures and training plate targets. These are legitimate stated algorithms, but a win cannot isolate the causal value of diffusion physics from the value of that auxiliary training information without an additional ablation.
2. Current measured endpoint temperatures are prediction inputs. This is conditional prediction of a colocated raw signal, not forecasting future weather, unseen boundary temperatures, latent melt, glacier retreat or professional-model superiority.
3. The signed coefficient combines scale, orientation and effective representation error. Calling its absolute value an independently validated intrinsic conductivity would be wrong. Temperature identifies diffusivity; plate-amplitude fitting adds information, under this assumed observation operator.
4. Two nearby sites on one glacier and correlated hourly samples do not provide independent glacier replications or a confidence interval.
5. Native logger averaging/sign/placement remains unresolved. Representative points must not be relabelled as exact interval-matched physical heat flux after results arrive.
6. A64/128fixed-parameter signal check does not establish unique parameter identification, phase-change validity or whole-model physical accuracy. Endpoint boundary data may include freezing temperatures; the linear sensible subcolumn is a diagnostic, not an enthalpy model.

## Disposition

**No remaining must-fix source bug identified for the narrowly declared representative-point diagnostic**, after the hourly-spacing repair. This is a bounded review, not a proof that every software/scientific error is absent. Verify source/data hashes and publish the freeze before any held-out score; preserve all outcomes and numerical failures. Historical originality and melt accuracy are not established by this experiment.

Seven reviewed tests pass. No hold-out target was scored by this audit.
