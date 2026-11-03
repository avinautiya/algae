# V2 cold-interface post-release audit

Run `audit_actual_cold_results.py` **only after the parent explicitly announces evaluation release**. It never opens camera files. No actual calibration, final prediction or result artifact was read while preparing this code.

CLI arguments are `--evaluation-released --normalized PATH --predictions PATH --seal PATH --results PATH --spec PATH --decoded-root PATH --physical-freeze PATH --out NEW_PATH`. The physical freeze is the V2 model's complete frozen JSON; its hash must match the normalized result header. Existing outputs are never overwritten.

The checker retains all39models:36thermal rows from2methods×3geometries×2holds×3capacity scenarios, plus constant-rate, buried-temperature and reconstructed-air-temperature baselines. The primary remains `forward__middle_gap_20mm__left__proxy`. All thermal predictions must satisfy Ypred=Ytrain×H1/H0 on identical camera endpoints and the native917 conversion.

Here H is **positive latent phase-change energy/Cscene**, not signed sensible interface heat. The audit independently sums trace values and checks latent energy, raw phase mass, unscaled sensible-plus-latent ledger, capacity identity, nonnegative phase mass and complementarity. Negative sensible cooling is supported by the explicit ice state; it is not clipped into latent heat.

The gain×L diagnostic is explicitly not an inferred physical heat capacity. The empirical gain relative to the physical component is gain×L/Cscene. The audit verifies that applying it reproduces the issued forecast while keeping the physical ledger attached to the **unscaled** raw simulation. It records assumed initial ice temperature, finite ice depth, fixed geometry and omitted wet-debris phase heat.

Primary errors, biases and ties are recomputed independently. Eighteen matched scenario summaries are descriptive counts, not independent samples or alternative primary winners. Positive-temperature baseline doses are reintegrated using a separate sign-crossing algorithm; the18h permitted air-source gap stays explicit. Shared-middle measurement-error ranking thresholds and fixed-issued-forecast thresholds remain distinct, conditional diagnostics with unknown source precision, not confidence intervals.

Six synthetic tests pass: complete39-model reconstruction, release refusal before artifact access, corrupted result detection, corrupted raw latent-energy detection, independent zero-crossing integration and exhaustive shared-reading corners around the rank threshold. The original V1 auditor and failed V1 physics gate remain unchanged.
