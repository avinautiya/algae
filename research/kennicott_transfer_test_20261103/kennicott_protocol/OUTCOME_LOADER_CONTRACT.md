# Native outcome release contract

`outcome_loader.py` has no import-time source access. No native outcome call has been made during implementation or its ten synthetic/mock tests.

Before publication, put `release_plan()` into the complete experiment spec under `native_outcome_release`. This binds the exact sourceCSV SHA256/CRC, operator and four native site/date/geometry rows. Publish the spec and normalized predictions before authorizing TRAIN. Calibrate the declared models from the one2020 row, publish all final calibrated2023 predictions, then authorize EVALUATION. A source mismatch, missing/duplicate row, changed dates, geometry, nonfinite outcome or changed snow/surface condition fails rather than choosing another row or operator.

```
from kennicott_protocol.outcome_loader import release_plan, load_after_freeze
# spec['native_outcome_release'] = release_plan() before freeze/publication.
# Only after the root explicitly authorizes the published receipt:
released = load_after_freeze('receipt.json', expected_role='train')
```

A receipt contains:

```
{
  "schema": "kennicott-outcome-release-v1",
  "access_role": "train",
  "freeze_published": true,
  "predictions_published": true,
  "publication_reference": "actual published commit or artifact reference",
  "prediction_stage": "normalized_training",
  "frozen_spec_path": "relative/path/spec.json",
  "frozen_spec_sha256": "exact file hash",
  "predictions_path": "relative/path/normalized_predictions.json",
  "predictions_sha256": "exact file hash",
  "approved_operator": "negative_source_mass_change"
}
```

For evaluation, use `access_role=evaluation`, `prediction_stage=final_evaluation`, and the separately published final prediction artifact. Relative paths resolve from the receipt's directory. TRAIN returns only P_AWS_1 July25–September6 2020; EVALUATION returns only P11/P12/P13 June8–September2 2023. Returned values are negative source mass change inkg/m², already equivalent to mm water equivalent. No900/917 multiplier is applied, no averaging by duration is performed, and negative ablation/net accumulation remains negative. This is a conditional net-loss comparison, not independently verified gross phase-change melt.

These are auditable logical role gates on a public shared CSV, not cryptographic blinding or verification that a publication reference genuinely exists. The root verifies the actual publication before calling. Fetching the native CSV necessarily retrieves public bytes containing multiple rows, but unauthorized role mass cells are not numerically converted or returned. No raw complete outcome CSV is saved by the loader.

From packet root: `python3 -m unittest kennicott_protocol.test_outcome_loader -v`. Tests verify role isolation with deliberately nonnumeric unauthorized targets, signed targets, surface/snow/geometry gates, duplicate rows, publication-before-fetch, stage mismatch, hash tampering and source CRC/SHA. They do not access real outcomes.
