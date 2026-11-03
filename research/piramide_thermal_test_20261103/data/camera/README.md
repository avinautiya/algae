# Sealed camera source and role-gated loader

Three camera MAT members were retrieved as one72,289-byte public archive range, individually decompressed and ZIP-CRC32 verified. `SEALED_SOURCE_MANIFEST.json` binds their source bytes. Complete archive MD5 is not asserted. Source: https://doi.org/10.5281/zenodo.11581073.

`audit_schema_only.py` inspects ONLY native schema/timestamps. It does not dereference timetable data cells or summarize outcome values. Timestamp-only lists permit endpoint choice by coverage. Native timezone is empty, matching the process record; UTC is not invented. This is programmatic outcome sealing, not an operating-system security boundary.

Native tracks columns are distance/ablation, variable units empty. Source README states distance in metres; its ablation-unit description conflicts with apparent scaling in the calibration script. The loader therefore never uses provided ablation. Finalpaper section2.2.3 (https://doi.org/10.1029/2025JF008360) describes calibrated exposed-stake length changes at the thermistor pit. Authorized operator: exact exposed-stake distance increase[mICE]×917kg/m³=mmWATER equivalent. Ten-centimetre markings establish a scale, not camera precision.

Workflow:
1. Freeze thermal model, normalized predictions, source manifests, both train/evaluation intervals and all comparison rules. `camera_evaluation` in the specification must match the template schema, including each interval's role.
2. Publish specification and normalized-prediction hashes. Issue receipt with access_role=train and prediction_stage=normalized_training. `load_after_freeze` returns ONLY frozen training endpoint values, never the entire camera record. Train-end/evaluation-start is one shared observation, so its measurement error is shared.
3. Fit any predeclared training-only gain using those training endpoints. Publish final predictions and hashes before evaluation access.
4. Issue receipt with access_role=evaluation, prediction_stage=final_evaluation, and the final published prediction artifact. The loader returns ONLY evaluation endpoints.

Receipt fields bind frozen specification/prediction bytes and the sealed camera manifest; endpoints/operator must be identical inside the frozen specification. Publication chronology is externally evidenced by the publication reference, not cryptographically established by this local code.

No interpolation occurs at absent endpoints; missing/nonfinite endpoints refuse scoring. Negative endpoint differences and returned-point increments are flagged, not clipped/reversed/reset. Diagnostics describe released endpoints only, not full native chronology. Camera manual-correction/reset semantics remain a gate: observed resets invalidate affected intervals rather than allowing favorable replacements. Camera uncertainty is unestablished; no statistical CI is invented.

Seven tests use synthetic frames/artifacts only: unit conversion, negative differences, absent/missing values, pre-freeze refusal, train/evaluation role separation and artifact tampering. No native outcome has been examined or scored during preparation.
