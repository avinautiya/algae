# Typed metadata gate, version2

The old checker at the pinned integration commit approves invalid values because it tests truthiness rather than types or physical domains. This isolated replacement rejects string booleans, boolean numbers, NaN/infinity, invalid fractions, empty evidence, unknown fields, unsettled hydraulics, missing numerical records and residuals above the supplied predeclared limits. It rejects duplicate JSON keys. Twelve tests pass.

Each field has `value` and an `evidence` record with HTTPS sourceURI, SHA256, locator, derivation and status. `METADATA_TEMPLATE_V2.json` is deliberately blank and fails. No field is inferred. The schema is intentionally incompatible with the old untyped template; migration must preserve unknowns instead of filling them.

**Crucial:** evidence strings/hashes are not authenticated or read by this checker. A dishonest self-assertion can still be syntactically complete. Therefore output is `metadata_complete` and `scientific_review_required`; `may_issue_environmental_prediction` is ALWAYS false. Exit0 only means metadata completeness, not field validation, model applicability, independent calibration or predictive superiority. Source review and the actual model-validity assessment remain human/scientific tasks.

The gate checks whether numerical results and predeclared limits are documented, not whether the limits are scientifically adequate or the stated convergence estimate is correct. Do not present it as a foolproof model-certification system.
