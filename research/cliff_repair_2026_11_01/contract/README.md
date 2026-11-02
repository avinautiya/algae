# Runnable metadata-gate repair

This isolated package proposes the replacement of `ice_integration/release/validity_contract.py` at repository commit `05e7312e3bbf248acf9f646267a5b8ea06a8b32e`. No current repository source, Claude worktree, calculation or frozen scientific result was edited. The four exact original source files are preserved in `legacy/`; their Git blob hashes are authenticated in `VERIFICATION.json`.

## What changed

- Adds typed values, explicit physical-domain checks, source-evidence records and the documented numerical-record group.
- Rejects the seven demonstrated legacy exploits; false strings, boolean uncertainties, negative/nonfinite numbers, empty containers, impossible fractions, missing numerics and unsettled hydraulics cannot pass.
- Rejects malformed evidence URLs, enormous numeric values, unknown fields, duplicate JSON keys and excessively nested JSON without a traceback/approval.
- Preserves `check(meta)` result keys `may_issue_environmental_prediction`, `unmet`, and `label_if_not` for callers, while adding `metadata_complete`, `issues`, and `scientific_review_required`.
- Never grants environmental scientific authority. Even fully typed synthetic metadata has `may_issue_environmental_prediction=false`. A schema checker cannot authenticate its evidence, independently review a source, establish identifiability or prove predictive accuracy.

The output `unmet` contains unique field names; `issues` retains distinct value/evidence reasons. Input bytes receive SHA256 fingerprints. The CLI accepts metadata up to8MiB, rejects duplicate keys and nonfinite JSON constants, and writes new outputs exclusively without replacing existing files.

## Execute

```sh
python3 verify.py
python3 validity_contract.py --template blank_v2.json
python3 validity_contract.py blank_v2.json
# Expected exit2: blank evidence and values remain unknown.
```

Explicit migration example (use fresh destination paths):

```sh
python3 validity_contract.py --migrate legacy/swiss_package_metadata_status.json \
  --output swiss_v2.json --migration-report swiss_migration.json
python3 validity_contract.py swiss_v2.json
# Migration succeeds as an operation; the check still exits2.
```

The supplied `swiss_package_metadata_v2_unverified.json` and `swiss_migration_report.json` demonstrate this migration. Original `true` placeholders for numeric floor uncertainty and survey dates remain wrong-typed `true`; no numerical uncertainty/date is invented. False independence remains false. Every evidence record remains blank, and numerics remain unknown. Unknown source fields/groups are retained in the migration report's full `legacy_document` and `unknown_fields`. Migration never upgrades an assertion to an observation.

## Integration instructions

In a NEW versioned integration change, place `validity_contract.py` and `validity_contract_v2.py` together in `ice_integration/release/`. Preserve the existing release package `__init__.py`; it need not change. Do not overwrite historical checker/results: archive them or retain their existing commit as this package does. Add the V2template and migrate old status documents explicitly. The tested module command is:

```sh
python3 -m ice_integration.release.validity_contract METADATA.json
```

For an actual check, exit0 means typed metadata complete ONLY; exit2 means incomplete/refused input. Template generation and migration exit0 mean those operations succeeded, not readiness. Existing callers relying on truthy legacy metadata or environmental approval must update; there is intentionally no permissive fallback. Legacy results remain reproducible with `legacy/validity_contract.py`, but that checker is demonstrably unsafe and must not authorize new claims.

Migration writes its report and metadata as two separate exclusive files. A filesystem failure may leave one partial output set, but never a success/approval; retain the error and choose fresh output paths. These writes are not a two-file transactional database.

## Verification and limits

`verify.py` runs22tests: the prior12domain checks plus10adapter/legacy-exploit/migration/CLI/package checks. It verifies four preserved source blob hashes. All pass. The package-entrypoint test builds an isolated repository-shaped package, so standalone success is not substituted for module-import compatibility.

Evidence sources are not downloaded or authenticated by `check`; their fields are assertions awaiting review. Predeclared limits may be scientifically poor even if syntactically valid, and declared hydraulic/grid results may be wrong. No self-asserted metadata becomes verified science here. Missing Swiss observations, failed model accuracy and historical novelty remain unresolved by this software fix.
