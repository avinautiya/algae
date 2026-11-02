# Field thermal validation protocols

No model fitting or predictive field scores were computed by this protocol task.

- `TEMPERATURE_VALIDATION_SPEC.md`, `TEMPERATURE_SPLITS.json`, `TEMPERATURE_PROTOCOL_FREEZE.json`: eleven real-sensor subcolumns, chronological train/test splits and four prespecified predictors. Negative temperatures are retained.2020b has only two measured probes and is excluded from this primary interior-target test.
- `profile_input.schema.json`: input/provenance contract. Semantic checks in its description are mandatory in addition to schema syntax.
- `ABLATION_TRANSFER_SPEC.md`, `ABLATION_PROTOCOL_FREEZE.json`: separate, conditional source-to-other-stake conductivity transfer. Preserves missing leading time, survey-date and depth/representativeness uncertainty. Covered-period one-sided consistency is distinct from full-period imputed point accuracy.
- `IDENTIFIABILITY_AND_FLUX.md`: why alpha-only temperature validation does not establish absolute melt and why seasonal heat can be insensitive to thermal dynamics.
- `thermal_identifiability.py`, `test_thermal_identifiability.py`: independent harmonic/reference identities; five tests pass.

Run the mathematical reference checks:

```sh
python -m unittest discover -s . -p 'test_*.py'
```

The data/measurement audit is in the parent `field_thermal` folder. Source hashes and exact UTC splits are frozen before model scoring. These are retrospective protocols, not externally registered blind studies. Do not choose a model, origin, clock/depth/imputation scenario or exclusion after inspecting test errors.
