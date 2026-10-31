# Glacier-channel predictive-model foundation — V2

Prepared 31 October 2026. This revision makes an end-to-end forward predictor the main deliverable and treats ice-step stability as its proposed scientific mechanism. It preserves the earlier foundation while correcting the emphasis on cold-ice memory and specifying a CFD/validation route.

**No full CFD solver or new environmental prediction was executed in this package.** The executable checks are analytical design/verification calculations. No remote repository, existing campaign or active Claude job was changed.

Read:

1. `PROJECT_V2.md` — scope, outputs, novelty gates, comparison study and environmental limits.
2. `CFD_AND_VALIDATION.md` — solver capability requirements, verification and evaluation procedure.
3. `MATHEMATICAL_UPGRADES.md` — closure counterexample, energy exchange, phase-change mass, conditional error bound and hydraulic causality.
4. `REFERENCES_AND_DATA.md` — additional primary literature and actual data-access status.
5. `CLAUDE_PROMPT.md` — revised next implementation prompt.

Run `python3 prediction_basis.py` to reproduce `design_check_results.json` and `temperature_counterexample.csv`. The script needs only Python 3.9+ standard library. The temperature profiles are constructed counterexamples, not measured glacier conditions or an observed improvement.

The earlier packet is at `../ice_step_stability_2026_10_31/`; its mathematical derivations and verification remain relevant within their stated scope. The V2 archive includes a copy of that compact foundation archive for convenience. Original third-party PDFs are excluded. Source files and package checksums are recorded in SHA256SUMS.json.
