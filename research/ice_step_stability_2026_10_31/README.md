# Ice-step research foundation

Prepared 31 October 2026. This package contains a targeted primary-literature audit, a scoped research decision, first-principles derivations and lightweight numerical verification. No production campaign was run and no remote repository or active Claude process was modified.

Start with **RESEARCH_DECISION.md**. The recommended question is whether dynamical stability, rather than an assumed maximum wavelength, can predict persistent water-carved ice steps and their melt concentration. Novelty remains provisional and the coupled stability result has not been computed.

Files:

* `RESEARCH_DECISION.md`: comparison of eight directions, closest gap, environmental relevance and falsification conditions.
* `MATHEMATICS.md`: conservation model, thermal response/passivity proof, flat diagnostic dispersion system, traveling profiles, moving-jump conditions and stability formulation.
* `SOURCES.md`: primary links, reading/access depth and remaining priority checks.
* `CLAUDE_HANDOFF.md`: substantial but bounded next implementation milestone.
* `verify_foundations.py`: standard-library verification calculations.
* `verification_results.json`, `thermal_refinement.csv`, `illustrative_timescales.csv`: reproducible outputs; illustrative quantities are not field results.

Run `python3 verify_foundations.py` from any directory. The checks establish agreement for selected mathematical subproblems, not environmental validity, full-model stability or originality.

`sources/` retains two open university-hosted manuscripts locally for audit. The compact ZIP excludes those third-party PDFs and text extractions; links remain in SOURCES.md. SHA256SUMS.json records the handoff file contents.
