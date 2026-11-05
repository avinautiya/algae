# Nature-style reporting and submission readiness

Checked5November2026against official Nature/Nature Portfolio guidance. Meeting reporting requirements does not establish editorial significance or acceptance. Requirements depend on the selected journal and article type; a life-science reporting form is not automatically applicable to this computational geoscience study.

Primary guidance:

- [Nature initial submission](https://www.nature.com/nature/for-authors/initial-submission)
- [Nature formatting guide](https://www.nature.com/nature/for-authors/formatting-guide)
- [Code publication guidance](https://www.nature.com/documents/GuidelinesCodePublication.pdf)
- [Nature Portfolio AI policy](https://www.nature.com/nature-portfolio/editorial-policies/ai)
- [Nature Communications policy guide](https://www.nature.com/ncomms/submit/policy) for a related journal's declarations; not a replacement for Nature's own article instructions.

| Item | Current status | Evidence or remaining action |
|---|---|---|
| Accurate prior-art positioning | Completed targeted audit; not exhaustive | PRIOR_ART.md; major debris/transfer concepts already published. |
| Scientific novelty/significance | Blocking for current superior-model claim | Robust calibration is established; new conditional transfer test fails. No Nature-level impact evidence. |
| Independent direct-melt test of revised model | Blocking | New external coefficients test completed, but target is sonic lowering, geometry conflicts and superiority fails. |
| Source/target operator descriptions | Documented with unresolved fields | OBSERVATIONS_AND_UNCERTAINTY.md; unresolved fields are limitations, not filled estimates. |
| Statistical sample/replicate definition | Documented |15development intervals, five campaigns;137zero-shift external days, one glacier/two seasons. Shared weather and endpoints disclosed. |
| Selection history and null results | Documented | Development, frozen Tasman comparison, failed attenuation/blending and external failure retained. No claim of selection-adjusted significance. |
| Equations and optimization | Documented/audited | PRIOR_ART.md and source; selected LP fits independently vertex-checked in earlier audit. |
| Data/code availability statements | Drafted, release incomplete | AVAILABILITY_AND_DECLARATIONS.md. Local files are not an archived public release. |
| Exact software/environment | Recorded for this module | ENVIRONMENT.json, requirements.txt; figure version pinned on regeneration. |
| Portable demo/audit | Executable | audit_release.py uses local TRAINING_COHORT.json and avoids ephemeral weather paths. Native scoring needs the verified ZIP. |
| End-to-end raw-source reconstruction | Partial | Original raw-to-feature pipelines still include an absolute/tmpKennicott path; portable feature-level audit is not full raw-source reproduction. Repair before release. |
| Immutable source manifest | Implemented | MANIFEST.json from build_package.py records included source/artifact hashes. |
| Versioned public code/DOI | Not completed | Author-approved repository release and persistent deposit needed. No publication or upload performed here. |
| Data redistribution licenses/citations | Partial | Source DOIs recorded; confirm license for each redistributed dataset. No blanket license inferred from open download access. |
| Figure source data and vector exports | Implemented for new external comparison | figures/ and figure_source.csv; earlier manuscript figures require final manuscript-specific inventory. |
| Main manuscript and supplementary text | Not completed in this task | Existing drafts are not a final paper incorporating this new failure. |
| Nature title/summary/length/figures | Pending manuscript choice | Apply target article limits when the actual manuscript is assembled. |
| Author list, affiliations and corresponding author | Author confirmation required | Do not invent roles or imply mentor approval. |
| Author contributions | Author confirmation required | Record actual intellectual, coding, validation and writing contributions. |
| Funding/acknowledgments | Author confirmation required | No unverified grant, institutional support or no-funding statement. |
| Competing interests | Author confirmation required | Each author must confirm. |
| AI assistance disclosure | Drafted; author verification required | AI contributed code, methodological suggestions and writing; disclose substantive use. Human authors must independently verify and own conclusions. |
| Ethics/permits | Applicability assessment required | Current analysis uses public geophysical measurements, no new human/animal study; confirm source permissions and journal-specific questions. |
| Cover letter and related manuscripts | Pending author input | Explain the actual contribution without a superiority claim; disclose overlap and concurrent submissions. |
| Journal-specific forms and all-author approval | Pending | Cannot be certified by this code package. |

**Readiness decision: not submission-ready as a newly superior general melt model.** A transparent transfer/calibration study can be drafted now, including negative external results. Reporting completeness must not be used to disguise missing novelty or validation.

The executable readiness validator reports blockers and exits2while any remains. It validates artifact existence/hashes separately from author attestations and scientific gates; a file's presence never proves scientific validity.
