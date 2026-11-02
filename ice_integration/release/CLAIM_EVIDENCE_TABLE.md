# Corrected claim–evidence table

Proposed additive correction, revision `corrected_evidence_v1`, based on pinned commit `05e7312e3bbf248acf9f646267a5b8ea06a8b32e`. Historical documents and results remain unchanged. See `CLAIM_REGISTRY.json`, `SOURCE_MANIFEST.json` and `check_claims.py`.

Evidence classes and statuses distinguish recomputed arithmetic, reported numerical verification, retrospective comparison, conditional demonstration and unsupported independent validation. Reported results have not been re-run by this correction. **No claim establishes field superiority or independent coupled-cliff validation.**

| ID | Claim | Evidence class | Status | Qualification |
|---|---|---|---|---|
| C01 | Point predictions reproduce V5 | numerical_equivalence | reported | Max relative difference4.2e-16 across5994values; no new physics or environmental validation. |
| C02 | Implemented heat/water accounting and restart checks | numerical_verification | reported | 15integration tests reported. Fluid-only closure does not certify phase-source water mass; atmosphere/cold ice/full geometry remain out of scope. |
| C03 | Plug-flow finite-heat limit | numerical_verification | reported | Integration check within1%; verifies its prescribed equations. |
| C04 | Offsets and support bookkeeping | numerical_verification | reported | Checks cover selected P3/P5 sources; unresolved clock conventions remain. |
| C05 | Consistent superiority over every baseline | retrospective_comparison | false | Candidate fails declared consistency; no recalibration to evaluation outcomes is authorized. |
| C06 | Nominal pooled error and scenario bias | retrospective_comparison | recomputed | Six nominal errors all positive, nominal MAE=bias0.584872726098m. Scenario marginal summaries must remain separate; engineering ranges are not confidence intervals. |
| C07 | Earlier local-geometry improvement | retrospective_development | reported | 23.6%within one family is a preserved development comparison, not consistent superiority. |
| C08 | Strict-native score | no_result | unavailable | No interval passes full native-support requirement; no score exists. |
| C09 | Station-specific calibration identifiable | training_diagnostic | unsupported | P5one low-SNRtraining interval does not identify a separate scale. P3leverage-equivalent count1.33948 is not independent statistical ESS. |
| C10 | Pressure head establishes full depth | input_diagnostic | unsupported | Clearance missing. Independent water-file station/time binding is hypothetical; mismatch exists only under that mapping, not proof which instrument is wrong. |
| C11 | Plan-view coupled demonstration | conditional_demonstration | reported | Imposed forcing, capped/settled and refined variants; melt-volume agreement does not certify contact convergence or field accuracy. |
| C12 | Wet/dry threshold sensitivity | numerical_sensitivity | reported | Bit-identical only for tested strip geometry and thresholds; no thin-film universality. |
| C13 | Mechanism predictor matches evolving geometry | model_to_model | qualified | Use latest raw-trace rescore and source-version precedence. One solid-strip loss history gives0.02day agreement; not three independent validations. Recovery contact and strength boundaries differ. |
| C14 | Natural porous-boulder behavior | closure_claim | unsupported | Ergun/drag/insulation closures unverified for real clast routing; omit from main environmental claims. |
| C15 | Rochemelon independent test | no_result | unavailable | Section/timing/thermal-path/geometry gates remain unresolved; no inputs invented. |
| C16 | Field cliff retreat/melt/lifetime accuracy | environmental_prediction | unsupported | No coupled cliff-model field validation supplied. |
| C17 | Independent validation against nature | independent_validation | unestablished | No independent-validation claim established in this package; retrospective comparisons are empirical evaluations but not independent confirmation. |

## Source/version rules

Use the raw-trace rescore for legacy ROI/recovery/history statements. Use solver versions/addendum for operator and boundary corrections, not as a completed full-v4 matched campaign. M3_FINAL is historical development evidence only where superseded. Different uncertainty calculations, endpoints and solver versions must not be pooled.

The source-reported removal classes are3maintained/1genuine reconnection/2failed, with2pre-loss histories and only1solid-strip loss comparison. The contact-persistence boundary is10.0–10.5days under the1-day rule; the strength threshold is about9.2days. A day10removal rewets for2.6days and loses contact again at day12.4. None is permanent recovery.
