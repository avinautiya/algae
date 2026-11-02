# Independent evidence-readiness review

Repository: `avinautiya/algae`; inspected commit `05e7312e3bbf248acf9f646267a5b8ea06a8b32e` on `claude/integration-m4`.

This is a read-only audit. No solver, frozen specification, historical output or Claude worktree was changed. Eleven source snapshots are preserved with exact Git blob hashes in `SOURCE_MANIFEST.json`. Run `python3 audit_evidence.py` to reproduce the arithmetic/provenance checks. The snapshots' hashes were all verified. This report does not certify literature priority or journal suitability.

## Decision

**No-go for a headline claim of consistent field superiority, independently validated cliff survival, or a Nature-level discovery. Conditional go for a bounded mechanism/approximation contribution, subject to prior-art distinction and the corrections below.** The current package has genuine verification and an informative retrospective field failure. Those are valuable but are not a shortcut to independent environmental evidence.

The largest readiness defect is not a solver crash. It is release prose combining legacy mechanism headlines, superseded recovery classes, different event definitions and different uncertainty sources into a stronger claim than the latest artifacts support.

## Concrete corrections, ordered by consequence

### 1. Stop citing the superseded mechanism report as the current evidence summary

`ice_integration/release/DECISION_MEMO.md` and claim13 of `CLAIM_EVIDENCE_TABLE.md` point to `REPORT_M3_FINAL.md`. Its early-removal dynamics, recovery labels and some model comparisons were superseded in `REPORT_M3_ADDENDUM.md` and the later `REPORT_M3_RESCORE.md`; `SOLVER_VERSIONS.md` identifies v4 as current. The stale source index must not bury the later raw-trace re-scoring.

- v1/v2 lifted the tail floor after bed updates; a late decline of the early-removal flow was an artifact.
- v3 fixed the tail but still discarded melt in ice cells reaching the clamp: geometric ice loss was18–21% below heat-implied loss in the stated coarse examples.
- v4 removes the clamp and closes that particular geometry ledger. The small reported changes in contact timing do not turn every v1/v3 result into a freshly verified v4 result.

Required correction: one current evidence index must name solver version, frozen spec, raw result set, event operator, physical footprint and superseded status for each quoted number. Preserve historical reports; do not silently rewrite their results.

### 2. The impressive loss-time result is one shared history, not three independent predictions

The legacy E04/E06/E07 trajectories share the same strip history before removal. Their identical2.478-day loss is not three independent demonstrations. The reduced2.49-day prediction is close to the medium-grid value for that history, but its0.01–0.02-day discrepancy is below the estimated discretization uncertainty.

The original ±0.37day is a conditional GCI estimate from development-grid loss times3.165/2.721/2.478days at dx0.1/0.05/0.025, with an observed order near0.86. It is not a statistical confidence interval, field uncertainty, or the uncertainty of the removal/recovery boundary. The v4 medium loss is2.4828days, not precisely the old2.4775day value.

Required correction: say “agreement within the numerical uncertainty for one shared pre-removal history,” report the number of distinct histories and keep post-removal tests separate.

### 3. The independent recovery-boundary comparison is less favorable than the release headline

The addendum reports a frozen reduced boundary9.50days versus dynamic-reference10.39days: error0.89day, nearly the entire predeclared ±1day tolerance. New removal-day classes match4/5, not an unqualified4/4 recovery success.

The current finest-grid values10.43/10.45days concern steady hydraulics evaluated on the **same medium-grid bed history**, with the fine solution warm-started. The0.02day shift is useful, but does not bound coupled bed-history error, interpolation error across0.5day removal samples, warm-start dependence or closure uncertainty. Flow fractions near the boundary change about20% and are nonmonotone with mesh. Classes within roughly0.5day remain unresolved.

Required correction: separate initial loss timing, recovery/nonrecovery, recovery fraction, dynamic removal boundary and fixed-history hydraulic crossing. Do not attach the old ±0.37day GCI to all of them.

### 4. Recovery classes must retain their corrected definitions

E03/E05/E09 maintained contact; they are not genuine reconnections. E06 loses contact and only partially recovers (R≈0.14). E07/E10 lose contact and fail the stated recovery criterion. The latest raw-trace rescore reports3 maintained/1 genuine reconnection/2 failed among six removal cases, with2/1/1 among the four nonporous cases having reduced forecasts. “Removal restores flow” and “restores the bare share” are different outcomes. Finite-horizon nonrecovery is also not proof of permanent impossibility.

The raw rescore adds a consequential distinction: bracket removals at9.5 and10days satisfy the one-day sustained-contact rule but fail the R≥0.05 strength criterion. Day10 reconnects for2.6days and loses contact again at12.4days. The contact-persistence boundary lies between10.0 and10.5days; the strength-defined boundary is approximately9.2days. A single undifferentiated “recovery boundary” conceals different physical questions. The4/5 bracket class result is for the declared R-based classification, not four demonstrations of permanently restored contact.

### 5. Confirmed numerical label error: nominal pooled bias is0.584872726m, not0.63m

The six signed nominal test errors in `figures/ERROR_BUDGET_NUMBERS.json` are all positive. Their arithmetic mean and MAE are both0.5848727260980938m. `ERROR_BUDGET.md` instead labels0.63m “nominal.” The value0.6309372253842787m is the linear-imputer marginal mean across27 nuisance scenarios, not the single nominal case.

Required correction: distinguish nominal bias, nominal MAE, scenario-averaged bias and the0.5237–0.8233m scenario range. The frozen failed-superiority verdict is unchanged by this prose correction.

### 6. Calibration influence is not a demonstrated causal explanation of test bias

P3:3 carries0.854863 of sum(basis²); leaving intervals out produces scales3.62265–9.58357. These establish leverage and instability under this fitting criterion. They do not prove which missing physics or input causes the positive test bias. The reciprocal concentration1/sum(weights²)=1.33948 is a leverage-equivalent count, **not a statistical independent-sample ESS** or a confidence calculation.

The published engineering experiment already retains all81 joint clock/temperature/width/imputer combinations. It supports the bounded statement that those **tested, jointly recalibrated scenarios** leave positive pooled bias. Marginal factor means should not be described as physical error-source contributions. Width is partly absorbed by per-scenario scale refitting. Other unquantified inputs and interactions remain outside that set.

Required correction: report sensitivity and calibration dependence; keep causal attribution unresolved. A small reconstructed-duration fraction alone is not a mathematical bound on the heat contribution of missing times, although the explicit tested imputer contrast is informative.

The “signed” budget also mixes reference directions. Its header defines+as increasing predicted incision, yet the viscous-dissipation row assigns−while saying adding dissipation adds melt. Cold-ice heat consumption uses the effect of adding an omitted process (negative predicted-melt change); the dissipation row instead appears to use the effect of its omission. The surrounding-ice-ablation row changes the observation rather than the prediction. Each entry needs an explicit baseline, counterfactual correction and delta(prediction−observation); otherwise the signs are not commensurable and cannot be added as a budget.

### 7. The independent-depth comparison does not establish which proxy is wrong

`FOLLOWUP_DIAGNOSIS_REPORT.md` correctly records the water-level file's P4/P5 ambiguity, hypothesized clocks, edge sampling, and floor interpolation over5–9day gaps during rapid incision. It then states the pressure proxy does not represent full depth consistently. The supported conclusion is weaker: the two reconstructions disagree **under the assumed station/time mapping**. Clearance, station assignment, clock, floor interpolation and water-surface sampling cannot be separated by this comparison.

Required correction: retain the discrepancy as a conditional diagnostic; do not calibrate clearance or attribute the mismatch to the pressure sensor from it.

### 8. Matched local melt is not total environmental melt accuracy

The addendum's identical-footprint comparison repairs the earlier −58%/−98% scope mismatch. It finds agreement within about0.3% for specified wet local regions under a common thermal closure. The latest raw-profile rescore uses a central0.5m² ROI and original medium exports: wet-case connectivity errors are within0.5%, but locked strips E04/E07 underestimate ROI melt by0.0020/0.0034m³ (−54%/−67%). Those absolute discrepancies exceed the stated5e-6m³ storage-precision scale by hundreds of times, even though they are small compared with the open-reach volume. They must be reported, rather than dismissed solely because the target melt is small. Storage precision is not total numerical or physical uncertainty.

The rescore explicitly labels raw-export provenance `metadata_identity_checked_not_run_authenticated`; the exports lack Q/T/commit identity. It refuses all10 v3/v4 rows because no separately versioned identity adapter exists. A complete v4 evaluation rerun is not done. Thus there is no authenticated all-case v4 matched-ROI accuracy table. The successful legacy comparisons remain important model-consistency results, not independent insulation validation or whole-cliff prediction.

The excluded entrance ramps melt more than the reaches and control hydraulics. Their exclusion is valid for a deliberately local scorer, but prevents extrapolation of local agreement into whole-system melt accuracy. The current-plan-view solver has no moving vertical cliff face/collapse model. Equal debris volume/footprint does not mean matched initial hydraulic resistance or identical initial flow distribution.

### 9. Uniform-incision intervention supports a bounded causal contrast, not a universal necessity theorem

The four counterfactual runs suppress loss when melt is spread uniformly. The math reviewer independently inspected `lower_bed`: it preserves the total assigned lowering **within each step of its own run** before redistribution. Once geometry feeds back into flow/heat, that operation does not guarantee equal integrated melt between independently evolved differential/uniform runs.

Required correction: use “redistributing the same step's melt uniformly eliminated the effect in the tested cases.” Do not assert an all-parameter necessity theorem or a physically realizable uniform-debris baseline. Archive full-system integrated heat and melt if claiming paired equality.

### 10. Eight contract failures are specific to this evidence package

The executed Swiss contract fails measured width, clearance, bulk-temperature relation, resolved clock, target-section discharge/conservation, discharge uncertainty, ice temperature and measurement independence. These are legitimate gates for this model's intended field claim. They do not prove all public datasets or future observations are incapable of validation.

The report's “new survey is the only unconditional go” should mean that the listed existing sources do not presently pass these gates. A proposed new survey is not automatically adequate until its actual measurement support and uncertainties pass the same checks.

## Current environmental evidence ceiling

- The integration candidate has pooled nominal MAE0.585m, worse than Gnielinski0.443m and constant-rate0.532m on six inspected intervals. The failed declared criterion must remain.
- The field split is retrospective and severely calibration-leveraged. A new parameter fitted to P5 evaluation outcomes would not repair validation.
- Two-grid coupled demonstrations show melt-volume agreement around1.5% in the settled subset but appreciable routing/time errors. Conservation does not certify hydraulic steady state or contact accuracy.
- Wet/dry sensitivity is a useful negative result for this strip geometry only.
- Solid and drag-layer obstructions are idealizations. The porous closure, under-rock paths, thermal closure and marginal hydrostatic crest adequacy are not established by their agreement with a reduced model sharing those assumptions.
- No present result demonstrates whole-cliff lifetime, Himalayan glacier transfer or professional-model superiority. This audit does not rule out future evidence or a scientifically useful narrower paper.

## Feasible decisive comparison before a large campaign

The best next synthetic discrimination is **not** another count of the same loss trajectory. Specify a small blinded set comparing the connectivity predictor against a stronger arrangement-blind approximation matched to the same initial head–discharge relation, initial flow fraction, thermal power and total debris amount/footprint. Then vary history/arrangement and score continuous flow, local matched-footprint melt, initial loss and removal-boundary location separately. Use v4 throughout; keep common closures explicit; include non-events, failed recovery and near-boundary cases.

Before calling that an environmental mechanism, bracket transmission/underflow and closure ranges from independent evidence, and bound the critical crest hydraulics against an appropriate independent reference. A favorable result would show added predictive information beyond initial resistance under specified physics. It would still need the separate prior-art gate and a plausible environmental parameter regime; it would not itself establish field accuracy or Nature-level significance.

## Exact source index

See `SOURCE_MANIFEST.json` for all URLs and blob hashes. Key pinned sources are:

- `ice_integration/release/DECISION_MEMO.md`, blob6b54eb6518f1b6c854f48b1717443d111eb209e6.
- `ice_integration/release/ERROR_BUDGET.md`, blob7c8c4d44eb2eb63c27b29ee265a229420aee8d1d.
- `ice_integration/release/figures/ERROR_BUDGET_NUMBERS.json`, blob122efbae308e6a9b920cbaa4fd930d38ad08b244.
- `ice_cliff_survival/reactivation/REPORT_M3_FINAL.md`, blobafb6e395a73eb86cb24cf781ec92ffb7fa66e12b.
- `ice_cliff_survival/reactivation/REPORT_M3_ADDENDUM.md`, blobe9bc6ef0ea65d5ff6e99f8218cb981d739f8391e.
- `ice_cliff_survival/reactivation/REPORT_M3_RESCORE.md`, blob32f69c3603310375ecdaee9be405b19e0be5a1d5.
- `ice_cliff_survival/reactivation/SOLVER_VERSIONS.md`, blob4cfc0510724f521853a1e03b186e95264293667a.
- `FOLLOWUP_DIAGNOSIS_REPORT.md`, blob7d2e25798c70320fd1b3dc9c3013b620852901a6.
