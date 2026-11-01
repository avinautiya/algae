# Finalisation and forecast comparison support

This package edits no model and runs no CFD. It provides a strict **self-reported
metadata checklist**, not independent scientific approval, and descriptive
paired forecast scoring that preserves event censoring.

From the parent directory:

    python3 -m unittest discover -s audit -t . -p 'test_*.py' -v
    python3 -m audit readiness audit/example_readiness.json
    python3 -m audit evaluate audit/example_evaluation.json

Every example is synthetic. A complete checklist only means the required
declared records exist and have passed statuses. Hashes identify records; the
tool does not retrieve or inspect their scientific contents. Environmental
prediction claims additionally require declared independent external validation.
No output approves historical novelty, glacier-wide impact or publication.

## Input and metrics

Readiness requires a pinned commit, frozen configuration, complete unique case
sets, no unresolved failures, and hashed passing records for hydraulic partition,
water/thermal conservation, spatial/time convergence, event timing, boundary
applicability and closure sensitivity. The current reported partition failure
must remain a failed gate until resolved or a new justified applicability scope
is frozen; never change its status merely to complete a checklist.

Evaluation uses one explicit operational definition for every model:
toe-flow fraction threshold, persistence duration, and assessment after return
to ordinary flow. `recovered` is a Boolean or null (unresolved label or model
abstention). Scored accuracy explicitly excludes abstentions; their counts are
reported. Pairwise comparisons use only common nonabstained labelled cases.
No independence, confidence intervals or p-values are assumed.

Event observations are either an interval containing the first sustained
recovery start or right-censoring at the last adequately observed time. A
forecast inside an interval has zero distance to the bracket; this is **not**
zero true timing error. Predictions before a right-censoring bound violate
observations; later predictions are unobserved, not correct exact times. No
MAE is manufactured for censored observations. Null event prediction denotes
no event. First contact and sustained recovery are different targets; do not
score a first-contact forecast as if it predicted sustained recovery.

## Tonight's order of work

1. Finish the existing hydraulic partition/reference diagnosis. Preserve failed
   cases and avoid a geometry change whose sole purpose is passing a tolerance.
2. Verify event scheduling, actual removal timestamps and the physical-time
   water/thermal/latent ledger; quasi-steady spin-up is a different clock.
3. Resolve spatial and morphology-step event convergence on development cases.
   Keep operational contact threshold and persistence sensitivity visible.
4. Freeze a bounded test set and predictions before test execution. Include a
   recovering case, a nonrecovering case and a never-disconnected control.
5. Compare complete clearing on evolved ice, the same clearing without prior
   differential incision, and the exact evolved bed restarted without debris.
   Preserve porous and lateral routes. This isolates substrate memory.
6. Score meaningful baselines: peak-only, storage-only, leakage-aware water
   bound and a geometry-aware time-resolved predictor. Do not equate beating a
   static baseline with superiority to all glacier models.
7. Produce figures from saved records: ordinary-flow recovery boundary, contact
   histories with actual removal/return times, bed changes, local delivered
   heat, forecast error, and numerical/closure sensitivities. Label synthetic,
   development and held-out results separately.
8. Write the claim that the result supports. If a gate remains failed tonight,
   finish a reproducible diagnostic release rather than label invalid data final.

The tool cannot make a nightly completion deadline scientifically safe. It can
make unfinished evidence explicit and prevent silent promotion of development
results into final environmental claims.
