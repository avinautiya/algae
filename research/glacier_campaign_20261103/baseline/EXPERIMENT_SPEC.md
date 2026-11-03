# Bounded experiment specification and viability decision

Sources read at commit f1f85ee496dd39617fedf0f58202b943b3af3168: readiness math audit, novelty audit and phase-source REPORT. These sources distinguish conditional covered quasi-steady runs from verified bare-ice physical-time steps. No unread newer branch is represented as audited. No Claude file or frozen result changed.

## Analytic gate before production

Twelve development cases are fixed in screen.py: toe area1/.5/.25m², bypass area1m²; discharge multiples1/2/4/8 of .07m³/s; initial gap.25m; equal critical-control K=.852447; commonh200Wm−2K−1,Tin.5K,21day horizon. The original geometry/thermal hypotheses are already inspected and these are NOT blind cases.

While g>=Hc(Q0)>0, equalK gives q_b>=Qpulse/2. The finite-water melt rate e(q/A) is increasing and bounded above by hTin/(rho_i L). Hence the instantaneous repair rate obeys, throughout the gap interval,

    E_toe − E_bypass <= hTin/(rho_i L) − e(Qpulse/(2 A_bypass)).

This upper bound is independent of toe area and gap within the declared interval. Multiplying it by21days bounds all possible repair during that constant pulse; it is not an extrapolation of the initial rate. The required repair is61.076mm. At2Q0 the bound permits no more than9.554micrometres/day and requires at least6392days;4Q0 and8Q0 require at least12784 and25566days. All twelve declared cases fail the21day durable-rescue feasibility gate. Equal-area cases also obey the stronger exact no-reset sign theorem.

**Do not launch a campaign expecting durable rescue in these cases.** Temporary hydraulic rewetting is possible. These conclusions do not rule out different full-SWE heat paths, different local transfer, atmospheric melt, changing route geometry, or natural floods. Adding such effects needs independent physical justification, not coefficient inflation to defeat this gate. The analysis itself is an elementary conditional limit, not historically certified novel mathematics.

## Supported verification lane A

Use native v5 physical-time mode only for bare connected temperate ice. Reproduce the19 phase tests and13 matched-ROI checks; short spatial/time refinement uses actual identical clock intervals, configuration/source hashes and separate mass/sensible/latent ledgers. Do not use spin-up seconds as evolution seconds. Covered solid/porous sources must refuse until a justified transport/storage closure exists. Long-duration contact convergence is not established by2second checks.

## Conditional mechanism lane B

Covered, morphologically accelerated legacy-compatible calculations can be labelled conditional quasi-steady demonstrations only. They cannot receive a native physical-time certificate. No phase-source result upgrades historical v4 outputs into authenticated v5 data. Distinguish maintained contact, loss, temporary reconnection, strength, post-window persistence and reloss. Keep missing v4 matched-support evidence as a failed gate.

## Fair comparators before any new development

Compare evolving topography plus geometry-dependent routing and thermal resistance, not only averaged crest and fixed geometry. Initial stage, route discharges, geometry, incoming water heat, initial per-route removed power and scored support must match. Record any fitted baseline parameter using only a nonzero-flow initial calibration state; a dry toe does not identify its heat coefficient. A baseline that changes initial crest while averaging shielding is confounded. Fit neither method to challenge outcomes.

If both methods use the same evolving controls and finite-water Stefan equations, their agreement is expected algebraic equivalence, not a novel improvement. Do not label a simplified comparator a reproduction of Strickland/Fastscape or a professional glacier model. Only use that label after actual equation/code reproduction and validation.

## Separate future test proposal — not launched

Only after a physically justified model passes the analytic viability gate, freeze9 withheld parameter configurations from a separately generated source manifest:3 initial gap ratios1.1/1.5/2 timesHc(Q0), crossed with3 pulse durations6h/24h/72h. Pulse amplitude and all heat/geometry coefficients must be independently justified and frozen BEFORE scoring, not selected to create recovery. A provenance check must establish these exact histories were never evaluated; otherwise label them inspected and create new ones. Include symmetric no-reset control and fixed-geometry temporary-rewetting control outside the primary win denominator.

Primary endpoint: sustained delivered heat over a declared48h post-pulse window, scored on identical ROI and clocks. Separately report first rewetting, total heat/melt and reloss. Freeze tolerances against independent numerical uncertainty BEFORE running, not against observed candidate errors. Predictions must include no recovery when bounded physics forbids repair. Use4 independent resolved cases only after cheap predictions and cost profiling, retaining failures/refusals.

New model support, endpoint reproducibility and source integrity are prerequisite gates; superior environmental predictions, historical originality, Nature-level significance and awards remain unestablished.
