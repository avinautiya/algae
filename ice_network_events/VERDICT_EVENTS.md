# Verdict — threshold-aware averaging for outlet activation/drying (bounded milestone, then stop)

**Short answer: a useful application / approximation, not a new mathematical result.** The criterion works as a conservative screen but fails the frozen precision rule; the event-aware segmentation (EA) did **not** beat a simple fine piecewise-forcing baseline; the most reliable cheap predictor was the already-existing reduced model driven by the time-resolved series.

## What was done
Read: Kingslake 2015 (all), Jarosch & Gudmundsson 2012 (all), Ogier 2021 (pp. 1–15), Schoof et al. 2023 JFM (pp. 1–10 only; rest and supplement unread). Threshold-averaging literature known at abstract/textbook level only (`docs/event_prior_art.md`). Dev: 7 configs (D01–D07, 35 network runs) + X09 stored outputs. Frozen before testing (`FREEZE.md`, commit 1004f58): methods, metrics, success rules, 10 test configs (E01–E10). Test: 80 network r17 runs, 12 resolved ref200 runs (relaunched once after an interface failure, `FAILURES.md`), RTM runs; ν×0.5/×2 on 3 configs, legacy-tolerance TR on 2.

## Results (frozen rules; `SCORES_TEST.md`)
* **Truth-positive configs (AVG misrepresents): E03, E04, E07 (3 of 10).** E03: AVG dry, TR wet 4.2 h; E07: AVG wet 35 h vs TR 22.9 h (wrong class); E04: AVG 5.3 h wet-time error and 17 vs 340 m³ in outlet 2.
* **Criterion A (R ≥ 1): recall 1.0 (3/3), precision 0.6 (3/5) → not established under the predeclared ≥ 0.8 rule.** False positives E05 (mean at the threshold, R = 516: R diverges as the mean margin → 0) and E06: the screen uses the initial-geometry I*, but the network stayed dry because sill lowering moved the threshold — the screen is conservative and blind to geometry drift. 5 true negatives incl. the three wet/dry-far cases and E08.
* **Mean forcing:** wet-hour error 4–12 h on positives; outlet-2 volume error median 68 % (positives), although the share of total delivered volume changes only 0.7–6 pp because outlet 2 is a small fraction of the total.
* **Event-aware segmentation (EA):** median share error 0.74 pp vs AVG 4.75, PW6 0.66, PW3 0.32. EA beats AVG (rule part 1 ✓.) but **not** the simple baselines (rule part 2 ✗) → **no improvement over simple piecewise forcing; PW3 (3 h blocks) is best.** Median relative volume error: AVG 0.68, PWn 0.33, EA 0.14, PW6 0.08, PW3 0.04. The same held on dev (PW6 ≈ EA).
* **RTM-TR (reduced model driven by the series):** median share error 0.14 pp (rule 3 ✓), wet-hour error ≤ 0.4 h on E03/E04/E07 against the network. **Counterexample kept:** E08 (R = 0, TR intermittent): RTM-TR wet 26.6 h vs network 34.1 h, share error 6.1 pp — the reduction misses a case where the network's sill erosion is large (256/278 mm).
* **Water-only model SE (fixed geometry) is not adequate:** wet time 7–23 h where the network time-resolved run has 0–8 h (D03, D05, E03, E05, E06) — sill lowering moves the threshold; this is the part of the criterion that cannot be closed without the thermal model.
* **Numerical/physical uncertainty (rule 4):** classes preserved under ref200 in all 4 configs; outlet-2 share changes by up to 1.6 pp (E04 TR 6.2→4.6 %). Heat-transfer ×0.5/×2: class flips for TR in E03 (×0.5 intermittent share 2.2 %, ×2 dry) and E05 (×0.5 gains a wet period) — **near-threshold wet/dry classification is not robust to the closure**; E04 (clearly intermittent) preserves class and the AVG-vs-EA ordering. Legacy tolerance gives identical shares to 7 digits in E03/E04 (the stall issue did not arise).
* **Cost (single core):** thresholds + tables 1.4 s (+0.8 s fixed tables), criterion 8 ms, SE 0.3 s, RTM-TR 0.4 s ⇒ ≈ 2 s total including preprocessing; network r17 run 12–16 s here, ref200 270–400 s. Averaging does not make the network run cheaper (same cost for TR and AVG).

## What information must be retained beyond mean water and heat
Inflow volume and enthalpy per run above/below I* at resolution ≥ τ_c (≈ 1–2 h here), their order, and the sill-drift rate. Counterexample and proof in `docs/event_theory.md` §1.6; X09 and E03/E07 illustrate it.

## Novelty assessment
* Rigorous parts (threshold I* = Q₁(s₂), comparison theorem, deficit bound, exact crossing time) are standard ODE/reservoir comparison arguments applied to this lake; no result beyond them is claimed.
* R is an asymptotic/fitted screen; its cutoff 1 passed recall but failed precision on the frozen test.
* EA and PWn have no demonstrated advantage over 3–6 h block means; a user can obtain the benefit simply by using sub-τ_c-resolved or 3 h piecewise forcing, or the RTM-TR reduced model.
* **Outcome: useful application/approximation (diagnosis of why means fail and a conservative screen), not a new mathematical result.** The paper-level overlap table finds no source that treats averaged vs time-resolved forcing for outlet activation, but the tools are known.

## Limits and failures preserved
3 positives only (low power; recall/precision are not statistically tight); 36 h synthetic series; two outlets; quasi-steady ratings; closure-dependent near threshold; RTM misses E08; screen blind to geometry drift; ref200 shows a 1.6 pp numerical uncertainty at the 5 % level. Not a real-glacier forecast.

## Is further development justified?
Only one narrow item: add a geometry-drift term to I*(t) in R (the false positives E05/E06 and SE failures all trace to it). Broader method development is not justified: simple 3 h piecewise forcing already achieves < 1 pp share error.

Files: `design.py`, `methods.py`, `series.py`, `sim.py`, `evaluate.py`, `score.py`, `cost.py`, `figures_events.py`, `figures/fig1–fig4`, `eval_*.json`, `scores_*.json`, `SCORES_TEST.md`, `FAILURES.md`, `docs/`.
