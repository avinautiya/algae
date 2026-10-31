# FREEZE — event-aware averaging study (written after development, before any test run)

Development evidence: D01–D07 (network r17, 35 runs, `eval_dev.json`) and sequence-study X09 (stored outputs, `x09_diagnostics.json`). Dev findings that shaped this freeze:
mean forcing reproduces wet hours in several intermittent cases but under-predicts delivered volume of outlet 2 by a factor 2–10; the water-storage-only model SE overestimates wet time (sill lowering moves the threshold); uniform fine blocks (PW3) and the reduced thermal model driven by the series (RTM-TR) are as good as the event-aware segmentation EA on dev. **No method parameter was changed after seeing dev results other than adding PW3 as a baseline.**

## Frozen definitions
* Series: 36 h hourly (`series.py`), initial stage 0.10 m above sill 1, outlet 2 initially dry/wet according to gap.
* Wet outlet 2: Q2 > 0.5 % of the mean inflow. Class: dry if wet < 5 % of the time, wet if > 95 %, else intermittent.
* Truth = network time-resolved run (TR, r17, `atol_energy` 1e-3 J). Candidates: AVG, PW6 (6 h blocks), PW3 (3 h blocks, test-only baseline), PWn (as many equal blocks as EA), EA (blocks bounded by the wet/dry runs of the water-storage model SE, runs shorter than tau_c merged), SE and RTM-TR (reduced model, RTM_seq with the series).
* Criterion A: `R = CUSUM_max(wrong-side inflow volume about I*) / (A |u_ss(mean inflow)|)`; unsafe iff R >= 1. Rigorous parts: comparison-theorem bounds (`guaranteed`), exact I* = Q1(s2) for fixed geometry. R and its cutoff are an asymptotic screen, not a theorem.
* Test configurations: `configs_test.json` (E01–E10), cases in `cases_test.json` (hash in `FREEZE.sha256`), drawn with seed 71 from draws 100–199, quotas by category from the RTM screen only (no network run).

## Predeclared metrics and success rules
Metrics (vs network TR): class agreement; |wet-hour error| (h); error of outlet-2 delivered share of total delivered volume (pp, `share_deliv`); relative error of outlet-2 volume (reported where TR volume > 1 % of total); final-share error (pp, reported separately, not a success metric).
"AVG misrepresents" (truth label) := network AVG class != TR class, or |wet-hour error| > 2 h, or |share_deliv error| > 5 pp.
1. Criterion A succeeds if recall >= 0.8 and precision >= 0.8 on E01–E10 with >= 3 positives and >= 3 negatives; otherwise "not established".
2. EA is a *useful improvement* only if, on truth-positive configs, median |share_deliv error| of EA <= 0.5 x AVG's **and** <= 0.8 x min(PW6, PW3), and EA is worse than AVG on at most one config. If it fails the second condition the verdict is "no improvement over simple piecewise forcing", whatever EA does against AVG.
3. RTM-TR is judged as a cheap predictor by the same share_deliv metric against network TR (success: median error <= 2 pp on truth-positive configs).
4. Robustness: ref200 (4 configs; TR/AVG/EA) and nu_factor 0.5/2 (3 intermittent configs; TR/AVG/EA) must preserve the class of every method and the AVG-vs-EA ordering for a conclusion to be called robust; legacy-tolerance (atol_energy=None) TR on 2 configs reports the numerical sensitivity.
5. Cost counts preprocessing: tabulation (4 tables), SE run, method forcing construction, and the simulator run.
Failures and counterexamples are kept in `SCORES_TEST.md`.
