# EXPERIMENT_STATUS — recovery milestone (M5)

**Verdict: NO-GO for a CFD/field campaign testing sustained recovery under the declared ideal closure.** The hypothesis "peak runoff can durably repair the bed-height gap within 21 days" is **falsified for the declared model** by a rigorous bound that the supplied screen and an independently written check both reproduce. No expensive simulation was launched; no parameter was inflated. Different physics (resolved heat paths, different local transfer, atmospheric melt, changing route geometry, natural floods) is **not excluded** by this result and would require independent justification and a new specification. The failed Swiss field-superiority verdict, all frozen results, the active branches and historical solver versions are untouched.

Basis `f1f85ee496dd39617fedf0f58202b943b3af3168`, reconciled with packet commit `10ee9706492160233e7b67cfcedbc29c5b2a3dd8` by fast-forward (worktree was clean). Branch `claude/recovery-m5` (isolated worktree `.worktrees/recovery`). Not merged.

## 1. Analytic viability (done before anything else)
| Item | Result |
|---|---|
| Packet MANIFEST.json | all 12 sha256 verified |
| `baseline/screen.py` rerun on a **scratch copy** (committed RESULTS.json untouched) | 12 cases, 0 viable; differs from committed JSON by ≤2.4e-13 relative (floating-point last digits only) |
| Independent re-derivation (`viability/independent_check.py`, no import of screen.py) | bound reproduced to 1e-16 relative in the bound; minimum repair time 3197 / 6393 / 12784 / 25566 d for Q0, 2Q0, 4Q0, 8Q0 |
| Full ODE of the ideal model, 12 cases + toe-area sweep 1e-3…100 m², dt 900 s and 450 s | repair in 21 d is 0 in every case (the gap grows); all ≤ bound; none reaches the critical gap (0.1889 m); dt difference ≤8e-12 m |

**What the bound is.** Gap g = z_toe − z_bypass; repair means the toe route lowers faster than the bypass. For equal-K critical controls and g ≥ 0, the bypass carries q_b ≥ Q/2. The finite-water floor melt rate e(q/A) = ρcqT(1−e^{−hA/(ρcq)})/(ρ_iLA) is increasing in q and bounded by the film limit hT/(ρ_iL) (28.2 mm/d at h = 200 W m⁻² K⁻¹, T = 0.5 K). Hence dg/dt ≥ −[hT/(ρ_iL) − e(Q/2, A_b)] throughout the interval, for any toe area and any gap. Multiplying by 21 d bounds *all possible* repair for a constant pulse; this is not an extrapolated initial derivative (the initial derivative is in fact negative: the gap grows). Required repair 61.1 mm; bound 0.40 mm (Q0) … 0.05 mm (8Q0).

**Assumptions that make it valid (all checked numerically in the independent script):** (i) two independent routes with equal critical-control coefficients, so the bypass share ≥ ½ for every gap ≥ 0; (ii) one fixed h and one inlet superheat common to both routes; (iii) each route receives its own incoming water heat ρcqT (no upstream depletion by the other route); (iv) floor melt only: no wall, atmospheric, dissipation, debris or conduction heat; (v) constant pulse discharge and temperature, no meltwater feedback on Q; (vi) quasi-static controls, no pool storage dynamics. **Why recovery fails here:** NTU = hA/(ρcQ/2) = 1.4e-3 ≪ 1, so both routes are film-limited and melt at almost the same rate; the bound is the film-limit times ≈ N/2 = 6.8e-4. Feasibility within 21 d would need NTU ≈ 0.22, i.e. 164 times the declared value. That is an analysis of which assumption controls the result, **not a proposal to raise h, T or pulse size**; larger NTU would have to come from independently justified resolved physics (long heat-limited reaches, small Q, different h per route), under a new specification.

**Analytic endpoint lane** (`analytic_lane/`, protocol frozen and committed in 6e6199e before scoring; 36 cases = 3 toe areas × pulse 1/2/4/8×Q0 × 6/24/72 h pulse, 48 h post-pulse window, threshold 1e-3, max sampling gap 900 s, packet scorer): 9 `initially_disconnected` (multiple 1), 27 `lost_without_reconnection` (rewetted during the pulse, peak toe fraction up to 0.28, lost afterwards), **0 sustained post-pulse recovery, 0 violations**; gap grows by 85–169 mm over the window; the fixed-geometry control gives the same classes (temporary-rewetting control, outside any denominator). The evolving-topography comparator with the same equations is algebraically identical to the candidate here, so agreement is **expected equivalence, not evidence of improvement**. Frozen-prediction note: my protocol predicted labels "temporary_or_incomplete_reconnection or initially_disconnected"; the scorer instead labels a pulse-rewetting from an initially dry route as `lost_without_reconnection` because its reconnection logic only counts reconnections after a first loss. The substantive prediction (never sustained, gap not decreasing) held; the label mismatch is recorded, and rewetting is therefore reported separately via `peak_toe_fraction_in_pulse`. Shared pre-loss histories: cases with the same toe area and multiple 1 produce identical sampled prefixes (shared histories, not replications).

## 2. Protected preflight (executed now)
`execution/preflight.py --mode physical_bare` on this isolated checkout, new output directory: exit 0, **66/66 pinned repair tests** (phase 19, contract 22, evidence 12, ROI 13), packet hashes and both installed phase-solver copies verified, source unchanged after execution, atomic solid/porous covered-source refusal probe passed (state arrays unchanged). `--mode conditional_quasi --covered`: passes and stays labelled conditional (no physical-time certificate). `--mode physical_bare --covered`: **refused** (exit 2, "drainage/storage closure required"). Existing tests executed now: `ice_integration` 15 + `research/readiness_2026_11_01` 20 = 35 passed; packet helper tests screen 4, scorer 10, preflight 5. Records: `verification/PREFLIGHT_*`. Exit 0 of the metadata checker means completeness only, never environmental validation.

## 3. Numerical verification (native physical-time, bare connected temperate ice only)
`verification/refinement.py` (20 s case, three independent controls; per-step maxima, not only final cancellation):
| Control | Levels | Max defect over all steps | Output difference to finest level |
|---|---|---|---|
| space dx | 0.25 / 0.125 / 0.0625 m | mass ≤6.2e-13 kg; sensible+latent ≤5.1e-9 J (≈3e-15 relative to ≈1.7e6 J) | melt volume −0.72 % / +0.04 % / ref; **toe discharge −45.6 % / −10.1 % / ref (not converged; the toe lane is one cell wide at dx 0.25)** |
| outer split step dt | 1 / 0.5 / 0.25 / 0.125 s | mass ≤7.2e-13 kg | melt ≤1.3e-3 %, toe discharge ≤1e-2 % |
| inner CFL | 0.4 / 0.2 / 0.1 | mass ≤6.1e-13 kg | melt ≤0.034 %, toe discharge ≤1e-3 % |
Fixed-ROI geometric melt minus heat-implied melt: |Δ| ≤ 8.2e-18 m³ in all 10 runs (absolute, tiny volumes, shown in fig2); positivity-added water exactly 0; per-route phase-source to sampled throughflow ratio 5.3e-4 (toe) and 3.3e-3 (bypass) at dx 0.25 (a sampled-section diagnostic, not proof of negligible feedback). **Limitation:** the delivered API ties hydraulic window, source epoch and bed update to one step, so the "morphology-update interval" cannot be varied independently of dt; it is covered by the outer-step control only. These are 20-second checks; they do **not** establish 21-day event accuracy or contact-time convergence (the earlier legacy-solver contact timing at dx 0.05 was ~10 % off; not re-established for v5).

## 4. Matrix: implemented / verified / assumed / unsupported
| Component | Implemented | Verified (executed now) | Assumed | Unsupported |
|---|---|---|---|---|
| Ideal two-route viability bound | yes (screen.py + independent check) | bound, ODE, adversarial toe-area sweep | assumptions (i)–(vi) above | other heat paths / resolved SWE / natural floods |
| Analytic endpoint lane | yes | 36 cases, no violations | declared ideal parameters (inspected, not blind) | any environmental inference |
| Native v5 bare physical-time solver | yes (delivered) | 19 phase tests, ledgers per step, 3-way refinement | temperate stationary ice, local melt | covered/porous meltwater, uniform-incision counterfactual (refused) |
| Covered quasi-steady runs | legacy only | preflight labels conditional | drainage transfer assumed | any physical-time certificate |
| Fixed-footprint exporter / scorer | yes (delivered) | 13 + 10 tests | float64, threshold 1e-3, gap ≤900 s | historical v4 bundles (not authenticated; gate stays failed) |
| CFD recovery campaign / evolving-topography comparator vs candidate | **not run (gate failed)** | – | – | – |
| Swiss field superiority | frozen | – | – | **failed; no new claim from reused intervals** |

## 5. Costs and limits
Preflight 21 s; refinement 1.5 s; analytic lane ≈ 8 s; independent check < 1 s; whole `reproduce.sh` 17 s on one thread (< 0.1 CPU-h against the 6 CPU-h cap; peak memory below 300 MB). No production job queued or running. A cost estimate for a full resolved v5 recovery campaign was **not produced**: the representative physical-time domain would need new specification after the gate, and spin-up/pseudo-time seconds must not be reported as environmental time.

## 6. Figures (opened and inspected; `figures/`)
fig1 viability (bound vs requirement; NTU dependence — the red line sits 3 orders above the declared-NTU bound); fig2 per-step conservation defects for the three controls; fig3 output refinement showing the non-converged toe discharge in space; fig4 endpoint classes and gap change. Ranges and differences are engineering sensitivities, not confidence intervals; GCI is not applied.

## 7. Decision and corrected direction
**No-go** on sustained-recovery testing in the declared closure; the temporary-rewetting outcome (toe rewets during a pulse, then is lost, gap unchanged or larger) is the *prediction*, and it was reproduced. Do not rename and repeat this hypothesis. Justified next steps, none launched: (1) establish independently whether any resolved thermal path (heat-limited long reach, route-specific h, partial depletion upstream, atmospheric or debris heat) can raise NTU by the required ~2 orders of magnitude under physical constraints — a new specification and a new viability gate; (2) resolve the toe-lane resolution problem (≥ 4 cells across the lane) before any partition result is trusted; (3) native v5 event-time and contact-time convergence on a justified bare case before any quasi-steady covered comparison.

## 8. Commands
`cd research/recovery_m5_2026_11_03 && ./reproduce.sh` (all of the above); figures `python3 figures/make_figures.py`. Baseline `screen.py` must be run on a copy (it overwrites RESULTS.json).
