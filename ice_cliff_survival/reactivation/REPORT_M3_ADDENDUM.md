# M3 addendum — identical-footprint melt, corrected recovery classes, removal-time bracket, differential-incision counterfactual

Scope reminder. Model class: 2-D depth-averaged free-surface hydrodynamics (hydrostatic) with coupled thermal incision. **Idealised obstacles only** (solid raised bed, or an Ergun-drag layer that is not resolved seepage). Nothing below is a prediction of natural boulder-deposit behaviour or of flood removal. Melt is local lane-reach loss, not cliff-wide or glacier loss. The reduced models and the CFD use the same series-resistance insulation closure, so agreement tests water delivery/connectivity/geometry, not the closure.

## 0. Corrections found during this work (read first)
1. **Solver defect found by the new static boundary check (v1/v2 affected).** `lower_bed` clipped *every* cell to `pool_floor`, so after the first bed update the free-overfall tail (floor −1.5 m − slope) was lifted to a flat −1.0 m: the outlet became a deep flat tail instead of a free overfall. Symptom: hydraulics never settled (storage up to 3 % of Q, f_toe drifting upward by tenths over 10⁴ s of hydraulic time in fixed-bed tests). Fix (solver v3): clip ice cells only. Effect on every run that does not remove debris early: flow fractions change by ≤2×10⁻⁴ (coarse v3 vs v1/v2, 8 of 10 E-runs); event times unchanged. **Effect on early removal (E05, E09): the v1 'decline' of f_toe from 0.67 to 0.19 after day 16 was an artefact; with v3 f_toe stays 0.664 (reduced model: 0.667).** v1/v2 outputs are preserved (`results/eval`, `results/v2`) and superseded by `results/v3`. Which v1 evaluation numbers must be treated as superseded: E05/E09 late-time f_toe, fig3 of REPORT_M3_FINAL (early-removal curve) and all hydraulics-storage diagnostics; classification and timing conclusions are unchanged. The v3 campaign repeated the v2 specification entirely (20 runs, coarse + medium) under a new freeze.
2. **Process error (recorded for provenance).** While parametrising the freeze script, a mis-patched script overwrote `frozen/production_spec_v2.*` and the incorrect hash was committed in f89eb09 (message with an empty hash); campaign jobs started from that rewritten spec were killed and discarded; the original v2 spec (sha256 458889f5…) was restored from 55ea17e (verified) in 28a3890; the v3 spec (sha256 4a471738…) was frozen and committed (142c0f0) before the v3 runs.
3. **Classification wording corrected.** v1's 'sustained reconnection' counted any f ≥ θ after removal, so it labelled runs that never lost contact as 'reconnected' (E03, E05, E09) and called E06's ≈0.1 flow fraction a reconnection. New classes (frozen in spec v2/v3), with R = mean f over [t_r+1 d, end] / bare-route share: maintained; lost-recovered-full (R ≥ 0.5); lost-recovered-partial (0.05 ≤ R < 0.5); lost-not-recovered (R < 0.05).

4. **Bed-operator clamp defect found later (solver v4; see `SOLVER_VERSIONS.md`).** v3 still clipped ice cells at `pool_floor`, so the entrance-ramp cells at the floor lost part of their melt: geometric bed loss was 18–21 % below the heat-implied volume (coarse 21 d). Reach and footprint melt (never at the floor) and the flow fractions are essentially unaffected: v4 differs from v3 by ≤5×10⁻⁴ (medium, no removal), ≤9×10⁻⁴ / 4.8×10⁻³ (coarse, none / removal 8 d) in f_toe and by +0.005–0.011 d in the loss day. v4 closes the heat–geometry volume to 1.4×10⁻⁹ (medium) and 1.8×10⁻⁹ (coarse). Ramp melt totals in §1 were heat-implied and are unchanged. The finest-grid boundary check was repeated on the v4 medium bed history (below).

## 1. Identical-footprint melt and heat (solver v3, coarse dx=0.05, 21 d)
Regions: footprint = deposit footprint cells of the toe reach (0.5 m², identical in CFD and reduced strip/corridor), open toe reach (0.5 m²; reduced 'bare corridor lane' for the corridor, not represented in the reduced strip model), bypass reach (1.0 m²). Entrance ramps (ice, 2.9 m², melting 0.81 m³ toe-side / 0.89 m³ bypass-side in 21 d, i.e. more than the reach) are in neither reduced model and are excluded. Melt numerical uncertainty (coarse vs medium, v1): ≤0.9 %.
* Bypass melt: reduced vs CFD +0.17…0.26 % in every run. Open toe lane of the corridor: +0.19 %.
* Footprint after prescribed removal (E05, E06): +0.15 % and −0.14 %, i.e. the earlier '−58 %/−98 %' melt errors were a **scope mismatch** (open reaches counted in the CFD figure), not model error.
* Footprint while the deposit is in place and no water reaches it (E02/E03/E04/E07): CFD delivers 0.4–0.6 W (melt 0.0023–0.0034 m³ in 21 d) through the tapered, partly wetted edge cells where the reduced model has exactly zero; absolute difference ≤0.003 m³ (≤0.5 % of bypass melt). Relative errors (−100 %, −50 %) there are not meaningful.
* Local heat delivered (mean over 21 d), CFD: bypass 99.6–99.8 W in all runs; footprint 0.0–1.6 W (deposit in place) vs 47.7 W (early removal); open toe reach 49.9 W (corridor, strip early removal) vs 26.7 W (strip lost: the open toe stretches lose their water supply when the strip locks out). Heat is reported next to flow contact: contact loss is accompanied by a 46 % drop of toe-reach heat delivery and a ≈ 96 % drop of footprint heat relative to early removal.
Porous cases have no reduced prediction.

## 2. Corrected recovery classes (v3 coarse, 21 d; v1 medium labels superseded)
| run | removal | class | R | decline | f_end |
|---|---|---|---|---|---|
| E01_bare_none | — | **maintained** | — | False | 0.732 |
| E02_corridor_none | — | **maintained** | — | False | 0.480 |
| E03_corridor_r12 | 12.0 | **maintained** | 0.657 | False | 0.480 |
| E04_strip_none | — | **lost** | — | False | 0.000 |
| E05_strip_r1 | 1.0 | **maintained** | 0.908 | False | 0.664 |
| E06_strip_r8 | 8.0 | **lost-recovered-partial** | 0.140 | False | 0.098 |
| E07_strip_r12 | 12.0 | **lost-not-recovered** | 0.000 | False | 0.000 |
| E08_porous_none | — | **lost** | — | False | 0.000 |
| E09_porous_r1 | 1.0 | **maintained** | 0.908 | False | 0.664 |
| E10_porous_r12 | 12.0 | **lost-not-recovered** | -0.000 | False | -0.000 |


E03: no loss occurred, but removing the corridor deposit at 12 d does not restore the bare share (R = 0.66): thermal memory persists even without contact loss. E05/E09 recover to R = 0.91 of the bare share. E06 (removal at 8 d) is a **partial** recovery (R = 0.14), not a reconnection to ordinary flow.

## 3. Removal-time bracket (NEW removal days, medium dx=0.025, solver v3, 15 d)
| removal day | f at +1 d (CFD / pred) | class CFD / pred | match | R | f_end | loss day | mean heat on former footprint W | footprint melt m³ |
|---|---|---|---|---|---|---|---|---|
| 8.5 | 0.0730 / 0.0513 | lost-recovered-partial / lost-recovered-partial | True | 0.095 | 0.0672 | 2.48 | 22.7 | 0.0941 |
| 9 | 0.0452 / 0.0233 | lost-recovered-partial / lost-not-recovered | False | 0.057 | 0.0393 | 2.48 | 20.9 | 0.0866 |
| 9.5 | 0.0218 / 0.0013 | lost-not-recovered / lost-not-recovered | True | 0.025 | 0.0150 | 2.48 | 19.1 | 0.0786 |
| 10 | 0.0044 / 0.0000 | lost-not-recovered / lost-not-recovered | True | 0.001 | 0.0000 | 2.48 | 9.5 | 0.0409 |
| 10.5 | -0.0000 / 0.0000 | lost-not-recovered / lost-not-recovered | True | -0.000 | -0.0000 | 2.48 | 1.3 | 0.0057 |


Reading: reduced boundary t* = 9.50 d (frozen); CFD boundary (f at removal+1 d crossing 1e-3, linear interpolation across the five new runs plus 8 d and 12 d) = 10.39 d, a difference of 0.89 d, inside the pre-registered ±1.0 d but late by almost the whole tolerance. Class predictions: 4/5 correct; the 9.0-d run recovers partially (R = 0.057, CFD) while the reduced model predicted none (its f(+1 d) 0.023 vs CFD 0.045). The classes are monotone in removal day. Local heat delivered to the former footprint falls from 22.7 W (8.5 d) to 1.3 W (10.5 d) while the footprint melt it accumulated shrinks accordingly.
**Finest-grid boundary check (fixed bed, steady hydraulics, bed history from the medium run; `boundary_check.py`, `fine_warm.py`).** A cold-start dx = 0.0125 solve (≈38 k cells) was abandoned after >24 min without finishing its first settle chunk (estimated >2 h per removal day); the finest grid was therefore run from a warm start (converged dx=0.025 free surface and discharges interpolated, then marched until Δf<1e-4 and |storage|<2e-4 Q; storage at exit ≤1e-7 m³/s). Steady toe fraction at removal (θ = 1e-3 crossing by linear interpolation):

| removal day | dx=0.05 | dx=0.025 | dx=0.0125 |
|---|---|---|---|
| 9.5 | 0.0343 | 0.0236 | 0.0287 |
| 10.0 | 0.0148 | 0.0064 | 0.0091 |
| 10.5 | 0.0021 | <1e-6 | <1e-6 |
| crossing of 1e-3 | 10.76 d | 10.42 d | ≈10.45 d |

**v4 repeat (bed operator without clamp, medium bed history `results/v4/V4M_strip_none`):** f at removal 9.5 / 10.0 / 10.5 d = 0.0241 / 0.0068 / <1e-6 (dx=0.025) and 0.0291 / 0.0093 / <1e-6 (dx=0.0125); crossings 10.43 d and 10.45 d (shift 0.02 d, criterion <0.5 d: pass). The v3 numbers below are superseded by these but differ by ≤0.04 d.

The boundary location is converged to ≈0.03 d between the two finest grids (pre-registered criterion: shift < 0.5 d: **pass**) and equals ≈10.4–10.5 d, i.e. the dynamic-run estimate (10.39 d) is confirmed. The earlier two-grid extrapolation toward 10.1 d is **not supported**. However the flow fractions themselves are not monotone in dx (0.0343 → 0.0236 → 0.0287 at 9.5 d; ±20 % between grids), so the classes within ≈0.5 d of the boundary (e.g. the 9.0-d run: partial vs none) are not resolved by the grid study; only the crossing day is. Limitation: the bed history is that of the medium run. The static f at 8 d (dx=0.025, 0.1055) reproduces the dynamic medium value (0.104), validating the static method after the v3 fix.

## 4. Counterfactual: differential incision disabled (same total melt spread uniformly over all ice cells)
| run | removal | contact lost | min f | f0 → f_end | toe-fp melt | bypass melt | crest toe → byp (end) |
|---|---|---|---|---|---|---|---|
| V3C_CF_corridor_none_uniform | — | False | 0.522 | 0.523 → 0.522 | 0.2665 | 0.5329 | -0.613 / -0.533 |
| V3C_CF_porous_none_uniform | — | False | 0.191 | 0.199 → 0.191 | 0.2671 | 0.5341 | -0.614 / -0.534 |
| V3M_CF_strip_none_uniform | — | False | 0.102 | 0.113 → 0.102 | 0.1908 | 0.3816 | -0.239 / -0.382 |
| V3M_CF_strip_r8_uniform | 8.0 | False | 0.107 | 0.113 → 0.735 | 0.2003 | 0.4006 | -0.481 / -0.401 |

No contact loss occurs in any of the four uniform-incision runs (strip: f 0.113 → 0.102; corridor and porous unchanged); removal at 8 d in the uniform strip restores full access (f → 0.735, bare share 0.735; reduced prediction 0.747). With differential incision the same strip loses contact at 2.5 d and recovers only 14 % of the bare share after removal at 8 d. **Differential incision is therefore necessary for both the loss and the memory in this model**; deposit shielding alone (the uniform runs keep the shield on the footprint) does not produce them.
The total melt volume is the same in both modes by construction; local heat delivery to the footprint is 0.19–0.20 m³ vs 0.0034 m³ (strip, no removal): uniform lowering also lowers the shielded cells, a deliberate property of the counterfactual.

## 5. What is and is not established
Established within the model: (i) the connectivity-aware reduced model reproduces the CFD contact-loss timing (0.01–0.02 d vs the medium grid; coarse-grid CFD is 0.23 d later), the corrected recovery classes, the identical-footprint melt (≤0.3 % wherever water reaches the footprint), and the counterfactual outcome; (ii) the removal-time boundary lies at 10.4–10.5 d in the CFD (converged between dx=0.025 and 0.0125 to 0.03 d; 10.76 d at dx=0.05) vs 9.5 d for the frozen reduced model — the reduced model is early by 0.9 d (9 %); (iii) recovery is partial, never full, for any removal after contact loss, and absent beyond ≈10.4 d.
Not established: natural boulder behaviour (the strip is an idealised obstacle; the porous layer is a drag closure), flood removal, 3-D/vertical flow, field insulation properties, behaviour beyond a 0.5-wide lane geometry with a 3 m entrance ramp, or any cliff- or glacier-scale melt.
Unresolved/limitations: flow fractions near the boundary are non-monotone in dx (±20 %), so classes within ~0.5 d of the boundary are uncertain; the boundary day itself is converged; classes near the boundary (9.0 d) are sensitive to a 0.02 difference in f; the Ergun closure constants are unvalidated; the entrance-ramp melt exceeds reach melt and is outside the comparison.

## 6. Cost (v3 campaign, measured)
20 runs (8 medium 15 d, 12 coarse 21 d), 4 workers: **wall 3821 s (64 min)**; medium runs 1526–1572 s each concurrently, coarse 240–270 s; peak RSS ≈ 155 MB. Fixed-bed steady solves (single core, warm-started): dx=0.025 ≈ 4–30 min per removal day, dx=0.0125 ≈ +7 min per day after the warm start; a cold start at dx=0.0125 was >2 h per day and was abandoned.

## 7. Next single physical test
Unchanged: fixed-bed flume head–discharge measurement of a surveyed boulder strip vs an equal-volume corridor (through-flow plus overtopping), because only that decides whether the drag-layer cases (E08–E10) or the solid strip describe natural debris — and therefore whether the removal-time boundary found here (≈10 d, set by the bypass having incised ≈0.2 m below the bump) applies to anything natural.

## 8. Reproduce
`freeze2 OUT v3` → commit → `campaign2 SPEC_DIR OUT 4 v3` → `analyze_v2 EVAL_V1 EVAL_V3 SPEC OUT v3`; `boundary_check V3M_strip_none_zb.npz OUT.json dx days`; figures `make_figures_m3` (figures_v3/ = fig1–6 from solver-v3 coarse E-runs; `figures/` fig1–6 are the superseded v1 set; fig7 = this addendum).
