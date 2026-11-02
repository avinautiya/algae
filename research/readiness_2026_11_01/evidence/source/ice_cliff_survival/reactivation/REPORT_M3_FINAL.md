# Final report — thermal memory and bypass incision: connectivity vs reduced forecasts (M3)

Model class: **2-D depth-averaged free-surface hydrodynamics (shallow-water, hydrostatic) with coupled thermal incision** (advected sensible heat; the same exchange flux lowers the ice floor; debris insulation by a series thermal resistance). It does not resolve vertical flow, grain-scale underflow or 3-D turbulence; the porous deposit is an Ergun-drag closure, **not** resolved seepage. Prescribed debris removal is not a flood-entrainment prediction. Results are local toe-lane melt, not cliff-wide or glacier-wide loss. All inputs are assumed test values.
Imported reduced package: `imported/ice_stream_reactivation_3e51620` (pinned 3e51620…, unmodified). Freeze: `reactivation/frozen/production_spec.json` (sha256 5b4938c6…), committed (7169e47) before any E-run was launched; solver commit b19e4f2.

## 1. Partition discrepancy — diagnosis (original failures preserved: `results/verify_swe_v4.json`, `diagnose_v1*.json`)
The strip error (+0.11, unchanged by refinement) was **not a strip defect**. Fixed-bed isothermal tests (`diagnose.py`) on identical symmetric lanes showed the discharge coefficient against the ideal broad-crested weir depended on the entrance step height (0.86 for a 0.1 m step, 0.57 for 1.5 m, converged in dx): the hydrostatic momentum-conserving treatment of an **abrupt** vertical step dissipates energy that the independent reference (energy-conserving, Cd=1) does not. Reference and solver therefore applied to different geometries. Remedy (geometry, not tolerance): entrance **ramps** (run:rise 3) and trapezoidal deposits (taper 0.1 m), the ramps being ice that incises with the lanes (no fixed sill). Result: symmetric-lane head error 0.0119/0.0065/0.0034 m at dx=0.1/0.05/0.025 (first order → reference), Cd 0.88→0.93→0.96; strip partition error 0.057/0.038/0.030, corridor −0.040/−0.034/−0.027, bare −0.026/−0.018/−0.012 against the independent per-column weir reference (gate |err|<0.05 at dx≤0.05: passed). Deposit side-slope sensitivity (taper 0.1 vs 0.2 m): errors similar (`partition_v5.json`), so steep deposit sides do not dominate; hydrostatic validity at the crest (depth/crest-length ≈0.1) is marginal and not independently tested against a nonhydrostatic solver. An early diagnostic sampled the stage 3 cells upstream of the entrance, which moves into the drawdown as dx shrinks (artifact; fixed-distance sampling used in v5b). Two further geometry defects (fixed pool sill, fixed tail sill) that made contact loss impossible were also removed.

## 2. Verification gate (final geometry; `results/verify_swe_v5.json`, `v5b`)
rest/well-balance 4e-15; water and sensible+latent+storage closure ≤1e-12 in every run; plug-flow heat limit error 2.5e-3/1.3e-3/6.4e-4 (first order); bed lowering vs boundary-flux heat 6e-12; outlet-length independence <1e-5. Convergence (development, strip): contact-loss day 3.165/2.721/2.478 d at dx=0.1/0.05/0.025 (observed order 0.86; Richardson ≈2.18 d; GCI₁.₂₅ ±0.37 d at the evaluation grid); bed-update interval 12/6/3 h changes it by 0.007 d; CFL 0.8/0.4/0.2 by 0; sustained reconnection after removal at 8 d at the removal time at all resolutions; cumulative toe melt changes 0.3 %. Unresolved: finest-grid (dx=0.0125) not run (cost ×8); the medium grid is not claimed converged, only bounded by the above U.

## 3. Evaluation campaign (frozen; `results/eval/`, `results/analysis/`)
10 configurations × (medium dx=0.025 evaluation, coarse dx=0.05 numerical check), 21 d, 6-h bed updates. **Actual campaign wall time 6197 s (103 min) with 4 workers; summed run time 5.6 CPU-h medium + 0.7 CPU-h coarse; peak RSS 153 MB per run; no failures.**
Honest status of the configurations: strip none/removal 8 d (and coarse removal 1/4 d) existed in development at coarse resolution, so E04–E06 are not held-out in configuration (only in resolution); E07 (12 d), E03, E08–E10 were not run before the freeze.
| run | status | loss day (θ=1e-3) | rewet episodes | sustained reconnection | f0 → f_end | post-removal mean f | delivered toe-water fraction | toe-reach melt m³ | bypass melt m³ | max water resid | max energy resid | wall s | RSS MB |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| E01_bare_none | complete | — | 0 | — | 0.735 → 0.738 | — | 0.7368 | 0.5912 | 0.5899 | 2.7e-14 | 3.4e-14 | 1580 | 153 |
| E02_corridor_none | complete | — | 0 | — | 0.532 → 0.491 | — | 0.4995 | 0.2988 | 0.5908 | 3.6e-14 | 2.8e-14 | 1460 | 153 |
| E03_corridor_r12 | complete | — | 0 | — | 0.532 → 0.495 | 0.495 | 0.5003 | 0.2991 | 0.5908 | 3.6e-14 | 2.8e-14 | 1459 | 153 |
| E04_strip_none | complete | 2.478 | 0 | — | 0.113 → -0.000 | — | 0.0055 | 0.1609 | 0.5910 | 1.1e-12 | 1.0e-12 | 2028 | 153 |
| E05_strip_r1 | complete | — | 0 | — | 0.113 → 0.192 | 0.609 | 0.5925 | 0.5777 | 0.5904 | 4.4e-13 | 3.8e-13 | 2822 | 153 |
| E06_strip_r8 | complete | 2.478 | 0 | 8.000 | 0.113 → 0.094 | 0.098 | 0.0674 | 0.4419 | 0.5912 | 1.1e-12 | 1.0e-12 | 1907 | 153 |
| E07_strip_r12 | complete | 2.478 | 0 | — | 0.113 → -0.000 | 0.000 | 0.0055 | 0.1621 | 0.5910 | 1.1e-12 | 1.0e-12 | 2022 | 153 |
| E08_porous_none | complete | 10.455 | 0 | — | 0.197 → 0.000 | — | 0.0414 | 0.2247 | 0.5911 | 6.9e-13 | 6.6e-13 | 1930 | 153 |
| E09_porous_r1 | complete | — | 0 | — | 0.197 → 0.189 | 0.608 | 0.5962 | 0.5777 | 0.5904 | 3.3e-14 | 3.6e-13 | 2822 | 153 |
| E10_porous_r12 | complete | 10.455 | 0 | — | 0.197 → 0.000 | 0.000 | 0.0414 | 0.2259 | 0.5911 | 6.7e-13 | 6.5e-13 | 1958 | 153 |

### Comparison with committed reduced forecasts
Numerical uncertainty: contact-loss day ±0.37 d; toe-flow fraction ±0.04
| run | method | loss pred / CFD (d) | loss outcome | reconnection outcome | mean abs f error | toe melt rel err |
|---|---|---|---|---|---|---|
| E01_bare_none | bare_reduced(evolving) | — / — | match | — | 0.011 | +0.16% |
| E01_bare_none | bare_fixed_geometry | — / — | match | — | 0.010 | +0.16% |
| E02_corridor_none | connectivity_reduced | — / — | match | — | 0.025 | -0.91% |
| E02_corridor_none | fractional_cover | 12.41 / — | MISMATCH | — | 0.383 | -38.71% |
| E02_corridor_none | fixed_geometry | — / — | match | — | 0.063 | +3.80% |
| E03_corridor_r12 | connectivity_reduced | — / — | match | True | 0.024 | -1.01% |
| E03_corridor_r12 | fractional_cover | — / — | match | True | 0.271 | +43.90% |
| E03_corridor_r12 | fixed_geometry | — / — | match | True | 0.149 | +44.09% |
| E04_strip_none | connectivity_reduced | 2.49 / 2.48 | match (err +0.02 d) | — | 0.002 | -98.93% |
| E04_strip_none | fractional_cover | 12.41 / 2.48 | match (err +9.93 d) | — | 0.111 | +13.79% |
| E04_strip_none | fixed_geometry | — / 2.48 | MISMATCH | — | 0.431 | +92.72% |
| E05_strip_r1 | connectivity_reduced | — / — | match | True | 0.047 | -51.07% |
| E05_strip_r1 | fractional_cover | — / — | match | True | 0.112 | +0.17% |
| E05_strip_r1 | fixed_geometry | — / — | match | True | 0.151 | +0.17% |
| E06_strip_r8 | connectivity_reduced | 2.49 / 2.48 | match (err +0.02 d) | True | 0.014 | -58.21% |
| E06_strip_r8 | fractional_cover | — / 2.48 | MISMATCH | False | 0.300 | +9.64% |
| E06_strip_r8 | fixed_geometry | — / 2.48 | MISMATCH | False | 0.568 | +9.68% |
| E07_strip_r12 | connectivity_reduced | 2.49 / 2.48 | match (err +0.02 d) | True | 0.002 | -98.94% |
| E07_strip_r12 | fractional_cover | — / 2.48 | MISMATCH | True | 0.223 | +165.52% |
| E07_strip_r12 | fixed_geometry | — / 2.48 | MISMATCH | True | 0.571 | +165.86% |

Classes realised: maintain contact (E01, E02, E03, E05, E09); lose contact (E04, E07, E08, E10); loss then **reconnection after removal** (E06, removal at 8 d: f_toe jumps at removal to ≈0.10, i.e. restored but only ≈1/7 of the bare-route share); **persistent loss after complete removal** (E07 solid strip, E10 porous strip, removal at 12 d: f_toe=0 for the 9 remaining days); **brief rewetting: no episode occurred in any run (class not realised, not invented)**. Early removal (1 d) never loses contact. The porous strip loses contact later (10.5 d vs 2.5 d) but is also locked out after late removal. Contact-loss time depends on the threshold (θ=1e-4…5e-2: 2.71/2.48/2.05/1.08 d for E04), which is why continuous f_toe(t) is reported; all conclusions on reconnection are threshold-independent here.

## 4. Fair comparison with the committed reduced forecasts
Inputs identical (Q, T, widths, floors, volume 0.1 m³, footprint 0.5 m², h_eff, k_debris); methods: connectivity-aware (imported package), fractional-cover (mean-crest, area-weighted shield), fixed geometry. **The thermal closure (series resistance) is the same in CFD and reduced models, so agreement tests hydraulics/connectivity/geometric evolution, not the insulation closure.** The reduced models have no seepage closure: no prediction exists for porous cases (E08–E10).
* Contact-loss time (θ=1e-3): connectivity model error +0.02 d (E04/E06/E07) ≪ CFD uncertainty 0.37 d; fractional-cover +9.9 d (loses contact at 12.4 d); fixed geometry never loses contact.
* Post-removal reconnection outcome (E03, E05, E06, E07): connectivity model correct in 4/4 including the 12-d persistent loss; fractional-cover wrong for E06 (predicts no loss) and right/wrong depending on case; fixed geometry mispredicts E06/E07.
* Continuous partition error (mean |Δf| over 21 d, uncertainty 0.04): connectivity 0.002–0.047 (E05 0.047, just above U_f); fractional-cover 0.11–0.38; fixed 0.01–0.57.
* **Local toe melt: not comparable.** The imported topology model counts melt only on the strip-lane control (its docstring excludes the open reaches), while the CFD toe-reach figure includes open stretches; relative errors of −98 %/−58 % for the connectivity model reflect scope mismatch and must not be read as model failure or success. Fractional-cover melt is within +0.2…14 % for several cases but wrong by up to +166 %.
* Every method is also compared against the **same** numerical uncertainty; where a method's error is below U it is not distinguishable from the CFD.

## 5. Verdict
Within this depth-averaged model and with a common thermal closure, **connectivity and evolving ice geometry improve the specified predictions beyond numerical uncertainty relative to fractional-cover and fixed-geometry approximations**: contact-loss timing, continuous flow partition and the persistent-loss/reconnection outcomes are reproduced by the connectivity-aware reduced model (loss time within 0.02 d) while the fractional-cover model misses loss timing by ~10 d. Equal volume and footprint give opposite outcomes by arrangement: the connected bare corridor keeps ≈0.49 of the flow for 21 d; the transverse strip loses contact at 2.5 d and stays locked out after removal at 12 d. Counter-evidence/limits: removal at 8 d reconnects (the lock-out requires removal after the bypass has incised ≳0.2 m); early removal fully restores access; permeability only delays loss in the drag model. The result does **not** establish that natural debris blocks flow, that a flood removes it, or any behaviour of real boulder seepage, 3-D flow, or whole-cliff melt.

## 6. Implemented / verified / unsupported
Implemented+verified: well-balanced SWE with wetting/drying, conservative heat/water/latent ledgers, plug-flow heat limit, partition vs independent weir reference within 0.04 at dx≤0.05, event-time convergence bounded. Implemented, not independently verified: Ergun-drag porous deposit, debris insulation closure, the nonhydrostatic adequacy near the strip. Unsupported: seepage/underflow, 3-D free surface, debris enthalpy beyond T_m, flood removal, glacier mass loss.

## 7. Next single physical test
A fixed-bed flume measurement of the head–discharge relation (overtopping plus through-flow) of a surveyed boulder strip versus an equal-volume longitudinal corridor with the same lane geometry, to calibrate/replace the drag closure — the only element that decides whether the porous (natural) case behaves like E04 (rapid loss), E08 (delayed loss) or something else.

## 8. Reproduce
`python -m ice_cliff_survival.reactivation.verify_swe OUT.json` (gate); `…convergence DIR OUT.json`; `…freeze OUT_DIR` then `…campaign SPEC_DIR OUT_DIR 4`; `…analyze_eval EVAL SPEC OUT`; `…make_figures_m3 EVAL SPEC OUT` (figures `reactivation/figures/fig1…fig6` PDF+PNG; fig7 regime summary omitted — only 10 runs, no supported regime map).
