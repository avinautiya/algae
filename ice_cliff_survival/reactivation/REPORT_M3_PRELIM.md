# M3 preliminary report — bypass/thermal-memory experiment: audit, solver, verification status, cost (NO final dataset yet)

Imported package: `imported/ice_stream_reactivation_3e51620/` (pinned commit 3e51620d2194430fedc4484b9c08df8ced5ad9d9, not merged). No evaluation run has been launched; no production freeze exists.

## 1. Audit — what actually exists
| capability | status |
|---|---|
| free-surface flow with Stefan ice motion (3-D/VOF) | **does not exist**: `interThermalFoam` = interFoam + passive heat, no phase change/ice motion; `flowStefanFoam` (M2) is single-phase rigid-lid. Building+verifying VOF+Stefan+moving mesh was not feasible here. |
| supported alternative implemented | `swe2d.py`: plan-view depth-averaged shallow-water solver (hydrostatic reconstruction, wetting/drying, HLL, Manning, optional Ergun drag), conservative depth-weighted temperature, bed-exchange sink h_eff·T that is the same flux lowering the ice floor (ρ_i L dz/dt = h_eff T). **Not 3-D CFD; no 3-D/under-rock flow; porous deposit = drag closure.** |
| energy/mass ledgers | water mass and sensible-heat storage/inflow/outflow/latent closed to ≤1e-13 (with storage terms); latent vs boundary-flux heat removed: 4e-12 |
| removal ledger | exported debris volume logged; debris assumed at T_m ⇒ zero exported enthalpy relative to T_m (declared) |

## 2. Defects found and fixed while building (all in-code, tested by the gate)
lateral ghost-layer leak of water at deposits (9 % mass loss; fixed); zero-gradient outlet on a flat deep tail produced backwater (replaced by a free overfall: 1.5 m drop + slope 0.5); a fixed pool sill and a fixed tail sill both capped incision (no contact loss possible) — replaced by a deep reservoir and an overfall so that the only controls are the incising ice floors; entrance step caused jump losses (kept as a physical step-up weir).

## 3. Verification gate — status (data `results/verify_swe_v4.json`; tolerances were declared in `verify_swe.py` before judging)
| check | result | verdict |
|---|---|---|
| lake at rest with partly dry deposit | max momentum <1e-10, stage deviation <1e-12 | pass |
| discharge/water/heat conservation incl. storage | ≤1e-13 | pass |
| plug-flow heat-depletion limit (analytic exp(−NTU)), lowering profile | error 2.7e-3 → 1.4e-3 → 6.8e-4 for dx = 0.1/0.05/0.025 (first order) | pass |
| bed lowering vs heat removed from water (independent boundary-flux measurement) | 4e-12 | pass |
| outlet-length independence (0.5/1.0/1.5 m tail) | f_toe changes <1e-5 (earlier geometry) | pass (rerun with final tail not repeated) |
| flow partition vs independent weir / friction-GVF reference, tol 0.05 abs | bare −0.073/−0.068/−0.065; corridor −0.062/−0.046/−0.042; **strip +0.114/+0.110/+0.107** | **FAIL for strip and bare (>0.05), slow convergence; corridor passes at dx ≤ 0.05** |
| separate spatial and temporal convergence of event timing | **not done** (needs the dev set below) | open |
| Stefan limit | analytic lowering-profile check only (quasi-static ice, no conduction); the M2 OpenFOAM Stefan verification does not apply to this solver | partial |
The strip/bare partition discrepancy is unexplained (candidates: entrance contraction/loss not in the 1-D references, residual hydrostatic-pressure error at the shallow strip crest). **Until it is resolved or bounded, nothing produced by this solver may be called final evaluation data.**

## 4. Development runs (coarse dx=0.05, 6-h bed updates, 15 d, no removal; `results/dev/`) — exploratory, NOT evaluation
| layout | f_toe d0 → d15 | contact loss (f_toe < 1e-3) | toe/bypass melt (m³) |
|---|---|---|---|
| bare | 0.678 → 0.685 | none | 0.423 / 0.423 |
| corridor (bare corridor, 0.25 m of 0.5 m lane) | 0.479 → 0.465 | none | 0.213 / 0.423 |
| solid transverse strip (equal area & volume) | 0.249 → 0.000 | day 5.6 | 0.151 / 0.423 |
| porous strip (Ergun drag, φ=0.4, D=0.2 m) | 0.239 → 0.000 | day 12.8 | 0.209 / 0.423 |
So in this model an unbroken bare corridor preserves access while the same debris arranged transversely loses it (solid strip: 5.6 d, porous: 12.8 d, i.e. permeability delays but does not prevent it). Contact loss is tied to the bypass lowering while the covered toe floor melts at ≈half rate. Consequences after removal have NOT yet been run. Because the partition gate is not met, these timings carry an unquantified error.

## 5. Measured cost (idle machine, 4 cores, `results/benchmark_v3.json`)
| resolution | cells | wall per simulated hour | spin-up | peak RSS | timestep limit (CFL 0.4) | hydraulic steps per 6-h update |
|---|---|---|---|---|---|---|
| coarse dx = 0.05 m | 1,680 | 0.26 s (6.3 s/day) | 7 s | 151 MB | 2.1–2.5 ms | 9.6 k |
| medium dx = 0.025 m | 6,560 | 2.01 s (48 s/day) | 50 s | 152 MB | 1.07–1.26 ms | 18.1 k |
A 15-day run: coarse ≈ 1.7 min, medium ≈ 13 min; finest dx = 0.0125 (≈26 k cells) ≈ 8× medium (≈100 min per 15 d, only affordable for short windows).
Bed-update interval: 6 h (cost ∝ 1/Δt_m); each update marches ≥ 24 s of hydraulic time (dominant cost).

## 6. Proposed smallest campaign (NOT launched; decision requires the partition issue to be resolved first)
Configurations: bare ×1, and {corridor, solid strip, porous strip} × {no removal, early, late} = 10 runs; early/late removal at 0.5× and 1.5× the strip no-removal loss time (declared from the strip development timescale, exploratory). Horizon 15 d.
Cost: 10 medium runs × 13 min = 130 CPU-min; coarse repeats 17 min; bed-update-interval series (12 h/3 h, coarse, 3 runs) 15 min; finest-grid 2-day checks (3 layouts) ≈ 40 min; **total ≈ 3.4 CPU-h ≈ 55 min wall on 4 cores (4 × 152 MB)**. If the medium–coarse timing difference is within the declared uncertainty, medium alone is the converged set; otherwise a dx=0.0125 subset is required (+≈ 5 CPU-h).

## 7. Not done / unsupported
Convergence series for event timing, matched removal runs, production freeze, evaluation dataset, comparison with fixed-geometry/fractional-cover/reduced-model forecasts, figures 1–7. Unsupported physically: seepage/under-rock flow beyond a drag closure, 3-D free surface, debris enthalpy export, flood removal.
