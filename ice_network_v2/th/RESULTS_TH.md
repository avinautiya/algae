# Matched-state upstream-thermal-history experiment — results (rigid-lid laminar, fixed geometry)

Protocol: `FREEZE_TH.md` (committed before any result). Solver: stock OpenFOAM v1912 `simpleFoam` + `scalarTransportFoam`; Re_H = 50, Pr = 7, Pe_H = 350; free-slip adiabatic lid (declared rigid-lid proxy); adiabatic upstream wall segments (inlet/wall corner removed from the measurement region; the only temperature-discontinuity corners sit 3 H (A) and 11 H (B) upstream of the matching plane x_c = −1 H). Code: `th.py`, `baselines.py`, `analyze.py`, `conductance_test.py`, `figures_th.py`; data: `results/*.json`, `analysis.json`, `conductance_test.json`; figures `figures/th1…th5`.

## 1. Matching (flow identical by construction; temperature scaled)
Both histories use the same simpleFoam flow field, so discharge, depth and the velocity profile at x_c are identical to round-off. History B is rescaled (inlet 1.078×) so the flow-weighted T_b at x_c equals A's (n=32: T_b,A = 0.96307, T_b,B = 0.89328 for unit inlet; Graetz series 0.96337 / 0.89420, scale 1.0773 vs 1.0781 resolved). Remaining mismatch = the tested quantity: T-profile shape (th1). Mismatch bounds: moving the matching plane to −0.5 H / −2 H changes δ_10 by +0.2 / −0.6 pp; (1 ± 10⁻³) scaling changes δ by ≤ 10⁻³.
Resolved upstream development vs classical Graetz series (film profile, no axial conduction): T_b −0.03 % (A, 2 H of cooling) and −0.10 % (B, 10 H); wall flux +0.26 % / −0.28 %.

## 2. Observables (n = 32 values, step vs straight)
| observable | step | straight | step − straight |
|---|---|---|---|
| δ_10H = (Q_A − Q_B*)/mean, bed 0–10 H | 0.2192 | 0.1971 | +0.022 |
| δ_30H, bed 0–30 H | 0.1350 | 0.1044 | +0.031 |
| spatial error e_q (0–30 H) | 0.1350 | 0.1044 | (q_A > q_B* at every x, so e_q = δ_30) |
Meaning: matched upstream bulk state but a short (A) versus long (B) thermal development changes the integrated heat delivered to the bed by 22 % (0–10 H) and 13.5 % (0–30 H) in the stepped geometry, and by 20 % / 10 % in the straight channel.

## 3. Error budget (computed; components in `analysis.json`)
| component (δ_30 step) | value |
|---|---|
| mesh: |n32 − n16| | 0.0064 (sequence 0.1423, 0.1413, 0.1350 for n = 8, 16, 32: **not monotone-convergent**; full spread 0.0073) |
| time step (dt 5/10/20 s) and final-window iteration change | 3e-4 / 0 (steady problem; dt only sets the pseudo-time path) |
| matching plane / scaling | 0.0082 |
| outlet length (28 vs 40 H), upstream length (16 vs 24 H) | 0.0054 |
| history-B variants (cold start −9/−11/−13 H): **physical** dependence, included by the frozen rule | 0.048 |
| **U (frozen rule)** | 0.049 (δ_10: 0.072); straight 0.037 (0.064) |
| U without the B-variant spread (numerical only) | 0.012 (δ_10: 0.012); straight 0.010 (0.008) |
The frozen rule counts the B-start variants in U; they are a physical sensitivity (the effect depends on how long the development was), not a numerical error, which makes the frozen U conservative. The step mesh sequence is not asymptotic, so |n32 − n16| understates the mesh error; the frozen rule was followed and the caveat is recorded.

## 4. Baselines
* **L (local equilibrium, history-blind, fully developed Nu_H = 1.885 for the film profile):** predicts δ = 0 by construction; error = the full 22 % / 13.5 %. Its own absolute error is as large: Q_L − Q_A = −19.8 % (0–10 H) and −26 % (0–30 H) for the step, −35 % / −22 % straight (it also omits the step face and separation); against B* the errors are −0.04 % / −15 % (step). The history effect is not the dominant error of L in this regime — the fully developed closure is wrong by 13–35 % either way.
* **D (classical developing-flow, history-aware, geometry-blind = straight resolved δ):** predicts δ_10 = 0.197, δ_30 = 0.104 for the step; error +0.022 / +0.031 absolute (10 % / 23 % of the effect), i.e. it reproduces ≥ 77 % of the history effect.
* Cost (n = 32): step 98k cells, flow 374 s + 76 s per history; straight 57k cells, 144 s + 40 s; L and the Graetz series: milliseconds.

## 5. Decision (frozen rule) and sensitivity
δ_X exceeds U (0.219 vs 0.072; 0.135 vs 0.049) so the history effect is **not** negligible; the step-specific excess (+0.022, +0.031) is within U_excess (0.097, 0.062): **decision (2): ordinary developing-flow history effect, reproduced by the straight-channel baseline D; no step-specific claim; no thermal-state model was developed and no withheld-history test was run (not triggered).**
Sensitivity (reported, not used for the decision): with numerical-only budgets U_excess = 0.014 (10 H) and 0.015 (30 H), the excesses (0.022, 0.031) would be 1.5× and 2.0× larger than that — the frozen conservative accounting hides a small (2–3 percentage-point) step-specific modification of the history effect. Whether it is real cannot be settled with the current non-asymptotic step mesh sequence.

## 6. Forecast-package review (`ice_incision_forecast`, commit 471415c; 22 stdlib tests pass)
Theorem 1 (exact map P = P0 e^(−s), meltwater-corrected) reduces to the passive CFD problem as Ste → 0: s(x) = ln(T_b(0)/T_b(x)) = G(x)/C0. Theorem 3 (forcing-uniform bounds) needs a(x,t) = γ(t) a0(x). `conductance_test.py` tests this with the CFD: the per-3H conductance increments of B* over A vary from 0.72 to 0.92 (step n=32; spread 19 %, straight 25 %); a common-shape forecast (A's shape rescaled to B's total) mis-allocates regional heat by +8–10 % in the first 9 H and −4 to −7 % in the last 12 H of 0–30 H. **The common-conductance-shape assumption is violated by changing thermal history** (also in the straight channel), so the theorem's forcing-uniform statement is not applicable as stated; a per-cell conductance band of ±14–16 % around the mean ratio would be needed (the theorem's *box-class* bounds remain valid if that wider class is used, but that is an unverified inference, and the sharp bounds then widen accordingly). The step n = 16 vs 32 tails differ (last bins 0.92 vs 0.85), another sign of unconverged shape.

## 7. Supported claim and limits
**Supported:** in a laminar rigid-lid channel at Re_H = 50, Pe_H = 350, two physically attainable upstream thermal developments matched in discharge, depth, velocity profile and flow-weighted bulk temperature differ in integrated bed heat over the 10 H after a 1 H bed step by 22 ± ≤7 % (frozen U; ±1.2 % numerical-only) — and over 30 H by 13.5 ± ≤5 % — a local fully developed closure cannot represent it, and a straight-channel developing-flow control reproduces ≥ 77 % of it. Known Graetz-type physics; **no new mechanism is supported.** A step-specific modification of 2–3 pp is unresolved.
**Not supported / limits:** laminar, rigid-lid, 2-D, steady, fixed geometry, one Re/Pe; nothing about free-surface flow, turbulence, melting feedback, or environmental heat transfer (Re_h 1e4–1e5, Pe ~1e5); only two (plus two sensitivity) histories; no withheld forecast test; mesh sequence for the step not asymptotic; thermal-layer resolution at Pe_H = 350 is 1–3 cells per Lévêque thickness at n = 8–16.

## 8. Itemized runtime estimate for the next milestone (from measured times; 4 cores)
| task | cells | est. wall (one core) |
|---|---|---|
| step n = 64 flow (to test asymptotic mesh sequence) | 393k | 374 s × ~4–8 ≈ 25–50 min |
| thermal solves, step n = 64, 4 histories (A, B, B9, B13) | 393k | 76 s × ~4 × 4 ≈ 20 min |
| straight n = 64 flow + 4 histories | 229k | ≈ 10 min + 10 min |
| outlet/inlet/dt variants at n = 32 (step) | 98k | ≈ 4 × 8 min ≈ 30 min |
| (only if rule (3) triggers) withheld histories: 6 start positions × 2 geometries at n = 32 | 98k/57k | 12 × ~1.5 min ≈ 20 min after flows (flows reused) |
Total ≈ 2–3 CPU-hours, ≈ 1 h on 4 cores. A free-surface or turbulent version is not affordable (≈ 20 CPU-h per 1e5-cell case, see MILESTONE2_REPORT.md).
