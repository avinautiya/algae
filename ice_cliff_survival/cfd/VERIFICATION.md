# 2-D heat + moving-ice solver (`stefan2DFoam`): equations, verification, profile

**What it is.** A finite-volume (OpenFOAM v1912 `finiteVolume`, custom app built per `ice_network_v2/foam/BUILD.md`) conduction solver for a liquid layer 0 < x < s(y,t) bounded by a moving ice interface (ALE: structured columns stretched in x), hot wall x = 0 with prescribed T_w(y), adiabatic sides, temperate ice (T_m = 0) beyond the interface. Equation `d/dt∫T dV + ∫T(u_w − w)·dS − ∫α∇T·dS = 0`; interface law ρ_i L m_n = q with **q the same discrete flux the temperature equation removes** (face `snGrad`), interface a graph x = s(y,t) advanced by face-averaged swept volume V_f = q A_f dt/(ρ_i L) (Heun predictor–corrector); meltwater enters at T_m (**no enthalpy flux through the ice patch** — found necessary by the ledger, below). Equal-density limit u_w = 0; with ρ_i ≠ ρ_w a uniform liquid velocity (1−ρ_i/ρ_w)⟨ṡ⟩ x̂ supplied through the hot wall closes the global liquid-mass balance (local lateral liquid flow is **not** represented, so with non-uniform ṡ this is an approximation). **Static debris shielding:** per-row conductance h_eff = k_deb/d (mixed BC, value fraction h/(h+kδ); debris heat capacity ignored; `geometry.shield_conductance` maps a static wedge to h_eff). **Not included:** flow, free surface, grains (mobile or otherwise), cold ice, atmosphere, collapse. Source `cfd/stefan2DFoam/stefan2DFoam.C`, driver `cfd/stefan2d.py`, suite `cfd/verify_cfd.py`, stored results `cfd/results/verify_cfd.json`; unit smoke test `tests/test_cfd_smoke.py` (runs if the app is built).

## 1-D limit (ny = 4 identical rows; exact Neumann solution with density jump as reference)
| case | front rel. error, nx = 20/40/80/160 (dt 2.5 s) | observed order | dt = 2.5/5/10 s (nx 40) | row spread |
|---|---|---|---|---|
| Ste 0.00631, ρ_i = 917 + jump | 3.0e-5 / 1.5e-5 / 7.6e-6 / 3.7e-6 | 0.98, 1.01, 1.04 | 1.53e-5 / 1.52e-5 / 1.55e-5 | 0 (exactly) |
| Ste 0.2, ρ_i = 917 + jump | 1.0e-3 / 5.2e-4 / 2.5e-4 / 1.1e-4 | 1.00, 1.05, 1.14 | 5.2e-4 / 5.0e-4 / 4.8e-4 | 0 |
| Ste 0.00631 and 0.2, equal density | same pattern (3.2e-5…4.0e-6; 1.1e-3…1.2e-4) | 0.99–1.13 | insensitive | 0 |
Reference uncertainty: root of the similarity equation to 1e-15 (negligible). The 1-D results reproduce those of the independent `stefan1DFoam` (e.g. 1.53e-5 vs 1.47e-5 at nx = 40, Ste 0.00631).
## Ledgers (two kept separate)
* **Discrete energy ledger** (storage change + interface conductive loss = hot-boundary supply, all from the numerical fluxes): |residual|/supply ≤ 8e-14 in all 1-D runs, ≤ 1.2e-13 in 2-D, ≤ 2.4e-13 with shielding.
* **Physical sensible-plus-latent balance** (latent ρ_i L V_swept from the geometry vs supply − storage): residual/latent = 4.2e-4, 8.3e-4, 1.7e-3 at dt = 2.5, 5, 10 s (Ste 0.00631; **first order in dt**, mesh-independent); at Ste 0.2: 2.4e-3, 4.9e-3, 9.7e-3. It is the Euler-vs-trapezoid inconsistency between the interface-advance and the flux integration, reported and not hidden.
* **Defect caught by the ledger:** the first version let the ALE flux through the ice patch carry the face enthalpy, so a shielded face (T_f ≠ 0) created liquid at T_f; the discrete residual was 3e-4…6e-3 (∝ Ste) with shielding and 1e-14 without. Zeroing the enthalpy flux through the ice patch (meltwater enters at T_m) restored ≤ 2.4e-13.
## Interface mass balance and remapping
* Swept volume from the **OpenFOAM mesh flux** equals the geometric volume change (geometric conservation law): relative residual ≤ 6e-15 (1-D) and ≤ 1e-13 (2-D).
* Energy-side swept volume Σ qA dt/(ρ_i L) / geometric = 0.9992 (jump) and 0.9991 (equal density): the 8e-4 offset is the first-order-in-dt physical residual above.
* Liquid mass: ρ_w ΔV_L = ρ_w V_in + ρ_i V_melt closes to 3.5e-14 and V_in/V_sw = 1 − ρ_i/ρ_w = 0.0830 to 12 digits — **this is an identity by construction** of the uniform inflow, not an independent test; the independent content is the GCL and the energy-side comparison above. Equal-density limit V_in = 0 exactly.
* Remapping: the ALE scheme has no re-meshing, so geometric conservation is the in-run remap statement (above). A separate conservative remap of the solved temperature field to a different column count (`cfd/remap.py`, exact overlap weights in ξ-space because s is linear in y within a row) preserves the sensible energy to ≤ 2.0e-16 relative (nx 20 → 40, 80, 7), whereas point-sampling interpolation changes it by 2e-7 … 5e-6; accuracy is first order and bounded by the source mesh (0.3–2 % deviation from the local 1-D profile on the 20-cell source). In-run remeshing is **not implemented**.
## 2-D mesh/time convergence (T_w(y) = T_0(1 + 0.3 cos(πy/L_y)), L_y = 5 cm, Ste 0.00631, ρ_i = 917; independent series)
| grid (nx×ny), dt = 2.5 s | swept volume V_sw [m³] | difference to previous |
|---|---|---|
| 10×8 | 8.705991e-8 | |
| 20×16 | 8.705351e-8 | 6.4e-12 (7.4e-5 rel) |
| 40×32 | 8.705086e-8 | 2.65e-12 (3.0e-5) |
| 80×64 | 8.704965e-8 | 1.21e-12 (1.4e-5) |
Consecutive-difference ratios 2.4, 2.2 → first-order in the mesh size; Richardson estimate of the remaining error at 80×64 ≈ 1.3e-12 (1.5e-5 relative). Time-step series at 40×32 (dt = 2.5, 5, 10, 20 s): V_sw = 8.705086, 8.705085, 8.705090, 8.705124 (e-8): spread 3.8e-13 (4e-6 relative) — **temporal error is below the spatial error** at every dt used. The local front position at mid-height converges 0.0027528 → 0.0027532 → 0.0027532 m. Sanity (not verification): the 2-D interface differs from the row-wise 1-D Neumann prediction (lateral conduction ignored) by ≤ 2e-6 m (7e-4 relative).
## Static debris shielding
Uniform shield, Ste → 0 quasi-steady analytic reference s/h + s²/(2k) = s_0/h + s_0²/(2k) + T_w t/(ρL): deviation −1.6e-5 (h = 30, Ste 5e-4) … −1.3e-3 (h = 500, Ste 0.00631); it scales with Ste (the analytic ignores sensible storage) and is mesh-independent over nx = 20–80. Limits: h = 10¹² reproduces the bare run to 1.5e-10; h = 10⁻⁶ advances the interface 5e-12 m (bare: 1.65e-3 m). Partially covered face (upper half h = 30 W m⁻² K⁻¹): covered/exposed advance ratio 0.0949 at three meshes (spread 1e-4), ledgers 1e-14.
## Profile (single core, this container)
| grid | cells | steps (dt 5 s, 3600 s) | wall | µs per cell-step | peak RSS | case disk |
|---|---|---|---|---|---|---|
| 40×32 | 1,280 | 720 | 5.9 s | 6.4 | 72 MB | 1.0 MB |
| 80×64 | 5,120 | 720 | 24 s | 6.5 | 74 MB | 2.1 MB |
| 160×128 | 20,480 | 720 | 123 s | 8.3 | 74 MB | 6.7 MB |
| 80×64, dt 2.5 s | 5,120 | 1,440 | 52 s | 7.1 | 75 MB | 2.6 MB |
Cost model: ≈ 7×10⁻⁶ s per cell-step (4 temperature solves per step). Restart/checkpointing: not implemented in the solver (cases take < 5 min); the driver is idempotent per case directory.
## Limits
Rigid fixed liquid layer, equal or uniformly fed liquid, one interface graph, static shield without heat capacity, temperate ice. The interface advance is first order in space; time error of the physical ledger first order. Nothing here verifies flow coupling, free-surface effects, or grain transport.
