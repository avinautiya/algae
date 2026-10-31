# V2 milestone report (verification + baseline prediction pipeline) — updated as results complete

Scope of this report: verified stock-OpenFOAM laminar benchmarks, a separate Stefan benchmark, the baseline `predict` CLI, a smooth-vs-step fixed-geometry comparison (laminar, rigid lid), and the Ogier data audit. **Not done:** free-surface hydrostatic/jump CFD, temperature in a free-surface solver, moving-ice coupling, turbulent runs, production campaign.

## 1. Fixed-geometry laminar verification (stock OpenFOAM)
Software: OpenFOAM ESI **v1912 patch 200626**, Ubuntu package `openfoam 1912.200626-2build3` (apt), gcc 13.3 (binary build), `icoFoam`, `scalarTransportFoam`. Minimal executions demonstrated by these runs (`ice_network_v2/cfd/poiseuille_graetz.py`; logs in the case dirs). `postProcess -func writeCellCentres` crashes in this build (IOstream/sha1 error) — avoided by analytic cell centres on the uniform mesh. No `wmake`/dev headers are installed, so **custom solvers/boundary conditions cannot be compiled yet** (`libopenfoam-dev` not available in the package index).

Case: 2-D plane channel, H = 0.01 m, L = 20 H, ν = 1e-6, Ub = 1e-3 (Re_H = 10), Pr = 7 (Pe_H = 70). Uniform inlet velocity, no-slip walls, outlet p = 0; **Poiseuille check:** velocity profile at x = 0.9 L (18 H; hydrodynamic entrance length ≈ 0.04 Re H ≈ 0.4 H, so developed) against u = 6 Ub η(1−η). **Graetz check:** the same velocity, inlet T = 1 (uniform), wall T = 0 (melting-point reference), outlet zero-gradient, scheme `Gauss linear`, thermally developing flow measured at x/(H Pe_H) = 0.257, where the leading eigenmode dominates. Reference for T: leading Graetz mode **including axial conduction** (eigenproblem with Pe_H = 70; neglecting axial conduction leaves a spurious ~2 % plateau — found and removed; `graetz_exact.py`). Cell Péclet number u·dx/α = 10, 5, 2.6 for ny = 10, 20, 40 (central differencing; still converged smoothly).

| quantity (at x = 18 H) | normalisation | ny = 10 (dy = H/10) | ny = 20 | ny = 40 | observed order |
|---|---|---|---|---|---|
| velocity-profile RMS error | Ub | 9.904e-3 | 2.540e-3 | 6.391e-4 | 1.96, 1.99 |
| temperature-profile shape RMS error (T/T_b) | local flow-weighted bulk excess | 7.130e-3 | 1.807e-3 | 4.582e-4 | 1.98, 1.98 |
| wall heat flux (via Nu_Dh = 2 q_w H/ T_b, 2nd-order wall gradient) relative error | reference value | 3.767e-3 | 4.982e-4 | 7.051e-5 | 2.92, 2.82 |
| Nu_Dh (ref 7.5413) | – | 7.5697 | 7.5451 | 7.5418 | – |
| thermal decay rate error | reference s = 7.5319 | 1.8e-5 | 1.4e-5 | 5.7e-6 | not clean (reference-FD floor) |

Heat flux is verified **separately** from temperature: matching T_b or the T profile does not imply the wall gradient is right; here the wall-flux error is the smallest of the three and converges faster (superconvergence of the quadratic wall extrapolation; reported as observed, not assumed).
Time and iteration errors: both problems are **steady**. Halving/quartering dt (0.1→0.05→0.025 at ny = 20) and tightening the pressure tolerance (relTol 1e-2→1e-4, 4 correctors) changed all error measures by < 1e-13 relative (identical to the printed digits): the converged steady state is independent of dt and linear-solver tolerance (final residuals ≈ 1e-13), so these are negligible against the 1e-3–1e-5 spatial errors. **This does not test transient time accuracy**, which a steady case cannot. Courant number mean 0.20, max 0.30 on the finest mesh. Conservation: outlet/inlet mass flux = 1 − 1e-10.
Limits: Re_H = 10 and Pe_H = 70 are laminar flume-scale values; **thin-thermal-layer resolution (Pe ≫ 1) is not tested here** (the first-mode thermal layer is H/2 thick, resolved by 10–40 cells). Environmental reaches have Pe ~1e5.

## 2. Separate Stefan benchmark (not OpenFOAM)
`cfd/stefan1d.py`: one-phase melting, liquid held at T_w, ice at the melting point; exact Neumann solution; fixed-grid enthalpy method, equal densities. **Ste = 0.2** (nondimensional, chosen to resolve the front on 50–200 cells; with c_w = 4217 J/kg/K and L = 3.34e5 J/kg this corresponds to ΔT ≈ 15.8 K — larger than glacier water; latent heat and c were not altered, but the *environmental* Stefan number ≈ 0.0013–0.0126 is **not yet verified**). Front-position relative error: 9.9e-5 (nx 50), 1.15e-4 (100), 1.8e-5 (200): **not monotone**, so no convergence order is claimed (envelope ≤ 1.2e-4); time-step halving ×2 changes the front by < 1e-8 (explicit scheme at its stability limit). This verifies the latent-heat bookkeeping of the enthalpy method only; the density-jump normal-velocity condition u_w·n = (ρ_i/ρ_w − 1) m_n (MATHEMATICAL_UPGRADES §3) is **not** exercised.

## 3. Baseline prediction pipeline
Command: `python3 -m ice_network_v2.predict ice_network_v2/cases/baseline_synthetic.json OUTDIR [--restart OUTDIR/checkpoint.npz] [--stop-at SECONDS]` (code: `predict.py`, hydraulics/heat step `reach.py` derived from `ice_network_reference/solver.py::channel_step`; reference file untouched). Depth-averaged 1-D shallow-water HLL + upwind temperature + Dittus–Boelter wall heat transfer with **one** interface power used for water cooling and melting (reach.py, q_ice_j); not CFD, no new closure.
Refusals: ice not temperate; missing inlet-temperature definition; later-geometry fields in the case; unknown outlet type. Inlet discharge = **prescribed forcing** (flag `inlet_discharge_is_prescribed_forcing`); width is an output only if `evolve_width` is true (otherwise the interface energy lowers the bed and no width is output — a declared approximation).
Example (`cases/baseline_synthetic.json`: 100 m reach, slope 0.02, width 0.5 m, Q = 0.05 m³/s, T_in = 0.5 °C flux-weighted, free outfall, 4 h after 900 s frozen-geometry spin-up):
bed lowering max 2.37 cm / mean 1.99 cm; melted ice 1.563 m³; heat to ice 4.788e8 J of 1.518e9 J supplied (31.5 %); outlet sensible export 1.055e9 J; delivered 721.41 m³ vs 720.0 m³ supplied (+1.41 m³ meltwater). Outputs: `timeseries.csv` (stage in/out, Q, T_out bulk, P_out, melted ice, bed lowering), `snapshots.npz` (x, t, bed lowering, depth, T, width), `summary.json`, `checkpoint.npz`.
**Residuals** (absolute / relative to named denominator): water 6.4e-12 m³ / 8.9e-15 (inflow volume 720 m³); energy 1.9e-5 J / 1.2e-14 (supplied sensible heat 1.518e9 J); ρ_i·V_melt vs ρ_w·V_meltwater equal to 1e-15; latent energy ρ_i L V_melt vs heat to ice 7.6e-15; geometry change vs latent volume 7.5e-13. These are rounding-level because the discrete scheme is conservative by construction — they verify bookkeeping, not physical correctness.
**Restart:** `--stop-at 7200` then `--restart` reproduces the full run **bit-identically** (snapshots, timeseries, outputs, 135,325 steps). Saved state: depth h, discharge Q, bulk temperature T per cell; bed z and width w (current); time t; cumulative ledgers (V_in, V_out, sensible in/out, heat to ice, meltwater, dissipation, melted ice); applicability-flag accumulators; initial water volume/energy (V0,E0) for the ledgers; record counter, step count, the time-series list and snapshot list. Not saved (re-derived from the case JSON): x grid, original bed z0/width w0, bank height, forcing tables, parameters; the spin-up is not repeated.
Profile (this machine, 4 cores): 21–27 s wall for 4 h simulated (200 cells, 135k steps, incl. numba JIT) ≈ 670 simulated h per wall h (single core); peak RSS 152 MB; checkpoint 0.7 MB, outputs 0.8 MB per case.
Applicability flags in `summary.json` (e.g. Re_Dh range, Fr range, laminar-closure invalid when Re_Dh < 1e4, depth-averaged flag, unsupported-process list).

## 4. Smooth ramp vs single step — fixed geometry, matched budgets (laminar, rigid lid)
`cfd/step_vs_ramp.py`, `step_ramp_post.py`. **Not a hydraulic jump and no free surface**: top is a free-slip adiabatic lid. Re_H = 50 (upstream depth 1 H = 0.01 m, Ub = 5e-3 m/s), Pr = 7 (Pe_H = 350). Same upstream length (12 H), inlet flux and inlet T = 1 (identical inlet enthalpy flux), wall T = 0, downstream depth 2 H, downstream 30 H, same block topology; ramp slope 1:2 (a 1:6 ramp failed checkMesh skewness, 7.5, and was not used; the 1:2 ramp has max non-orthogonality 63°, average 28°). `simpleFoam` (laminar, linearUpwind) then `scalarTransportFoam` to steady state; wall flux from first-cell gradient (first-order).
| mesh (cells) | outlet T_b step / ramp | ramp − step (T_b) | reattachment x/H step / ramp | mean wall flux 0.5<x<10 step / ramp | wall-heat ledger residual (÷P_in) |
|---|---|---|---|---|---|
| n = 8 (4,608) | 0.80305 / 0.80360 | +5.5e-4 | 6.68 / 6.28 | 0.795 / 0.759 | −1.6e-4 / −2.6e-4 |
| n = 16 (18,432) | 0.79800 / 0.79861 | +6.1e-4 | 6.75 / 6.32 | 0.788 / 0.751 | −3.1e-4 / −3.6e-4 |
| n = 32 | pending (running) | | | | |
Reading so far: the ramp lowers the near-field wall flux by ≈ 5 % and shortens the separation by ≈ 0.4 H at both meshes (differences stable under n = 8→16: mean-flux change from mesh 0.007 vs step–ramp difference 0.037), whereas the downstream enthalpy differs by only 6e-4 (0.08 %) — **smaller than the 5e-3 mesh change in T_b**, i.e. not resolved. So at fixed budgets and Re = 50 the geometry redistributes wall heat locally but does not measurably change heat export. Third resolution and a time-independence note (steady) to follow; q uses a first-order gradient so absolute q has O(dy) error.

## 5. Observational route
See `data_audit/DATA_READINESS.md`. Retrieved: the 337,848,233-byte archive (MD5 verified), 454 entries; geometry (5 stations, ±0.5 m) does not support a forward bed-lowering test; discharge forcing exists; temperatures are near-bed point sensors; Nusselt numbers are derived.

## 6. CFD architecture (design only; nothing below is implemented)
* Air–water surface: VOF (`interFoam`/`compressibleInterFoam` exist in this install) with physical surface tension for the thin-flume reference; hydrostatic-rest and conjugate-depth tests still to be run.
* Water enthalpy/temperature: `compressibleInterFoam` carries T for two phases (to be verified for water-only thermal accuracy against the Graetz result above, in the water phase with alpha-weighted properties).
* Solid ice (temperate, base case): no conduction needed: the interface is held at T_m and the Stefan condition is ρ_i L m_n = q_w with q_w the water-side conductive flux; cold ice would need a conducting ice region (`chtMultiRegionFoam`-style coupling).
* Phase change: m_n > 0 recession; ice material velocity ≈ 0; interface (mesh) velocity V_Γ·n = −m_n; water normal velocity u_w·n = (ρ_i/ρ_w − 1) m_n (not the geometric wall velocity); the added water volume ρ_i/ρ_w · (melt) must not be counted twice by ALE and the source.
* Geometry motion: body-fitted ALE with boundary displacement from m_n and a conservative remap; the generic `solidificationMeltingSource` fvOption is a porosity/enthalpy-method source, not a sharp-front free-surface Stefan boundary and is not expected to satisfy the requirements above.
* Blocker: custom boundary condition/solver needs `wmake` and headers (absent); building OpenFOAM from source (or the OpenFOAM.com repository packages) is the next toolchain step.

## 7. Unresolved issues
1. n = 32 step/ramp meshes and hence a third-resolution envelope; wall-flux gradient is first order.
2. Free-surface hydrostatic rest, conjugate-depth jump, VOF + temperature coupling, moving-ice coupling: not started (need the toolchain note above).
3. Environmental-Stefan (Ste ~ 1e-2) benchmark and the density-jump normal velocity not exercised.
4. Stefan benchmark and the baseline use different codes (Python) from the CFD (OpenFOAM); the reach model's closure (Dittus–Boelter at Re ~1e5) is untested by CFD.
5. Observation operator and sensor/time audit for Ogier data not built.
6. No production campaign launched; no novelty claim made.
