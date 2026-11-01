# Corrections to milestone 1 (a680fae) made in milestone 2

Nothing from M1 was overwritten: `REPORT_M1.md`, `cfd/VERIFICATION.md`, `cfd/results/verify_cfd.json` and the M1 commit stay as published. Each defect is fixed in code, with regression tests in `tests/test_m2_fixes.py`.

## 1. `geometry.debris_normal_thickness` (wrong formula)
M1 returned (horizontal width) × sin θ = (z_f−z)(cot β − cot θ) sin θ. The normal-to-face distance from the face point at height z to the pile surface (ray–line intersection; face through (−z cot θ, z), normal (sin θ, cos θ), pile surface through the face point at z_f descending at angle β) is
**d = (z_f − z) tan(θ − β) / sin θ**, with 0 < β < θ ≤ 90°. Derivation: P + d n on Z = z_f − tan β (X − x_f) gives d (cos θ + tan β sin θ) = (z_f − z)(1 − tan β cot θ), i.e. d = (z_f−z) sin(θ−β) / (sin θ cos(θ−β)). Test: `ray_thickness()` in the tests solves the linear intersection directly (no identities), for several inclinations, θ=90°, crest and invalid input.
Repro θ=62°, β=33°, z_f=0.012, z=0.003: M1 0.00801133650308 m → corrected 0.00565014447451 m (kept as `debris_normal_thickness_v1_projected_horizontal`, deprecated).
Effect on the only M1 example that used it (`heat-ice-2d --shield wedge`, rerun by `cfd/wedge_correction.py` → `cfd/results/wedge_correction.json`): row-0 thickness 10.07 mm → 7.10 mm; mean interface advance after 1 h 2.384 mm → 2.425 mm (+1.7 %). The M1 *verification* of shielding used uniform conductances and is unaffected.
Scope: h = k/d remains a 1-D normal effective-resistance approximation (massless debris, normal conduction only); the uniform-conductance tests do not validate any wedge geometry.

## 2. `validate_records` accepted NaN/inf
`ExchangeRecord(0,1,nan,1,0,0,0,0,0)` passed because NaN comparisons are false. All numeric fields are now checked finite; tolerances are validated (finite, ≥0). `diagnose_records()` returns residuals and never applies an energy/latent tolerance; `validate_records()` is the acceptance gate that applies them. Signed `storage_change_J`, `sensible_*` and `dissipation_J` are kept; negative `q_wall_W` remains representable diagnostically but fails the gate (refreezing unsupported). Latent residual is normalised by max(|heat|,|latent|).

## 3. `neumann_fixture` mixed densities
The equal-density, zero-liquid-flow similarity solution was paired with ρ_i=917, ρ_w=1000, so emitted heat was 0.917 × its own conductive heat. Now one density is used, ρ_i≠ρ_w arguments raise, and `ExchangeRecord` carries the producing solver's densities. Heat is the exact integral of the analytic front flux; tests integrate that flux independently (quadrature) and compare with the records and with ρ L ΔV. No ledger is rescaled.
`cfd/stefan2d.solve_lambda` was checked, not assumed: `TestSolveLambda` builds T(x,t) from the returned (λ,b), and verifies the advection–diffusion PDE with u=(1−r)ṡ, both boundary values, and the Stefan condition by finite differences for Ste ∈ {0.01,0.2,1}, r ∈ {1, 0.917}. It is **not** defective (PDE residual ≈2e-5 relative; Stefan condition ≤1e-3 limited by the FD stencil).

## 4. `ToeRun` dt / horizon
dt must be finite, >0 (and not a bool); constructor raises before any loop. `advance()` rejects non-finite horizons and backward horizons; an unchanged horizon is a no-op.

## 5. Restart identity bound to executable definitions
`closures.py`: serialisable registry closures (`ConstantSources`, `ConstantCapacity`, `LinearCapacity`, `TableSources`). Checkpoint identity = version, config fingerprint, dt, label, provenance, and per closure: type + parameters + (for tables) the hashed forcing content up to the checkpoint time. Forcing extension policy: a table may be extended after the checkpoint time (intervals starting before it must be identical, `extension='forbid'` past the table end). Opaque callables need `unsafe_opaque=True`; their checkpoints need `unsafe_restart=True`, are stamped `restart_verified=False`, and cannot be adopted by a registry run. Function names/repr/ID strings are not used (`forcing_id` removed from the API; callers of the old signature must migrate). Malformed or incomplete checkpoints raise.
