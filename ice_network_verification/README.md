# Independent verification for the ice-network project

This package supplies analytical limiting cases and an export checker independently
of the production `ice_network` solver. It does **not** demonstrate novelty, field
forecast skill, a valid open-channel heat-transfer closure, or CFD validation.
No production solver files are modified. Core functionality uses Python 3.9+ and
the standard library. The optional production-kernel audit additionally uses NumPy.

## Run

From the repository root:

```sh
python -m unittest discover -s ice_network_verification/tests -v
python -m ice_network_verification generate /tmp/ice-analytical
python -m ice_network_verification evaluate /tmp/solver-exports --report /tmp/ice-check.json
```

`generate` writes parameters, units, requested times, reference CSVs and SHA-256
hashes. Reference filenames end in `.reference.csv`. A solver adapter must create
`<case_name>.csv` with its own numerical results at the requested times. Simply
renaming analytical files is a harness identity check, not solver verification.
`evaluate` returns exit code 1 for missing cases, failed comparisons or malformed
exports. It will not silently omit an unimplemented case.

The committed CSVs in `records/analytical_references/` are **derived references,
not measured observations or production-model results**.

## Seven limiting cases and assumptions

| Case | Independently derived reference | What it checks |
|---|---|---|
| `linear_reservoir` | `h=h0 exp(-G t/A)`, `Q=G h` | Storage/outflow integration |
| `finite_drainage` | `h=max(0,sqrt(h0)-C t/(2 A))^2`, `Q=C sqrt(h)` | Emptying and nonnegative storage |
| `symmetric_outlets` | `h=h0 exp(-(G1+G2)t/A)` | Equal branch splitting |
| `unequal_outlets` | Same head; outlet volumes in ratio `G1:G2=1:3` | Junction accounting beyond symmetry |
| `pre_overtopping` | `h=h0+I t/A` up to `t=A(h_s-h0)/I` | First sill contact and no premature outflow |
| `thermal_front` | `T=Tin exp(-beta x/u)` after `t=x/u`; zero before | Thermal transport timing and cooling |
| `heat_to_melt` | `Tout=Tin exp(-UA/(rho cp Q))` | Sensible heat, latent heat, and density conversion |

Linear conductance and ideal reservoir outlets are deliberately prescribed
constitutive laws. They do not test the production Manning/GVF closure. Use a
test-only edge-law adapter if needed; do not change production physics merely to
match these cases. If the solver API cannot represent them, report that explicitly.
They are useful unit-level references, not seven configurations of a real glacier.

The temperature-front case has zero diffusion, prescribed uniform velocity and a
linear sink, not a complete melting channel. Its value at the exact front is the
post-front limit. Numerical diffusion and the discontinuity must be considered in
a refinement study; a failed pointwise gate does not alone diagnose an incorrect
scheme. The heat-to-melt case has fixed geometry, a constant effective conductance,
temperate ice, no channel water storage, no heat source, and no refreezing. It
verifies the integrated heat budget, not evolving geometry. Here cp is explicitly
4180 J/kg/K; match manifest constants rather than silently using another value.

For each case, exported scalar fields and units are specified by `manifest.json`.
Times must match within 1e-9 s and be strictly increasing. No interpolating or
resampling is done by the checker. Candidate quantities are nonnegative in these
particular cases. NaNs, infinities, missing outputs, duplicate times and negative
quantities fail closed.

## Error metrics and numerical gates

The checker reports RMSE, mean signed error, maximum absolute error, and RMSE
normalized by the maximum reference magnitude. Normalized RMSE is null for an
identically zero reference. Default pointwise gate:

`abs(error) <= 1e-8 + 0.01 * abs(reference)`.

The absolute floor is in **each field's SI units**. CLI `--atol` and `--rtol`
override it. These are explicit numerical smoke-test tolerances, not observational
uncertainty or a scientific success threshold. Use tighter, field-specific checks
in a solver adapter where appropriate. Passing this gate is not convergence:
compare errors on several resolutions against the same independent reference.
`checks.observed_order` calculates pairwise reference-error orders for a specified
uniform refinement ratio. It does not extrapolate a numerical truth.

## Optional global ledgers

Export all columns of a ledger to request its audit. A partial ledger fails.
Without a ledger the report says `not_submitted`; passing solution fields does
not imply a budget audit occurred. Cumulative ledgers are differenced from the
first exported row; export initialization at time zero when auditing a full run.

Mass columns (kg):

`liquid_mass_kg, mass_in_kg, mass_out_kg, melt_mass_kg`

The checked residual is:

`delta liquid - delta in + delta out - delta melted`.

Liquid storage is the whole modeled domain; internal edge fluxes are not external
in/out. `melt_mass_kg` is the ice mass converted into liquid, not lost ice volume.
No refreezing is supported by this ledger contract. The checker does not infer
correctness of an independently exported ice-geometry inventory; compare it with
`rho_i * ice_volume_lost_m3` separately.

Thermal-energy columns (J):

`sensible_energy_j, latent_energy_j, boundary_heat_j, advected_energy_in_j,
advected_energy_out_j, dissipated_to_heat_j`

The checked residual is:

`delta sensible + delta latent - delta boundary_heat - delta advected_in
 + delta advected_out - delta dissipated_to_heat`.

Sensible heat uses the melting-point datum consistently. Boundary heat is signed
net conductive/radiative/etc. heat entering the domain, **excluding** heat already
counted as advection or mechanical dissipation. This is a thermal inventory, not a
total mechanical-energy test. The checker cannot detect an origin counted twice
if every exported field is adjusted to agree; independently verify source terms.

Reports include absolute residuals, relative residuals and whether inventory change
is nonzero. A zero-flow identity can pass, but is labeled as such and cannot
substitute for nonzero-throughput checks.

## Independent production-kernel audit

`audit_heat_kernel.py` inspects the two actual production thermal-function bodies
and compares them with an independently evaluated 60-digit solution of the
constant-coefficient segment ODE. This is a numerical audit of a frozen law,
**not independent evidence that the law represents open-channel physics**.

```sh
python -m ice_network_verification.audit_heat_kernel /path/to/solver-checkout --report /tmp/kernel.json
```

The input checkout must contain
`ice_network/src/ice_network/{constants.py,thermal/channel_heat.py}`. The audit
records file hashes. It executes the original function bodies without Numba
decorators; therefore it does not verify JIT behavior or production performance.
Review downloaded code before executing it. Only recognized function definitions
and literal constants are loaded, not the original package initialization.

See `CLAUDE_HANDOFF.md` for the pinned audit findings and integration instructions.
