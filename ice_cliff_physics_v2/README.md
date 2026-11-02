# Conservative water–debris–ice cliff physics

This additive research package upgrades the earlier reduced cliff forecast. It is an executable model with independent analytic and numerical checks, not a certification of field superiority, historical originality, journal suitability or competition outcomes.

## Concrete improvements

- Finite water heat exchanges with wet debris and bare ice without counting the inlet heat twice.
- A normal-diffusion slab replaces the uncontrolled single-temperature debris approximation when requested. The RC model remains an explicitly labelled reference.
- Separate wet/dry columns retain thermal profiles; only material actually reclassified by a moving wet boundary transfers enthalpy between classes.
- Upper/lower cover and a fixed-elevation covered footprint are explicit. Matched-footprint comparisons do not silently change the area being scored.
- Height-band geometry supports nonuniform retreat and geometric undercuts, with true ice volume separate from the basal hydraulic proxy.
- Lossless insertion of water/cover boundaries removes thin-patch averaging errors. Actual adaptive grid counts and caps are recorded.
- Local transfer closures report geometry, Reynolds/Froude numbers and domain refusals. A bypass needs its own hydraulic input for a nonfixed correlation.
- Generated meltwater is tracked per route. A globally small source does not automatically certify negligible feedback in a nearly disconnected branch.
- Zero-external-inflow melt requires an explicit downstream free-drainage contract and capacity; no infinite ratio is silently set to zero.
- CLI outputs are atomic. Package imports do not overwrite another model's global module cache.

## Execute

From this folder:

```sh
python verify.py
python model.py example_temperate_profile.json result.json
```

From a repository containing this folder as `ice_cliff_physics_v2`:

```python
from ice_cliff_physics_v2.model import forecast
result = forecast(configuration)
```

The earlier `example.json` still runs; its existing cumulative totals match the pinned previous implementation exactly. It is a synthetic regression fixture, not a glacier forecast. The new example is also synthetic and declares its net heat boundary, cover position, temperate phase and ledge assumptions.

The standard-library core requires Python3.9+. Plotting uses NumPy/Matplotlib. Field-temperature diagnostics, when included, additionally use SciPy and have their own frozen protocol. Do not confuse the field diagnostic with the complete coupled forecast.

## Verification and controlled studies

`verify.py` records isolated test output in `VERIFICATION.json`. `thermal/` contains an independent finite-volume slab, classical Robin eigenfunction reference, meaningful space/time comparisons and numerical stress tests. `face_geometry/` checks nonconvex volume against polygon area and independent analytic evolution. `reviews/` retains pre-fix failures and post-fix adversarial checks.

`study/SPEC.md` and `study/COUPLED_CASES.json` define the new synthetic comparison. The pre-execution Git commit and frozen source/input hashes accompany every case record. All completions and refusals remain visible. The mechanism packet also fixes 12 constitutive boundary cases; these verify known heat-transfer limits and quantify proxy errors, not a newly discovered natural phenomenon.

Core statements to keep separate:

1. Numerical conservation and analytic agreement establish implementation/approximation properties.
2. Synthetic scenarios establish conditional behavior under their supplied boundaries.
3. Temperature hindcasts validate only the supplied thermal observation problem.
4. Full cliff melt, survival, collapse and glacier-scale transfer require additional matched evidence.

## Essential limitations

The water model uses quasi-steady common-stage ratings and an imposed downstream friction-slope closure, not resolved free-surface CFD. In profile mode horizontal ledges are explicitly adiabatic; undercut diagnostics do not predict fracture or collapse. Wet-debris modes currently require temperate ice and no unresolved face cold debt; freezing, cold-ice conduction and pore-water phase change remain unsupported. Wet/dry class transfer uses donor-average normal profiles, not fully resolved material advection. Fixed/fractional cover is prescribed, not a validated grain-transport model. Atmospheric ports require independently prepared net fluxes; equal absorbed heating above and below water is a controlled boundary, not automatically the same weather. No default coefficient is claimed universally field-valid.

The previous Swiss field failures are not erased or reclassified by these improvements. Existing results and Claude's active work remain unchanged.
