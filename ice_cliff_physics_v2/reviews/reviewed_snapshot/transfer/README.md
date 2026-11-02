# Explicit local transfer closures

This additive module replaces a shared fixed coefficient only when a caller
explicitly selects an alternative. It does not select a winner, infer flow from
temperature, fit a field coefficient, or claim validated open-channel physics.

```python
from local_transfer import evaluate

diagnostic = evaluate(q, width, depth, config,
                      legacy_fixed_h_W_m2_K=old_physics_h)
h = diagnostic['h_W_m2_K']
```

An omitted configuration, `{}`, or `{'method':'fixed'}` preserves the supplied
legacy fixed coefficient for positive flow. At zero flow the returned coefficient
is zero because this interface describes **advective exchange**, not natural
convection or a stored-water reservoir. Geometry may have zero depth only at zero
flow. A positive flow over a dry cross-section is refused.

Configuration schema version1 supports:

- `method`: `fixed`, `db04`, `db03`, `gnielinski_smooth`.
- `fixed_h_W_m2_K`: explicitly prescribed coefficient, default legacy value.
- `nu_m2_s`, `water_k_W_m_K`, `prandtl`: supplied constant bulk properties.
- `gain`: default1; any other value requires nonempty `calibration_provenance`.
- `domain_policy`: `strict` by default. `allow_high_re_extrapolation` requires
  `extrapolation_reason` and stamps warnings; it never enables a low-Re fallback.

Defaults are nu1.8e-6m²/s, k0.56W/(mK), Pr13.5. These values are declared constant
approximations, not an inferred thermophysical-property model. Unknown fields,
nonfinite values and unsupported method names fail.

The flow diagnostics use rectangular area `w*d`, wetted perimeter `w+2*d`,
hydraulic diameter `4*A/P`, velocity `Q/A`, Reynolds `u*Dh/nu` and Froude
`u/sqrt(g*d)`. A numerical gate is not a proof of hydrostatic validity or a
validation of applying a smooth-tube law to rough ice and a free surface.

## Numerical-domain policy and sources

The DB closure uses `Nu=.023*Re^.8*Pr^n`; `n=.4` is the explicit published-style
comparator, whereas `n=.3` is the cooling-fluid comparator. A water stream warmer
than a melting ice wall is a cooling-fluid problem. Neither choice is silently
substituted for the other. The conservative default gate is Re10,000–120,000 and
Pr0.7–160, following the range warnings in the primary
[MOOSE implementation](https://mooseframework.inl.gov/docs/doxygen/modules/DittusBoelterFunctorMaterial_8C_source.html).
This is narrower than some engineering implementations and the earlier Swiss
comparator, intentionally requiring an explicit extrapolation decision for larger
glacial Reynolds numbers. It is not a claim that all uses above120,000 are false.

Smooth Gnielinski uses its own smooth-tube Darcy expression and no externally
inferred friction coefficient. Its gate is Re3,000–5e6 and Pr0.5–2,000. The lower
limit is a conservative choice above the2300 threshold in the primary
[MOOSE Gnielinski implementation](https://mooseframework.inl.gov/docs/doxygen/modules/classADWallHeatTransferCoefficientGnielinskiMaterial.html).
No source code was copied from MOOSE. These equations are independently
implemented mathematical comparators; no fitted Swiss multiplier is included.

`finite_heat_chain` is an exact serial, zero-C-surface heat-budget helper with
stable `expm1`/exponential evaluation. Every downstream section receives only the
upstream remainder. A surface at nonzero temperature or with debris storage needs
the separate wet-debris kernel. At zero discharge there is no advected heat;
warm standing-water energy is outside this helper.

Run `python -m unittest discover -s . -p 'test_*.py'`. Seventeen tests cover
geometry/mass identities, independent known-value DB evaluation, explicit
calibration, regime refusals, zero flow, serial subdivision equivalence and
finite-energy bounds at extreme thermal exchange strengths. Passing these is
implementation evidence, not field validation or a novelty claim.

## Integration rules

Evaluate toe and bypass separately from their own supplied width, depth and flow.
Do not use toe hydraulic diameter for the bypass. A failure must retain which
route/method/regime failed; do not retry with fixed h and report that as the
requested model. Preserve the configuration and its provenance in restarts and
results. For an ensemble, declare its methods before scoring and report all
refusals, including shared-case and full attempted-case counts.
