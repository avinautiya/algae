# Environmental cliff-extension package

Additive contribution, 1 November 2026. **47 component checks pass** with the
system Python. This package leaves the existing solver, active jobs and frozen
results unchanged. It is not a completed or field-validated Himalayan glacier
model, and these components are not evidence of a new mathematical theorem.

## Delivered

- **Real Himalayan geometry:** selected published Lirung/Langtang profile files,
  provenance and error audit; 391,911 points, 83 profile epochs, 11 profiles,
  9 cliff IDs. Every profile is plotted. See `data/HIMALAYAN_INPUT_AUDIT.md`.
- **Environmental energy:** signed atmospheric balance helper, cold-content
  reserve, direct bare-ice energy and finite debris heat storage, with separate
  conservation ledgers. Fixed-area lumped columns, not a cold-ice PDE.
- **Obstruction transmission:** explicit parallel/under-rock ratings that cannot
  flow below intact ice, common-stage allocation, exact vertical aperture-head
  integration, an alternative viscous fully wetted slot and diagnostic outputs.
  Actual pathway geometry and coefficients remain required inputs.
- **Hydrostatic diagnostics:** aspect ratio, Froude and curvature acceleration.
  Unknown curvature is preserved. None is labelled a validation pass.
- **Moving geometry and volumes:** energy-constrained convex patch retreat and
  a nonconvex vertical-step extrusion with independent face/floor melt. Corner
  removal is counted once. Uniform retreat cannot produce undercuts or collapse.
- **Executable predictor:** `integration/predict_cliff.py` maps explicitly
  supplied net latent-energy packets to normal retreat, floor lowering, exposed
  areas and separate melt volumes. It refuses invalid/exhausted geometry. Its
  synthetic example is not an environmental forecast.
- **Conditional recovery bound:** a proved sufficient no-contact inequality for
  a stated uniform-toe model, with explicit invalidation conditions.
- **Paper material:** `PAPER_METHODS_INSERT.md`, `ENVIRONMENTAL_SCOPE.md` and the
  focused prediction/novelty gate in `geometry/ENVIRONMENTAL_NOVELTY_GATE.md`.

## Run now

From this directory (or substitute its repository prefix):

```sh
python3 verify.py
python3 integration/predict_cliff.py integration/synthetic_energy_packets.json integration/synthetic_forecast_result.json
python3 forcing/synthetic_demo.py --output forcing/synthetic_ledger.json
python3 transmission/synthetic_scan.py
```

Core tests and examples need only Python's standard library. Matplotlib is
optional for PNG/PDF plots. `data/summarize_profiles.py --plot` requires it.
Data acquisition requires network access but tests do not.

Retrieve the published profile data if not included in your distribution:

```sh
python3 data/audit_zenodo_profiles.py --output data --max-headers 100
python3 data/fetch_profiles.py --root data
python3 data/summarize_profiles.py --root data --plot
```

No full 60 GB dataset download occurs. Bounded tar-member downloads total about
27 MB; the full archive checksum is not claimed verified. Source attribution,
license and local member hashes are included.

## Integration that remains necessary

The next solver milestone must connect compatible hydraulic/thermal face
fluxes to the moving geometry, refine the recovery boundary and compare face
melt on identical support. It must account for changing-area storage and avoid
assigning one water-heat packet to both floor and face. A diagnostic or a closed
ledger alone cannot establish environmental predictive accuracy.

Atmospheric cliff balance and water effects are established prior art. The
specific proposed contribution is **whether predicting hydraulic recovery
improves a consequential cliff-face melt prediction over prescribed contact**.
The new scan deliberately probes whether nonzero transmission erases that
mechanism. Neither a favorable scan nor added complexity establishes priority.

Collapse mechanics, differential face shape, glacier motion and complete
nonhydrostatic free-surface coupling are not supplied here. Do not remove those
limitations from the paper until an implementation and independent evidence
actually address them. See `CLAUDE_HANDOFF.md` for the bounded integration order.
