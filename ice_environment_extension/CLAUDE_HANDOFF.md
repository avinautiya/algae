# Integration handoff — preserve active work

This is an additive package. Import into a new directory; do not replace the
existing production solver, frozen specifications, checkpoints or raw results.
It was prepared against public branch head `142c0f056f6bf2995559172c4286e8ae97073a98`
only as a publication base. No production outputs were regenerated here.

## Verified now

`python3 verify.py`: **47 tests pass** (2 data, 12 convex geometry, 6 energy,
9 transmission, 18 integration). Core checks use the standard library.
`VERIFICATION.json` contains the exact suite output. Independent geometry review
confirmed the stepped polygon and face-then-floor removals are disjoint.

Actual data: 391,911 corrected Himalayan profile points; 11 profiles, 83 epochs,
9 cliff IDs. Source and error audit in `data/`; full local archive includes raw
selected members, while a source-only repository distribution may omit them.
Never claim source full-archive MD5 verification from partial range downloads.

## Use the correct component for the physics

- `forcing/energy.py`: fixed-area bare/debris columns. Cold reserve is an energy
  debt, not a resolved cold-ice PDE. Area evolution needs its own conservative
  remapping/storage term. No phase refreezing or finite-domain exhaustion here.
- `transmission/transmission.py`: passive ratings and common stage. Aperture
  head is integrated over height, preserving the actual floor. Use explicit
  coefficients. Thin-opening rates do not include a long path's viscous loss.
- `transmission/filled_slot.py`: an alternative fully wetted, narrow smooth slot
  with viscous and entrance losses. It refuses partially filled ends. Check its
  Reynolds/development/roughness conditions before use; no turbulent validation.
- `integration/heat_delivery.py`: deplete a single water-heat packet serially
  between floor and face; never apply the incoming sensible energy twice.
- `integration/stepped_cliff.py`: exact finite vertical-step extrusion, with
  fixed crest and homogenized face retreat. This is not undercutting, fracture
  or arbitrary terrain. `geometry/cliff_geometry.py` is convex-patch-only.
- `integration/predict_cliff.py`: a ready CLI consuming net latent packets.
  It predicts face/floor melt separately and audits their combined polygon
  volume. The example is deliberately labelled synthetic.

## Next bounded integration milestone

1. Preserve the current production freeze and finish its existing gate. Do not
   reclassify any old plan-view run as a full cliff result.
2. Add a separate environmental run mode. Define a common extruded floor/face
   footprint, exact coordinate frame and the water route reaching each patch.
   Feed compatible hydraulic heat fluxes to cold/debris storage, then only net
   latent heat to geometry. Report source, storage, outgoing sensible heat,
   latent heat, phase-water mass and geometric ice removal separately.
3. First verify a single constant-area bare face, then a shielded face, then
   combined face/floor. Compare both splitting orders and halve the coupling
   step. Evolving area must not create debris/cold energy by multiplication.
4. Add nonzero-transmission alternatives using physically compatible geometry.
   Under-rock water below a debris crest is allowed; below intact ice is not.
   A lower lateral/englacial route is a new path, not a leak coefficient. Do not
   extrapolate the tiny-gap inviscid scan as a field prediction.
5. Test **4–6 cases bracketing a recovery boundary** at two resolved meshes,
   with at least one finer case near a class change. Keep dynamic and
   prescribed-contact models identical except their contact operator. Compare
   floor and face melt on the same support. Refuse a headline if closure spread
   or numerical uncertainty covers the effect. Freeze rules before evaluation.
6. Use Himalayan profiles to constrain observed geometry ranges and assess
   whether a vertical step is representative. Do not infer debris masks,
   discharge, water temperature or clearing/collapse times from unlabelled
   point clouds. Keep survey/cliff dependency in any split.

## Paper gate

The methods insert can document implemented equations and verification. The
main environmental claim may only change after coupled results exist. Do not
delete collapse, nonhydrostatic crest or site-input limitations because a helper
function or sensitivity range now exists. Collapse needs a separate mechanics
model or a clearly labelled prescribed perturbation, with detached solid ice
kept distinct from latent melt and runoff.

The candidate new result is a consequential error from **prescribed hydraulic
contact in a cliff-face melt forecast**, robust to transmission and environmental
forcing. Standard energy balance, polygon area and hydraulic ratings are not
the novelty. If the effect vanishes, record that null; do not tune site ranges
to preserve it. See `geometry/ENVIRONMENTAL_NOVELTY_GATE.md`.
