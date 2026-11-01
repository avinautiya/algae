# Additive integration handoff

Import this entire directory separately and run `python3 verify.py` (16 tests).
Do not replace the current v4 CFD solver, raw outputs or frozen specifications.
The vendored sources are inputs to this reduced integration model, not patches
to production. The example and development outputs are synthetic.

## Completed here

- Coupled executable CLI with actual geometry/latent-volume checks, one water
  heat allocation, signed atmosphere, extensive debris/cold states, explicit
  new-area enthalpy, and bulk warming of fresh cold ice.
- Distinct upstream stage and downstream contact depth, with a necessary head
  check and explicit refusal of submerged debris or excessive phase-water
  hydraulic feedback.
- Six attempted development variants, three coupling timesteps, reproducible
  figures, and preserved refusal reasons. Not a frozen independent campaign.
- Pinned, verified actual Lirung station meteorology and quality audit.
- Literature-backed comparison protocol distinguishing public CliffEBM from
  the complete 2016 evolving-geometry model.

## Next bounded work

1. Read README.md, AUDIT.md, REVIEW.md and the development log before using the
   outputs. Keep independent hydraulic-rating/backwater and lumped-face caveats.
2. Use current v4 reference outputs to benchmark the routing/contact component
   first on geometry it actually represents. Do not claim the plan-view floor
   solver independently validates a vertical cliff-face melt forecast.
3. Add a resolved or independently verified downstream contact-depth operator
   before environmental use; arbitrary Manning slope is an assumed input, not
   data inferred from a Himalayan image.
4. Prepare net face fluxes from the real AWS data using documented terrain,
   radiation and atmospheric assumptions or the professional EBM. Keep station
   QC choices versioned. Source actual or explicitly scenario-bound water Q/T.
5. Freeze matched-physics dynamic-contact vs prescribed-contact predictions and
   score common-support cumulative FACE melt. Independent reference or field
   data are needed to call one better. A difference between two reduced models
   alone is not forecast improvement.
6. Run 4–6 boundary-bracketing cases with explicit transmission uncertainty only
   after the compatible coupling gate passes. Preserve near-boundary unresolved
   classes and count shared histories once. Do not optimize cases for a large
   or favorable effect.

No collapse law, neural network or new theorem has been claimed. Those would
require separate derivation and evidence. The scientific target is a robust,
quantified environmental prediction error from prescribed contact—not added
complexity by itself.
