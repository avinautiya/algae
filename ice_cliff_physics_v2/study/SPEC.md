# Pre-execution coupled verification study

This synthetic study tests implementation and approximation errors. It does not
test field accuracy, historical originality, award readiness or superiority over
professional glacier models. Analytic limits and component tests have already
been examined. The following coupled outcomes have not been run when this spec
is written.

## Fixed cases and comparisons

`COUPLED_CASES.json` defines 30 cases. The first 22 use two areal heat capacities, two prescribed
net covered-surface heat fluxes, fixed material conductance and shared hydraulic
forcing. Each material/source pair includes the lumped reference, 16/32 normal
cells, outer-step refinement, and a labelled no-water/debris-exchange control.
Two additional runs refine the face-height mesh. There is no parameter search,
gain fit, model selection or case removal based on results.

Eight further bare-cliff geometry cases compare a planar face with aligned
height profiles at two supplied water temperatures, with separate outer-step
and initial-height-grid settings. These test whether a bulk melt total conceals
different geometric recession. They do not test mechanical collapse. No
threshold for an impressive-looking undercut is selected after execution.

All comparators start with zero stored thermal energy and identical temperate
geometry. A fixed elevation interval defines the covered footprint in every
case. Its projected area must remain identical; total face melt is not a
substitute for covered-footprint melt. New exposed bare area has explicit zero
incoming cold debt. The absence of water/debris exchange in the control is an
intentional boundary intervention, not a claimed professional model.

The net heat source is imposed directly on all covered material, including wet
material. This is a controlled heating boundary, NOT an assertion that the same
weather produces equal absorbed heat above and below turbid water. Atmospheric
and frictional physics beyond the specified ports are not being validated.

## Metrics and acceptance

- Record every completion or refusal, runtime, input hash and source fingerprint.
- Require water/latent/geometry bookkeeping closure; report maximum time-series
  energy residual normalized by cumulative supplied energy, not only its final
  signed sum. Numerical tolerance is 1e-8 relative and 1e-9 m3 for geometry.
- Report covered and bare face melt separately, fixed covered area, outgoing
  water heat, thermal-column storage and geometrical undercut diagnostic.
- Compare RC/FV only on common material, source, forcing and footprint. Report
  normal-mesh and outer-time differences separately. No extrapolated field
  confidence interval is inferred from numerical refinement.
- Adaptive insertion at current exposure boundaries is lossless. Record actual
  height-band counts, not just initial grid size. Time-step refinement changes
  boundary-insertion history as well as time discretisation; do not claim a pure
  temporal convergence order from that combined test. Independent analytic
  geometry and conduction checks provide the separated component orders.
- Compare coupled versus disabled-exchange controls without renaming flow
  connection as thermal recovery. Classify warming/cooling from the results;
  the standalone fixed-coefficient steady formula is a separate limit, not an
  exact assertion for evolving geometry.
- An undercut is a geometric output; no fracture, collapse time or cliff
  lifetime is predicted. New horizontal ledges have explicitly zero heat flux.

## Independent constitutive comparison

The separate mechanism packet fixes 12 wet-boundary networks before the resolved
comparison. Compare the normal-diffusion solver to its independent steady
resistance formula, the RC reference, wet-treated-as-bare proxy and dry-exchange
control on the same wet footprint. Retain all cases and all signed water heat.
Known insulation and advective cooling are not relabelled as new phenomena.

## Versioning

Save code/spec/input hashes before execution and publish the freeze commit. A
numerical defect found later requires a documented new version and retained
previous output. Never repair a run by silently relaxing applicability gates.
