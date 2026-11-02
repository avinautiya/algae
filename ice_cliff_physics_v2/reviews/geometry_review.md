# Independent geometry integration review

Reviewed the parent package/model.py and exercised the integrated forecast without
editing parent files. This is a synthetic verification audit, not a field test.

## Accounting and implementation checks: pass

A deliberately nonplanar initial face with fronts [0,-0.25,0.2,0.1]m over four
height bands was combined with a partially submerged fixed cover [-0.06,0.34]m,
resolved wet/dry debris, positive atmospheric heat and floor lowering. The
physical remaining profile differed from its planar routing proxy by about
-0.08765m³, so the check was not a disguised planar case.

- Global energy residual:2.07e-8J over600s.
- Physical geometry minus phase volume:1.03e-15m³.
- Moving lower-fraction cover variant also passed: energy residual-2.45e-8J,
  geometry residual-1.66e-15m³.
- Four-class latent allocation is captured before geometry advancement and its
  sum matches the face phase total. Positive energy on absent classes is refused.
- The parent's final conservation check uses actual profile volume, not the
  planar proxy. The proxy is used only for reduced basal routing geometry.
- Fixed-height cover retains its area while new lower bare face is exposed;
  column storage/import accounting remains consistent in the exercised case.

No blocking volume/energy integration defect was found in these tests.

## Finding G1: timestamp ambiguity in outputs (medium; simple label repair)

The row time_s is end-of-step. Geometry, profile and port_states are end-of-step,
whereas contact_stage, wet areas and covered_z_interval are beginning-of-step
values. hydraulic_sample_time_s already marks the hydraulic staggering, but the
surface-mask fields do not explicitly state their own time.

For the moving-cover example, the final row reports the old interval
[-0.0801760,0.3198468]m while recomputing the interval from end geometry gives
[-0.0801956,0.3198298]m. This is correct explicit stepping, not an energy defect.
However scoring wetting or mask transitions using time_s would add one step of
lag. Recommended repair: add surface_partition_sample_time_s and
state_sample_time_s, and use hydraulic_sample_time_s for routing event times.
Optionally output the next-step interval separately. Do not overwrite old values
with new geometry while retaining the old energy calculation.

## Finding G2: height resolution is independent of thermal depth resolution

A second synthetic audit ran exactly the same initial nonplanar shape, fixed
cover and6h forcing on successively split height bands. Results:

| Bands | Basal retreat m | Final floor z m |
|---|---:|---:|
|4|0.004906757|-0.087040761|
|16|0.003949767|-0.087040811|
|64|0.004240813|-0.087040794|
|128|0.007029445|-0.087040650|
|256|0.007029445|-0.087040650|

The initial uncovered basal patch is only2cm high. Coarse bands average its
water heat with covered/air portions, so basal retreat is not converged at64
bands even though energy ledgers pass.128/256 agree exactly for this case. This
is a documented band-average discretization effect, not a conservation defect.
A campaign claiming basal/capture accuracy must refine face height separately
from normal debris cells and temporal substeps. An initially mask-aligned grid
could reduce the required count, but would be a declared discretization policy;
no parent feature was added in this review.

The floor difference is much smaller than the basal-front difference, so a
floor-only convergence check would miss this issue. Near a competition boundary,
affected route length/contact quantities must also be checked.

## Interpretation limits that remain

- fixed_cover_interval means a fixed Eulerian height mask, not tracked grains,
  material transport or an unchanged Cartesian surface while ice retreats.
- Horizontal ledges/ceilings are explicitly adiabatic. Their area is tracked but
  they have no modeled exchange, material storage, air trapping or collapse.
- A basal routing proxy cannot resolve flow within an overhang. A common supplied
  hydrostatic stage is an engineering wetting/accessibility assumption.
- Profile mode prohibits spatial face cold/refreezing but the floor/bypass ports
  can still accumulate their own lumped cold debt under negative imposed heat.
  Descriptions of 'no cold reservoirs' should be read specifically for the face
  allocation contract, not as proof of a solved cold-ice field.
- Source-derived pressure depth is not necessarily full water depth for a sensor
  suspended above the bed; this review uses explicitly synthetic geometry.

Files geometry_counterexample_results.json and geometry_height_refinement.json
retain the actual numerical evidence. No scenario was tuned to produce a
favorable field claim, and no frozen parent campaign was altered.

Reproduce the height audit with `python reproduce_geometry_audit.py --package PATH_TO_PACKAGE --out NEW_RESULT.json`. The frozen complete synthetic input is `geometry_counterexample_config.json`; all five reported refinement rows reproduced exactly before delivery.
