# Matched-footprint melt repair

Pinned source: `f3602e35e1595bcdbee4208a354a34944f1ae974`.
This analysis is independent. No production source, frozen prediction, raw
result or running job was changed, and no CFD was launched.

## Main finding: existing exports can repair this without rerunning CFD

The mismatch is specific and recoverable:

- `cases.py` reports `cum_melt_toe_reach` over the full2m x0.5m toe reach.
- Imported `topology.run(..., spanning_strip)` evolves and counts ONLY a
  central1m x0.5m control. It explicitly excludes the remaining open reach.
- The corridor topology model counts both longitudinal branches over the whole
  2m reach. Its existing scalar already has the same reach footprint as CFD.
- `reduced.fractional` and `fixed` use a2m x0.5m uniform toe control; their
  scalar melt counts the whole reach, including diagnostic melt when geometry
  is fixed. Do not use final bed changes as their melt operator when fixed.

Therefore the reported roughly−99% strip connectivity melt error is not a fair
model score. However, not every existing melt score is mismatched: corridor,
fractional and fixed full-reach volumes already share the declared footprint.

The CFD runner saves `toe_floor_x` at EVERY morphology record: mean ice-bed
elevation across the full toe width, for each streamwise reach cell. These
profiles suffice to reconstruct actual geometric melt over any full-width
axial subregion. Full2-D fields are not necessary for this particular repair.

## A single physical ROI shared by all methods

Use the central1m of the toe reach, across its full0.5m width:

    absolute x3.5..4.5m, toe y0..0.5m, area0.5m².

This is the whole strip-control footprint. The reduced observation operators:

| Reduced model | Scalar support | Common central ROI volume |
|---|---|---|
| connectivity strip | central1m x0.5m | unchanged scalar |
| connectivity corridor | both branches over2m x0.5m | 0.5 × scalar |
| fractional/fixed/bare | uniform2m x0.5m | 0.5 × scalar |
| porous | no pinned reduced prediction | abstain |

The factors are exact restrictions of the reduced models' own axially uniform
control volumes. They do not claim that CFD melt is axially uniform. CFD is
integrated from the actual saved profile, with exact cell-overlap weights.
An arbitrary transverse subset of the corridor cannot be reconstructed from
its aggregated scalar without per-branch melt; this adapter deliberately uses
the FULL width. An arbitrary open portion of a strip model is unsupported and
must not be filled with invented melt.

## Precision and physical accounting

The CFD profiles were rounded to5 decimal places. Independent initial/final
profile rounding implies a conservative deterministic volume bound of

    area ×10^-5m =5×10^-6m³

for the0.5m² ROI. This is storage precision only. It is not a confidence interval
or the spatial/temporal/closure uncertainty. If initial elevations are known
exactly one could tighten this bound, but the supplied tool does not assume it.

The observer measures actual `zb_initial-zb_final`, preserving the distinction
between geometry change and nominal latent-energy melt. The pinned
`Sim.lower_bed` clips bed elevations to `pool_floor` after computing requested
lowering. It can therefore count requested latent volume that differs from
actual lowering. The tool reports full-reach profile-vs-cumulative disagreement
outside rounding precision; investigate it instead of silently equating them.
Reach cells may remain above that limit even when entrance ramps do not.

## What is ready

`observation.py` provides:

- exact fractional overlap weights on nonuniform rectangular finite-volume grids;
- masked actual geometric melt from before/after fields;
- latent-energy/volume conversions;
- existing mean-profile reconstruction with rounding bound;
- reduced-control support mapping with refusal of unmodelled regions;
- pure future `export_sim_interval` hook comparing actual lowering against heat.

`repair_eval.py` reads the ORIGINAL frozen predictions and existing medium CFD
JSON exports and writes new scope-matched scores over the common ROI. It does
not fit models, rerun predictions, modify the freeze or overwrite raw data.
Missing records produce `unavailable`, never zero error. It also decomposes
actual CFD whole-reach melt into central-control plus remaining open reach.

Run identity checks require the raw `dx`, `layout`, `removal_day`, `days` and
physical `dtm` to match the frozen case. Missing keys or disagreements refuse
analysis. Truncated endpoints, a missing time-zero profile, or an initial floor
inconsistent with the pinned geometry do not receive scores. Output must be
outside the entire raw evaluation directory and cannot replace the spec.

These checks do NOT authenticate the execution merely from the filename or
matching metadata. The raw schema omits full Q/T forcing, complete geometry,
solver commit and configuration hashes. A production hash manifest and original
execution provenance are still needed to link records to the frozen run. Input
SHA256 values identify supplied bytes; they do not certify the run that made them.

## Availability limit

The public pinned Git tree contains source and the frozen reduced predictions,
but not the referenced M3 `results/eval/E*_medium.json` files. No corrected real
M3 volume or ranking can therefore be calculated from the supplied snapshot.
The required raw records are already held by Claude: copy/export those JSONs,
then execute the command below. No fresh CFD campaign is needed for this repair.

## Execute

From the repair packet directory:

    python3 -m unittest discover -s melt -t . -p 'test_*.py' -v
    python3 -m melt.repair_eval /path/to/existing/results/eval snapshot/frozen/production_spec.json --output matched_melt.json

Fifteen tests cover independent fractional overlap, control/open decomposition,
energy conversions, model support, raw-data absence, synthetic scope matching
and floor-clipping disagreement, metadata mismatches, truncated/incorrect initial
profiles and immutable raw-input paths (including hardlink aliases). Test examples
are explicitly synthetic.

For future2-D exports call the pure hook immediately after `lower_bed` and BEFORE
another hydraulics solve, using `before_zb=sim.zb.copy()` saved before lowering.
The supplied interval is physical morphology time, not hydraulic spin-up time.

Any improved common-ROI score remains a conditional comparison with shared
thermal closure and idealized reduced geometry. It does not demonstrate total
cliff or glacier melt accuracy or validate debris heat transfer.
