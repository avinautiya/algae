# Height-resolved conservative face geometry

The previous `SteppedCliff` represents one planar vertical face. This optional,
independent geometry operator represents a piecewise-constant horizontal front
`x_front(z)` and can therefore represent geometric undercutting without inventing
a fracture rule. It does not replace the parent coupled model automatically.

The initial body has a horizontal floor from the upstream edge to the lowest
face position, a fixed back edge/bottom/crest, and a constant extrusion width.
Below the floor, ice occupies the full horizontal domain. Above the floor, each
height band contains ice between its front and the fixed back edge.

## API

```python
from face_profile import Band, Profile, advance, face_partition

g=Profile.uniform(upstream,back,bottom,face,floor,crest,width,n_bands=32)
g_next,ledger=advance(g,band_net_latent_J,floor_net_latent_J)
areas=face_partition(g_next,supplied_water_stage,
                     debris_intervals=[(covered_bottom,covered_top)])
```

`advance` needs one net latent energy for every existing band. It removes the
corresponding face rectangles first, then lowers the floor using the **new** basal
footprint. The removals are disjoint. The lowest band extends downward to include
the newly exposed face; that new area receives no uncharged face melting within
the same split step. Density, latent heat and all geometric/data inputs are
validated. Exhaustion of the finite ice body and changes below coordinate
floating-point resolution are refused rather than clipped.

The ledger separately reports energy-implied melt, a stable geometric removal,
remaining-volume subtraction, and subtraction roundoff. The polygon can be
checked by an independent shoelace-area computation. The planar uniform limit
reproduces the old step geometry for any initial band count.

`face_partition` gives exclusive bare-water, bare-air, debris-water and debris-air
areas within each height band, using exact intersections with supplied height
intervals. These are **Eulerian height masks**, not advected material parcels.
Water stage is supplied; the operator does not solve wetting, entrapped air,
hydraulic access, backwater or a nonhydrostatic interface.

## Important newly visible surfaces

Nonuniform fronts create horizontal ledges. Downward-facing ceilings are
geometric undercuts; upward-facing ledges are terraces. Their areas and
orientations are explicitly reported. They are **not** silently assigned the
vertical-face heat flux. A future full thermal integration must supply their
boundary conditions and storage/material remapping. The example prescribes zero
ledge latent flux explicitly. This is not an assertion that natural ledges are
thermally inactive.

`maximum_horizontal_undercut_m` is geometry only. No tensile strength, bending,
fracture, stability threshold, collapse time, cliff lifetime or debris failure is
inferred from it. The finite back-connected ice body remains connected by
construction.

## Verification and reproducible example

Run:

```sh
python -m unittest discover -s . -p 'test_*.py'
python run_example.py --plot
```

Numerical tests use only the standard library. Optional PNG/PDF plotting uses
Matplotlib. Fifteen tests check independent polygon area, disjoint face/floor
removal, the planar limit, signed-datum translation, exclusive height masks,
finite-domain refusals, prescribed nonuniform-flux height refinement, and
independent continuum solutions for both the height-resolved and reduced uniform
models.

The demonstration is a **synthetic prescribed-net-latent-flux geometry test**:
2000 W/m² below a fixed stage, 60 W/m² above, 300 W/m² on the horizontal floor,
zero on generated ledges. It is not a solved water heat budget, atmospheric
surface-energy model, CFD experiment or field prediction. Its numerical front
approaches a closed-form piecewise moving profile; the reduced uniform model
approaches its own distinct closed-form solution. A volume difference between
them is therefore a geometry-representation effect under this forcing, not a
discovered field improvement or a claim that the planar algorithm is internally
inconsistent. The exact derivation is in `ANALYTIC_BENCHMARK.md`.

No historical-originality claim is made for height bands, prescribed-flux front
kinematics, the area formula or the elementary analytic solution. The contribution
here is an executable conservative extension and a falsifiable verification
problem that exposes a limitation the planar geometry cannot represent.
