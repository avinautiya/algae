# Conservative planar cliff patch geometry

## Implemented and checked

`cliff_geometry.py` implements a **convex ice cross-section patch**, not the full
cliff–stream domain. The parent floor solver remains separate. The concave
L-shaped solid domain at a cliff toe is not represented by this primitive.

```python
from cliff_geometry import retreat_by_energy, segment_exposure
result = retreat_by_energy(vertices, inward_normal=(1, 0),
                           energy_J=net_latent_J, width_m=width)
```

Coordinates are horizontal distance and absolute elevation in metres; width is
perpendicular to the section. For a vertical ice wall with ice on its right,
inward normal `(1,0)` gives horizontal retreat. For ice below a horizontal
floor, `(0,-1)` gives floor lowering. These are **different moving surfaces**.
Do not feed floor lowering into the face or count its melt twice.

For convex polygon P, inward unit normal n and initial supporting plane
c=min(n·x), define P(d)=P intersect {n·x >= c+d}. Solve the monotone area equation

    rho_i L W [area(P)-area(P(d))] = Q_latent.

The shoelace area uses translated coordinates to reduce cancellation. Bisection
enforces the supplied latent budget; this is an energy-constrained geometric
operator, not a derived heat flux, flow law, or fracture law. If n does not
coincide with an existing edge it creates a corner cut; callers must explicitly
check that this is their intended physical interface. Positive energy exhausting
the patch is refused, as are nonconvex/self-intersecting inputs and negative
energy. Refreezing, cold-ice evolution and collapse are not implemented here.

Before calling this operator the thermal module must remove cold-content
warming and other storage from heat input. A positive total heat input below
the cold reserve means **zero** net latent energy and therefore zero retreat;
passing that total input directly would be a physical error despite a closed
geometric ledger. Debris sensible heat must similarly be charged once upstream.

`segment_exposure(a,b,water_level_m,debris_top_m)` splits a prescribed straight
face into three disjoint lengths: covered, bare submerged, bare exposed. The
mantle is an elevation mask; wet debris is covered, not counted as submerged
bare ice. Multiply each length by patch width for areas. Do not infer masks
from unlabelled point clouds or assume this predicts exposure/transport.

## Coupling requirements before environmental prediction

1. Register disjoint floor and face patches and their interface areas; use one
   water-side heat removal for each corresponding ice-side latent gain.
2. Supply separately modelled atmosphere, water, debris and conductive fluxes,
   net of storage, in joules integrated over their *current* physical areas.
3. Recompute water contact and debris coverage after moving geometry. Remap
   water volume, momentum and enthalpy conservatively; a changed face alone
   cannot be fed back as an arbitrary bed-height edit.
4. Check the combined ice area/volume partition independently. Two overlapping
   convex patches cannot be summed. A whole concave section requires a proper
   nonoverlapping mesh/level-set/remeshing implementation, not a convex hull.
5. Report mass as rho_i removed_volume. The displaced water is not generally
   equal to that ice volume because densities differ.
6. Treat prescribed mechanical removal as a separate measured/bounded event
   with its own solid mass ledger. No slope threshold or unsupported fracture
   criterion is supplied. This primitive cannot predict cliff lifetime.

Run verification: `python3 -m unittest discover -s research/environment_extension_2026_11_01/geometry -v`.
Twelve independent cases cover analytic rectangle/triangle geometry, normal
direction, width scaling, orientation, map coordinates, exposure partition,
zero energy, exhaustion and invalid geometry. They verify code, not field melt.

## Prior art and novelty boundary

This is standard polygon clipping, area accounting and the Stefan latent-energy
condition; there is **no new mathematical theorem claim** for these primitives.

- [Buri et al. 2016](https://doi.org/10.1002/2016JF004039), a physically based
  3-D cliff evolution model, precedes moving cliff geometry; the accessible
  [author archive](https://eprints.whiterose.ac.uk/id/eprint/124673/8/Buri_et_al-2016-Journal_of_Geophysical_Research__Earth_Surface.pdf)
  is the primary full-text source to consult for an equation-level comparison.
- [Evatt et al. 2024](https://doi.org/10.1017/aog.2024.34) reassess cliff dynamics
  and the surrounding debris-lowering reference frame. Moving-frame corrections
  and cliff-slope evolution should not be presented as new here.
- [Petersen et al. 2024](https://esurf.copernicus.org/articles/12/727/2024/)
  provide observed stream controls on cliff evolution/survival. Stream support,
  abandonment and reburial are not first discoveries of this project.

The worthwhile *testable* extension is whether incision-driven hydraulic
disconnection/reconnection adds out-of-sample predictive value to a properly
energy-balanced cliff model, beyond prescribed-contact and thermal-only
baselines. That is not demonstrated by this geometry package. It needs the
coupled contact calculation, common melt support, defensible transmission and
forcing ranges, convergence at the decision boundary and independent profile
evaluation. Observed Himalayan profiles can constrain geometry/change but do
not supply unmeasured stream temperature, contact or debris labels.
