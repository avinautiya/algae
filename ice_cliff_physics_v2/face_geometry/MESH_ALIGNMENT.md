# Lossless sharp-class mesh alignment

```python
from mesh import split_at_elevations

aligned = split_at_elevations(profile, [water_stage, cover_lo, cover_hi],
                              max_bands=4096)
```

The helper inserts all supplied, sorted unique finite elevations that lie strictly
inside an existing height band. Every child copies the parent's front coordinate
exactly. Original shape, material mask, horizontal ledges, floor and heat budgets
are unchanged. The continuous polygon is identical; computed volume is identical
up to floating-point summation roundoff. This is remeshing of the same geometric
state, not ice melting, area creation or a thermal-state remap.

The maximum permitted count is checked before constructing any replacement
profile. Exceeding it raises with the required and allowed counts. No elevation
is dropped, no coarsening is performed, and the input frozen Profile is unchanged.
Outside-face elevations and already present interfaces are no-ops; an entirely
unchanged mesh returns the original object. Distinct representable elevations
are not epsilon-merged, so very thin interfaces remain auditable.

Call alignment on the **current old-step profile before heat allocation**. Both
the water stage and cover boundaries must be the same supplied coordinates used
for thermal partitioning. Then each band's exposure class is sharp, rather than
a band-average mixing a thin basal bare-water strip with the covered/air region.
The integrated adapter's existing API is unchanged.

Alignment eliminates that classification alias for classwise-uniform latent
fluxes. It does not resolve heat variation inside a class, a continuously sloping
initial face, turbulent flow, ledge heat, entrapped air or material transport.
Moving stages preserve earlier front differences; historical interfaces are not
coarsened away. Therefore count can grow with the number of forcing steps. Record
the declared cap and actual maximum count, bound campaign cost accordingly, and
refine timestep as well as the initial geometric representation. A cap refusal
must remain a recorded domain/runtime failure, not disappear from scoring.

Eight new tests verify nonplanar polygon/volume/ledge identity, exact initial
step preservation, no-op behavior, preflight-cap refusal, invalid-input rejection,
retention of distinct nearby floats, deterministic many-split geometry, and
identical classwise-flux profile advancement across aligned starting grids.
All thirty profile tests pass (fifteen geometry, seven integration, eight mesh).
