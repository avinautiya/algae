# Spatial heat-budget exclusion of route opening

This is a classical constrained-path argument used as an independent forecast
bound. It is not a claim of a new mathematical theorem or a new field discovery.
It extends the previous fixed-width two-route energy bound to a resolved ice
surface. Its practical role is to reject impossible recovery predictions before
running expensive coupled hydraulics.

## Precise discrete proposition

Let a rectangular finite-volume grid represent removable ice after perfect
instantaneous debris clearing. Cell `v` has elevation `z_v` and area `a=dx*dy`.
Assume the water surface can never exceed a specified upper bound `H` during
the forecast interval. Define the minimum excavation depth

    d_v = max(z_v - H, 0).

The inlet comprises the first row and the target comprises the last. Only
face-connected flow is admitted. Immutable walls are excluded from the graph.
Suppose a route of open cells forms by vertical melting alone. If `S` is its
excavated cell set, it must contain a simple inlet-to-target path `P`. Thus

    melted volume >= sum(v in S) a*d_v
                  >= min(P) sum(v in P) a*d_v = V_path.

Nonnegative vertex costs make Dijkstra's algorithm applicable. Conversely,
melting precisely the minimising path to `H` attains this minimum in the
discrete geometric optimisation. Equality only gives zero-depth connectivity,
not positive flow. Hydraulic head losses, finite water depth and realistic
distribution of heat can only require more than this optimistic construction.

If a declared corridor must occupy at least transverse width `w` in EVERY row,
minimise excavation independently in each row while dropping connectivity.
Sort the required depths from smallest to largest, fill the cheapest `w` metres
of transverse width (allowing a fractional last cell), and multiply by `dx`.
Summing gives `V_width`. Every admissible corridor satisfies this relaxed
constraint, hence

    V_melt >= V_width.

If full-cell opening is ALSO required, the stronger discrete result is

    V_melt >= max(V_path, V_width).

The code defaults to `V_width` alone. `resolved_cell_route=True` explicitly
opts into the discrete full-cell condition. It never silently applies a
whole-cell energy requirement to a narrower continuous corridor. In either
mode `E_latent >= rho_i*L*V_bound`.

If an independently justified upper bound on all usable heat is smaller than
this energy lower bound, no admitted corridor can open within the interval.
Equality remains inconclusive. The condition is necessary, not sufficient.

## Assumptions that determine whether the certificate is useful

- `H` must bound the water level EVERYWHERE, including transient backwater and
  hydraulic jumps, not merely equal the mean inlet stage. A kinetic-energy or
  storage-driven level excursion can invalidate an unjustified stage bound.
- Include initial sensible heat and every admitted heat source in the upper
  heat budget. Cold ice increases required energy, so ignoring its warming is
  optimistic. Mechanical removal, collapse or fracture invalidates a
  melting-only exclusion unless separately bounded.
- Elevations must describe ice, not a removable grain layer. Taking perfect
  debris clearance as free is deliberately optimistic; it isolates substrate
  geometry from remaining sediment. Grain transport energy is not latent heat.
- Walls must really be immutable. Do not mark porous deposits as walls merely
  to obtain exclusion. Existing wetted cells require no excavation and can be
  included as zero-cost vertices.
- This bound does not include seepage or subgrid channels. For a continuum
  problem with no minimum corridor width, a vanishingly thin path costs
  arbitrarily little. `V_path` is a whole-cell resolved-grid statement and generally
  weakens with transverse refinement. A physical minimum width must be
  justified before claiming a finite continuum energy requirement. In
  particular, a 0.1 m corridor through a 1 m-wide cell need not melt that whole
  cell. The default continuous-width bound permits that narrower opening.
- The width condition is a projected transverse width, not an automatically
  derived channel width. It is NOT implied by a prescribed discharge without
  depth, velocity and hydraulic constraints. It must be an explicit operational
  target or an independent hydraulic bound. Zero-width defaults avoid silently
  asserting it.
- No hydraulic contact, water delivery, undercutting or cliff survival follows
  merely because the energy budget exceeds this necessary lower bound.

## Observable consequence worth testing

At the same remaining debris amount (including zero), different post-abandonment
ice topographies can have different recovery requirements. Exported sediment
volume therefore need not be an adequate state variable for a recovery model.
That proposed environmental consequence needs Claude's verified spatial runs.
The included synthetic example only verifies the geometric calculation.

## Independent checks and use

The tests enumerate every subset of multiple 3x3 grids, independently check
connectivity by flood fill, and recover the exact minimum excavation. Other
tests check the finite-width relaxation, material scaling, level monotonicity,
walls, strict exclusion and invalid inputs.

From the parent folder:

    python3 -m unittest discover -s recovery_certificate -t . -p 'test_*.py' -v
    MPLCONFIGDIR=/private/tmp/ice_novelty_mpl python3 -m recovery_certificate.demo

Before applying to actual exports, export the ice-bed elevation, valid wall
mask, coordinate spacing and an independently defensible water-level upper
bound. Sweep the water-level bound rather than select one favorable value.
Use optimistic elevation uncertainty (`z_lower`) and maximum usable energy for
a robust exclusion; this package does not infer those uncertainty bounds.
