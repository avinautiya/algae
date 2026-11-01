# Spatial recovery-screening CLI

This self-contained tool evaluates a **conditional necessary energy condition**
for opening a route through an ice bed after perfect debris removal. It does not
solve hydraulics, forecast sustained ordinary-flow contact, or predict glacier
melt. It makes no claim of new shortest-path mathematics or external validation.

## Run

    python3 -m unittest -v test_spatial_forecast.py
    python3 spatial_forecast.py example_input.json --output example_output.json
    MPLCONFIGDIR=/private/tmp/ice_spatial_mpl python3 spatial_forecast.py example_input.json --output example_output.json --plot example_plot.png

Requires Python 3 and NumPy; plotting additionally requires Matplotlib. There
are no imports from Claude's active solver or any other research package.

## Schema version 1

`example_input.json` documents the complete shape. Required fields are:

- `schema_version`: integer 1.
- `elevation_m`: rectangular finite numeric matrix of post-clearance **ice**
  elevations. First row is inlet; last row is the target.
- `dx_m`, `dy_m`: positive streamwise and transverse cell widths.
- `level_bounds_m`: `lower`, `upper` parameter bounds on a genuine **global local
  free-surface upper bound**, never merely upstream mean stage.
- `usable_heat_bounds_J`: nonnegative `lower`, `upper` bounds on all usable heat
  over the forecast period, including stored heat and all external sources.
- `assumptions`: every documented prerequisite explicitly asserted true.

Optional fields:

- `elevation_bounds_m`: elementwise `lower`, `upper` matrices enclosing nominal
  elevations. Bounds must apply to the reconstructed subcell terrain, not just
  noisy cell-centre points. Defaults to exact nominal matrix.
- `minimum_cross_section_width_m`: required projected transverse route width
  at every streamwise row. Defaults to zero, which provides no nonzero continuum
  energy exclusion. Do not infer a physical channel width from discharge alone.
- `blocked`: Boolean matrix of truly immutable inaccessible cells. Mobile or
  porous debris is not an immutable wall.
- `ice_density_bounds_kg_m3`, `latent_heat_bounds_J_kg`: positive `lower`, `upper`
  material bounds; defaults 917 kg/m³ and 334000 J/kg.
- `resolved_cell_route`: false by default. True explicitly requires a route of
  entirely opened face-connected cells. This is a **discrete feasible class**.
- `provenance`: JSON metadata object carried into the output.

Unknown fields, duplicated JSON keys, NaN/Infinity, ragged arrays, numeric bools,
reversed uncertainty intervals, non-Boolean walls and numeric overflow are
rejected. An invalid CLI input returns code 2 and a JSON error without traceback.
Assumption flags are user assertions, **not automated physical verification**.

## Exact meaning of the calculation

At water-level bound H, cellwise minimum lowering is max(z-H,0). Within each
row, the default algorithm fills the cheapest transverse pieces until the
declared width is reached, allowing a partial final cell and disconnected
pieces. Summing these relaxed row costs underestimates the volume required by
any actual corridor meeting that width. It deliberately drops connectivity,
positive water depth, friction and heat-allocation constraints.

The whole-cell path cost is also returned as a **diagnostic**, using Dijkstra
with nonnegative cell-volume weights. It is not included in the default energy
bound: a narrow continuous trench need not melt an entire cell. Only the
explicit `resolved_cell_route=true` mode uses max(width cost, path cost).

## Deterministic uncertainty enclosure

For the declared independent parameter box, lower terrain, higher permitted
water level, lower density and lower latent heat minimise the required energy
bound. Their opposite endpoints maximise it. The rowwise and path minima are
monotone in each nonnegative excavation cost, so these endpoints enclose the
**computed bound**, not actual required heat or actual flow.

- `excluded_for_all_declared_box_values`: even the optimistic required-energy
  lower bound exceeds the maximum usable heat. This is a conditional exclusion.
- `exclusion_depends_on_declared_values`: the screen changes classification
  within the supplied box. It is not a probability.
- `not_excluded_for_any_declared_box_value`: the necessary condition excludes
  none of the supplied values. It does **not** establish recoverability.

Equality is inconclusive. Infinite infeasibility is represented using standard
JSON `{ "value": null, "infinite": true }`, never nonstandard Infinity. Output
includes canonical-input and exact-source SHA256 hashes for reproducibility.

## Physical failure modes

The bound is invalid if mechanical collapse, erosion, uplift/displacement,
unmodelled lateral/englacial paths or changing boundaries open a route outside
the stated melting-only feasible class. Velocity head, run-up and backwater can
make local water levels exceed inlet stage; a global bound must be justified
independently. Width requirements and cellwise terrain bounds need physical
support. Subgrid channels cannot be excluded merely because a discrete model
does not represent them. Positive-flow or sustained-contact targets need further
hydraulic modelling and are not supplied here.

The example and figure are **synthetic method checks**, not final experimental
data or a validated glacier discovery. Twelve tests cover narrow trenches,
uncertainty, walls, zero heat, equality, invalid inputs and randomized monotone
enclosures. Existing frozen campaigns and active solvers remain untouched.
