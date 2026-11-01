# Parallel transmission contract

`transmission.py` uses metres, seconds, and cubic metres per second. Supply
`Route(name, floor_z_m, debris_height_m, width_m,
weir_coefficient_sqrt_m_per_s, ...)` explicitly. The weir coefficient has units
m^(1/2)/s and is an assumed or calibrated closure, not a measured default.
An aperture also requires an explicit nonzero dimensionless discharge coefficient;
its default zero disables it. Alternatively specify whole-route Darcy conductance
in m²/s (including area and path length). Do not use both under-rock closures.

The actual floor bounds every path. Overtopping begins above floor plus debris
height; under-rock transmission begins only above the actual floor. Receiver head
can suppress flow. No ice-floor penetration or englacial route is introduced.
The submerged overtopping expression is a monotone effective rating, not a
resolved submerged-weir solution. Aperture and overtopping addition assumes
vertically distinct paths within one corridor; it does not resolve their interaction.

Each Route is a laterally disjoint corridor. Partition the available corridor
width among Route objects once; do not assign the entire valley width to each.
The aperture and overtop width in a single Route refer to vertically distinct
openings, not two lateral copies of the corridor.

For the aperture, hydrostatic upstream-minus-downstream pressure at elevation
z gives velocity `Cd sqrt(2g max(0, eta-max(eta_d,z)))`. Thus discharge is
`Cd w integral sqrt(2g max(0, eta-max(eta_d,z))) dz`, from the actual floor
to `min(floor+aperture_height, eta)`. Below downstream stage the head is
constant; above it, integrating `sqrt(eta-z)` gives the 2/3 power-law primitive.
The code evaluates the latter with `log1p`/`expm1` to preserve tiny-aperture
precision. This is an assumed hydrostatic orifice closure, not resolved flow.
It represents a thin opening. A long under-rock gap also has viscous/frictional
head losses, which this coefficient cannot automatically represent. Tiny-gap
flows must not be presented as credible field predictions: aperture Reynolds
number and path-length-to-aperture ratio matter. The scan reports a Reynolds
indicator using assumed kinematic viscosity 1.79e-6 m²/s; no high-Re validation
or low-Re extrapolation is established. The Darcy option is a separate effective
conductance requiring independent justification, not a cure for unknown geometry.

`solve_common_stage(routes, inflow_m3_s, downstream_stage_m, ...)` requires at
least two uniquely named routes and returns stage, each route discharge, residual,
and iteration count. Scalar bisection conserves prescribed steady inflow to the
requested absolute flow tolerance. `max_stage_m` supplies a physical upper bound;
unreachable discharge raises an error. Without it, bracketing expands upward.
These ratings are nonnegative and monotone in upstream stage, with zero flow
below the floor or receiving stage. They are not bounded in stage unless the
caller supplies that bound; individual flow is bounded by total inflow at solution
up to numerical tolerance. This is a steady rating allocation, not transient CFD.

`hydrostatic_diagnostics` returns depth/length, u²/(gR), Froude number and signed
bed slope. Unknown curvature stays `None` with a warning. Indicators and warnings
are for review only: no universal or guessed validation threshold is imposed.
The caller must assess vertical accelerations, abrupt geometry, nonhydrostatic
effects, and environmental applicability independently.

Run: `cd transmission && /usr/bin/python3 -m unittest -v`.

`/usr/bin/python3 synthetic_scan.py` writes a synthetic JSON table and SVG.
The aperture grid is a mathematical sensitivity scan, NOT a defensible site
range. A serial plug-flow heat closure assumes a toe exchanger with UA=200 W/K,
water initially 0.5 K above melting, and a melting-temperature boundary. It
reports modeled toe-floor uptake `rho cp Q Tw (1-exp(-UA/(rho cp Q)))` and
remaining advected heat separately. The toe floor itself is this exchanger.
Residual heat is an upper bound available downstream; actual downstream uptake
would require its own heat-transfer model.
Neither this assumed closure nor a positive pathway establishes field relevance.
