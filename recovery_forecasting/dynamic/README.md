# Dynamic reduced forecast

This is an independent baseline, not a replacement for Claude's spatial solver.
It adds finite water AND sensible-heat storage to the previous ideal two-route
forecast. It may predict a supplied time-dependent hydrograph, whereas the
static scalar water-budget tool only tests existence of a feasible pulse.

## Equations and exact accounting contract

State: reservoir volume V and sensible heat U, and two ice-control elevations z.
Fixed reservoir hypsometry: eta = pond_floor + V/A. Temperature is superheat
above the melting point: T = U/(rho_w cp V). There is no ice inside this fixed
reservoir footprint; the melting outlet controls are downstream of it.

    dV/dt = I - Q1 - Q2
    dU/dt = rho_w cp [I*T_in - (Q1+Q2)*T]
    Qj = width_j sqrt(g) (2/3)^(3/2) [eta-zj-debris_j]_+^(3/2)
    Pj = rho_w cp Qj T [1-exp(-h_eff Aice_j/(rho_w cp Qj))]
    dzj/dt = -Pj/(rho_i L Aice_j)
    1/h_eff = 1/h_water + debris_thickness/k_debris

The zero-transfer and dry-route limits are handled explicitly. Outlet sensible
power is rho_w cp Qj T - Pj. Meltwater Pj/(rho_w L) enters downstream at the
melting point, diluting the reported downstream temperature; it does not return
to the reservoir or the control section. Upstream and downstream hydraulic
feedback from that small added flow is outside this reduced approximation.

The ODE tracks inlet/outlet water, inlet/downstream sensible energy and each
route's latent energy. Reservoir water closure, total sensible+latent energy
closure and latent-volume consistency must pass before returning a forecast.
These algebraic budgets are necessary verification, not independent validation.
The implicit Radau integrator handles the separation between fast reservoir
adjustment and slow ice lowering; it does not require resolving every hydraulic
relaxation time throughout a multi-day reduced forecast.

## Important restrictions

- Ideal hydrostatic critical controls with no entrance loss, backwater coupling,
  channel momentum storage or arbitrary lateral flow. Earlier partition
  discrepancies mean this model cannot independently validate spatial routing.
- A uniform deposit changes the crest and effective resistance. Removal is an
  exact prescribed event; grain supply, porous gaps and sediment transport are
  not predicted. Patchy coverage must not be represented by this uniform law
  without a separate justification.
- Fixed reservoir area and floor, separate melting outlet footprints, passive
  well-mixed water, temperate ice. No cold ice, atmosphere, frictional heating,
  mechanical collapse, moving cliff face or ice creep.
- If a melting control reaches the fixed reservoir floor, the run refuses. It
  does not introduce a hidden sill to keep the model running.
- Time-varying forcing and prescribed clearing are integrated in exactly aligned
  intervals. Output includes those boundaries, not just the regular sample grid.
- `toe_fraction` is route-0 reservoir outflow share. A contact/recovery threshold
  and persistence duration must be declared separately. No cliff-survival or
  validated environmental probability is emitted.

Run from the delivery folder:

    python3 -m unittest discover -s dynamic -t . -p 'test_*.py' -v
    python3 -m dynamic.cli dynamic/sample_input.json --output /tmp/dynamic_forecast.json

The example and demonstration figure are explicitly synthetic. Do not quote
them as a glacier measurement, a superiority benchmark or an originality claim.
