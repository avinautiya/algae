# Zero-external-flow atmospheric melt

The earlier coupled model computes atmospheric phase water and then compares it
with external inflow. At zero external inflow this ratio becomes infinite, so a
physically possible atmospheric-melt-only case is refused. That refusal honestly
marks missing hydraulic-feedback scope; replacing infinity with zero would hide
the missing physics and is not a correct fix.

## Small, defensible extension

Permit a separately named **one-way free-drainage boundary mode**, only when the
caller supplies evidence/assumptions that local meltwater drains below all modeled
thermal-contact patches without pooling, backwater or an upstream hydraulic path.
Under that contract:

1. At Q_external0 the prescribed upstream routing flows remain0. Atmospheric
   energy still warms/melts the appropriate ice/debris interfaces.
2. The generated water mass is `rho_ice * V_melt`; its outgoing liquid volume is
   that mass divided by rho_water. It leaves at the melting point, with zero
   sensible-water heat referenced to that temperature.
3. The outlet mass ledger is `V_out = V_external + V_phase`. Phase water must not
   be sent backward through the obstruction/common-stage solve.
4. No additional melt is charged to the0-C generated water. Do not recycle latent
   heat into a second melt calculation.
5. Report `phase_water_m3`, source-specific generated discharge, outlet capacity
   margin and mode. The phase/external ratio is undefined at zero external flow,
   not zero and not evidence of a negligible feedback.

An explicit independently justified downstream capacity can be used as a
necessary check, but it does not by itself prove the absence of pooling or
thermal-contact changes. The mode remains a boundary approximation, not a solved
runoff film or wetting/drying model. A dry model face exposed to net interface
energy is compatible with thin meltwater drainage only if the supplied flux
already represents that interface condition; station shortwave alone is not it.

## Refusals that should remain

Retain refusal if meltwater can collect on the modeled floor, wet a face that
the atmospheric model treated as exposed, reach the upstream pool, or exceed the
declared outlet domain. Those cases require at least local storage/routing and
an updated wet/dry observation operator. There is no general conservation-only
argument that removes their hydraulic feedback. For positive external flow the
earlier small-feedback bound should remain unless this stronger boundary contract
is explicitly adopted and its domain checked.

Tests for the parent extension should include zero-inflow positive atmospheric
melt, zero-energy zero-flow, complete liquid/ice mass balance, no invented water
heat, a declared pooling/capacity refusal, and unchanged legacy-mode refusal.
