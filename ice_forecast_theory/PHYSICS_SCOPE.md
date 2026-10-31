# Environmental physics and scope ledger

No finite research model accounts for every factor. A credible model states which processes it solves, which it bounds and which exclude an application. This ledger is an audit to carry into the later physical study, not a claim that all rows are implemented.

## Newly quantified issue: geometry and freshwater buoyancy

For near-freezing fresh water, warmer water below the density maximum can be denser than meltwater at 0 C. With z increasing upward,

```
N^2 = -(g/rho) d rho/dz
Ri = N^2 h^2/U^2
Ra_destabilizing = -N^2 h^4/(nu kappa)
```

Cold, light water at an ice **bed** beneath warmer, denser water gives a negative mean N^2; swapping bed and ceiling reverses the endpoint-density stratification. This sign matters when transferring results between an ice ceiling and an open channel above ice. It is established physics and appears explicitly in prior work, including Couston et al. (2021).

`regime.py` uses the declared UNESCO pure-water polynomial for a transparent screening calculation. It is not a replacement for a validated EOS/property model in production, especially at nonzero salinity, different pressure or across the density maximum. Its endpoint-density approximation can miss a nonmonotone interior profile. The software deliberately rejects unsupported salinity inputs.

Bulk velocity can conceal the relevance of buoyancy in slow recirculating water. For a **synthetic** 0.05 m water layer at 0.2 C, changing a supplied characteristic speed from 1 m/s to 0.005 m/s increases |Ri| by 40000. This does not prove such a recirculation exists or predicts a convection threshold. A local speed, vertical thermal contrast and residence time must come from resolved results or observations.

Do not apply the classical 1708 Rayleigh threshold or a universal Richardson cutoff to a sloping, sheared, open, melting jump without deriving the corresponding stability problem.

## Process ledger

| Process | Current foundation / required treatment | Scope consequence |
|---|---|---|
| Incompressible liquid motion | CFD verification is owned by Claude's active task | No claim based only on depth-averaged reference agreement |
| Air–water interface and jumps | Verify mass flux, interface transport and appropriate pressures | Rigid-lid results do not establish free-surface behavior |
| Temperature and wall heat flux | Verify the actual liquid thermal equation and boundary fluxes | Temperature-profile agreement alone is insufficient |
| Melting and phase-change mass | One conservative exchange; distinguish interface/material velocities | Equal-density reference does not test the density jump |
| Freshwater buoyancy | Screen sign, local/bulk ratios; assess in physical cases | Constant positive expansion coefficient may be inappropriate near freezing |
| Salinity/dissolved solids | Use suitable EOS/transport or exclude | Pure-water screening cannot model brackish water |
| Temperature-dependent viscosity, conductivity, cp | Declare properties and sensitivity | A nominal laminar benchmark does not establish environmental similarity |
| Turbulent momentum/heat transport | Selected 3D closure/resolution comparisons | RANS, LES and DNS have different evidentiary limits |
| Near-wall roughness and separation | Geometry/roughness uncertainty, not only an effective Manning coefficient | Heat and drag need not respond identically |
| Surface tension and contact line | Assess local curvature and thin-film scales | Long bedforms alone do not justify removing capillarity |
| Ice conduction/cold content | Temperate base case; add when relevant with initial/boundary data | Warm-water heat is not necessarily all latent heat |
| Atmospheric and radiative exchange | Measured or modeled net heat, surface definition and sign | An imposed power is not a validated surface-energy model |
| Inlet enthalpy and temperature definition | Flux-weighted temperature and observation operator | Near-bed point sensors are not automatically bulk temperatures |
| Water/ice thermal storage | Retain transient inventory where relevant | A steady heat balance is not a universal instantaneous melt bound |
| Mechanical energy and dissipation | Track retained heating, kinetic/potential/pressure work | Numerical loss is not automatically physical latent heating |
| Lateral wall erosion/width | Evolve and verify, or state fixed width | Do not forecast cross-sectional change from a bed-only model |
| Wetting, drying and control movement | Convergent event/interface treatment | Smooth tangent gains may fail at topology changes |
| Finite-reach boundaries | Specify causal hydraulic influence and sufficient thermal development | A periodic heated domain is not an open reach |
| Lake/network storage | Existing conservative coupling plus compatible hydraulic boundaries | Imposed discharge is forcing, not a predicted lake release |
| Ice creep/deformation | Compare deformation and forecast timescales, then model or exclude | Short-horizon neglect must be justified |
| Fracture/englacial drainage | Separate mechanics and verification if present | A surface channel model does not predict hydrofracture |
| Debris/sediment/slush | Exclude clean-water study or add verified transport/thermal physics | Roughness alone is not sediment transport |
| Weathering crust/permeable ice | Assess lateral inputs, seepage and changing contact geometry | Solid impermeable bed assumptions may fail |
| Geometry and measurement uncertainty | Initial/boundary ensembles and observation operators | Fine numerical resolution does not create fine field knowledge |
| Seasonal/glacier-scale evolution | Separate coupling, scale and conservation tests | A channel improvement is not a demonstrated regional melt improvement |

## Quantities computed by the offline screen

Re, Pe, Pr, Fr, We, Bo, sensible/latent Stefan number; mean thermal N^2 and signed Ri/Ra; transit, vertical diffusion and buoyancy times; advected sensible power, potential-power scale and specified external power. A listed-power latent-volume equivalent is only a unit-consistent conversion, not a universal melt ceiling. Stored heat and other mechanical/thermal terms need the full control-volume ledger.

These diagnostics prioritize controlled sensitivity checks. They do not assign an automatic valid/invalid label to the whole physical model or replace the checks in the active milestone.
