# Independently justified physics extension — bounded decision

## Decision

**GO** for an executable energy-adequacy gate and explicit independent input requirements. **NO-GO** for claiming a new successful melt predictor, historical originality, field superiority, or durable recovery on the present data. Six independent limiting/refusal tests pass. No Claude file or frozen result changed and no production campaign was launched.

## What the current no-go theorem actually assumes

For equal hydraulic rating coefficients and a positive toe-minus-bypass bed gap, bypass discharge is at least half the imposed supply. For common fixed water temperature and transfer coefficient, toe melting is bounded above by its film limit and bypass melting below by the finite-water rate at half supply. The resulting small positive upper bound cannot repair61.076mm in21days: its maximum is0.20064mm for the representative2Q0 pulse. This bound does not include other ice-directed heat, different local water temperatures or exchange coefficients, or additional geometry-changing processes.

The new `forcing_gate.py` generalizes the bound while preserving its logic. Write E_t <= max(0,h_t,max*T_t,max + q_ext,t,max)/(rho_i*Lf) and E_b >= max(0,q_water(Q/2,A_b,h_b,min,T_b,min)+q_ext,b,min)/(rho_i*Lf). Their difference bounds gap repair throughout the positive-gap interval. All supplied bounds must apply jointly over that interval. Passing this necessary-condition screen means only NOT EXCLUDED, never predicted success. Known extra heat can invalidate the restrictive old theorem without validating the broader model.

## Mechanism1: differential externally supplied ice heat

Spatial radiation, shading, terrain longwave and turbulent exchange are established controls on cliff melt; they are not a new mechanism. Buri2016 explicitly solves distributed energy exchange and moving cliff geometry. Kneib2022 compares the static energy-balance model without new calibration and reports changing contributions of radiation/turbulent exchange. Adding these parameters can improve realism but must be compared against those strong existing energy-balance models.

The representative screen requires at least10.276W/m² of additional TOE-ICE heat advantage to avoid this specific21day exclusion. This is a necessary lower limit, not an observed or recommended parameter. We must NOT set an extra20W/m² solely to make recovery happen. In particular, measured atmospheric heat on a vertical cliff is NOT the same quantity as heat delivered per horizontal toe-floor area. View factors, local surface normal, water/debris attenuation, thermal storage/cold debt and the geometry mapping are needed. Covered-water transport remains unsupported by the repaired native physical solver.

Required inputs: independently measured/projected shortwave and longwave, local normals and occlusion, albedo/emissivity, surface/debris temperatures, meteorological sensible/latent terms, actual receiving footprint, and a justified mapping from external heat to that ice. Required uncertainty: joint bounds on the DIFFERENCE between routes, not two independently chosen favorable marginal extremes. No such paired floor-heat forcing has been established for the current Swiss comparison.

## Mechanism2: local exchange/temperature asymmetry

A common upstream water temperature does not guarantee equal local bulk temperatures after different thermal paths, and a single h need not apply to different local depths/velocities. These are established finite-advection and heat-transfer effects. The new gate permits independent h/T bounds for each branch instead of assuming symmetry. For fixed h=200 and bypassT=.5, an increase of about.0514K in the toe temperature upper bound would remove this exclusion; alternatively the film-coefficient upper bound would need to increase by about20.55W/m²/K atT=.5. These only relax an upper bound; they do not demonstrate actual incision or recovery.

Required inputs: bulk temperatures at route entries/contact points, section geometry and sensor clearance, discharge at the target, distances and wetted perimeters, roughness/transfer-domain justification, and uncertainty in upstream cooling and route partition. A coefficient fitted to the same target melt event is NOT independent evidence. A free station-specific gain is unjustified when one calibration interval is below survey uncertainty. The existing project already implements finite heat supply; this extension must not be presented as the first such model.

## Scientific adequacy and identifiability

Melt identifies delivered energy per area, not independently h, bulkT, contact fraction, conduction loss and cold debt. Matching melt by tuning any of them does not identify the physical mechanism. For finite-water exchange, water-side temperature drop and flow plus known geometry can separately constrain delivered energy; local velocity/depth or independent flux measurements constrain h. Additional parameters without these observations make the inverse problem weaker.

A fair experiment would predeclare measured route-specific bounds, predict the sign and magnitude of gap drift, and compare against the existing common-h model AND a distributed energy-balance/evolving-topography baseline with matched input support. If the required inputs are absent, only sensitivity analysis is authorized and no field-superiority result follows. Numerical positivity and conservation are not enough to close that gap.

## Primary sources read

* Buri etal2016, *A physically based3-D model of ice cliff evolution over debris-covered glaciers*, https://doi.org/10.1002/2016JF004039 — distributed atmosphere and geometry/reburial precedent; not a baseline invented for this project.
* Kneib etal2022, *Sub-seasonal variability of supraglacial ice cliff melt rates and associated processes from time-lapse photogrammetry*, https://tc.copernicus.org/articles/16/4701/2022/ — section3.7 and results/discussion support the need for whole energy balance and geometric interactions; native melt data do not isolate our toe-floor port.
* Petersen etal2024, *Stream hydrology controls on ice cliff evolution and survival on debris-covered glaciers*, https://esurf.copernicus.org/articles/12/727/2024/ — stream undercutting/debris evacuation known; observations do not uniquely identify their relative contributions.

The report derives a new testable input requirement for this implementation. It does NOT claim a new physical law, new mathematical priority, STS success, or a Nature-level discovery.
