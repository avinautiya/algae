# One decisive CFD test before a larger campaign

Do not alter an active Claude job or assume this package has already established a field discovery. This specification is an independent experiment proposal for the existing flow–heat–moving-ice development.

## Question

Can a finite temporary insulating/obstructing deposit cause a hydraulic separation that survives its complete removal? Does an open corridor around the same solid volume prevent that separation?

## First milestone deliberately excludes granular-motion closure

Prescribe immobile coarse obstacles before removal, then export all of them at a declared time. Record their mass and enthalpy as an external export. This isolates thermal/topographic memory without inventing a packed-bed transport coefficient. It does not predict whether a natural flood would actually remove them.

Run the same forcing and initial geometry with early and late removal. Include a bare control, a non-melting substrate control, and at least two layouts of the same grain volume: one allowing an unbroken bare corridor and one obstructing the full cross-stream route. Publish all layouts, including a layout that fails to produce lockout.

The reduced model's illustrative scales are Q=0.07 m³/s, water superheat=0.5 K, width=0.5 m, reach length=2 m, initial toe ice floor=-0.08 m relative to bypass floor=0, initial coarse solid stock about 0.096 m³ and 0.2 m grain size. Exact grain counts and geometry must be integer-realizable; adjust the reduced comparator to the actual CFD geometry rather than claim a fractional-particle match.

## No fitted success criterion

After complete removal, return to identical ordinary upstream discharge. Measure actual flow to the cliff-side route, crest elevations, wet ice area, integrated toe and bypass melt, and outlet enthalpy. A failure of ordinary-flow access must be established from the resolved free surface and available paths, not assigned when mean stock crosses a threshold.

Compute the activation prediction from independently extracted control geometry. The quantitative reduced forecast should be frozen before reading the final comparison runs. Development runs may determine mesh requirements and expose a bad closure; changes after seeing results require a new version and new tests.

## Necessary numerical and physical checks

- Separate grid and time convergence for the hydraulic split, ice-bed relief, timing and integrated melt.
- Open-path/bare control, equal material properties, verified water mass and sensible-plus-latent energy budgets, and conservative obstacle enthalpy export.
- Distinguish steady-flow approximation error from interface time stepping. Include changing width or show the observed width change is negligible.
- Establish whether lateral heat around grains defeats the area-weighted shield approximation.
- Vary additional heat to the inactive floor. Existing reduced tests show that the effect can disappear when this source becomes large enough.
- The exclusion certificate must use a bound on **all** toe-lowering energy, including initial grain heat and other external fluxes, and must allow any bypass rise that the solver permits.

## Outcomes and claims

If the resolved system reconnects through a bare corridor, the mean-crest model is falsified for that layout. If a blocked layout alone locks out, the result is geometry-dependent and should be reported as such. If late clearing repeatedly fails after the grains are gone, and the prediction follows measured ice relief and survives the controls, that supports the specific reactivation-memory contribution.

The output would improve a stream-access closure for local glacier melt modeling. It would not, by itself, establish improved regional Himalayan glacier forecasts or field-calibrated cliff lifetimes.

## Added topology and permeability gate

A solid transverse strip in `topology.py` is a nonporous limiting geometry, not a default assumption about boulders. Test open corridors, through-debris flow, and under-rock gaps explicitly. Start with solid patches to verify the theorem's limiting conditions, then test porous/grain-resolved arrangements. A generic impermeable occupancy mask would manufacture blockage and is unacceptable as environmental validation.

The supported environmental gap is quantitative prediction of access/recovery from connectivity and evolving ice relief, rather than a claim that preferential thermal channels, heterogeneous debris, or stream–cliff feedback were previously unknown. See TOPOLOGY_PRIORITY_CHECK.md for primary precedents and limitations.
