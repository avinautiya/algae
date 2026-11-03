# V2 cold-ice interface basis and independent review

## V1 remains a no-go

The original frozen no-refreezing conversion was not applicable to the forward-fit temperature result. `V1_PHYSICS_NOGO_INDEPENDENT.json` hashes the untouched normalized artifact and independently records all six forward geometry×hold failures. Negative/positive heat in its half-hour traces is0.520–0.591%; the six derivative-OLS scenarios remain eligible. The artifact states camera outcome fields were not loaded and a melt gain was not fitted. No camera outcomes were opened by this review. The small negative fraction is not permission to clip it or retroactively pass V1.

V2 is a distinct physical specification. Its construction must precede opening either its training or held-out camera outcomes. V1 outputs and its failed gate remain preserved.

## The physical balance

Use a common unit plan area, with z downward, debris above an ice interface and an explicit cold-ice layer below. Let qd enter the interface from debris and qi leave into the ice. With positive melt mass flux m in kg m⁻² s⁻¹,

qd−qi=L m,

TΓ≤0, m≥0, TΓm=0.

Each material obeys sensible heat conduction Cvol∂tT=∂z(k∂zT). Ice initial temperatures and its deep boundary must be≤0. The surface forcing may be signed. When net interface heat is positive at0°C it produces melt; when it is negative, the interface cools and cold content develops instead of creating hypothetical refreezing water.

**No measured water availability is needed to implement this deliberately no-refreezing branch:** it never consumes liquid water. Meltwater is assumed to drain at0°C. This does not establish that the real debris/ice system has no retained water or refreezing. Wet-pore phase heat, capillary retention and liquid refreezing remain outside the model. Measured soil moisture makes that omission an environmental uncertainty, not an accounted-for process.

## Conservative backward-Euler active set

For centre-to-interface half-cell conductances Gd=2kd/Δzd and Gi=2ki/Δzi,

qd=Gd(Td−TΓ), qi=Gi(TΓ−Ti).

1. Solve both material domains implicitly with common trial TΓ=0. If qd−qi≥0, retain those solutions and set melt mass over the step ΔM=Δt(qd−qi)/L.
2. Otherwise set ΔM=0 and **resolve the entire coupled system** using the series-resistance interface link G=GdGi/(Gd+Gi). Recover TΓ=(GdTd+GiTi)/(Gd+Gi) from the newly coupled cell states.

The second step cannot be replaced by inserting a harmonic interface temperature into the first branch's unchanged cell solutions. Temperature changes feed back into both implicit domains.

For positive capacities/conductances, elimination of the bulk cells gives interface imbalance S(TΓ)=A−B TΓ with B>0. The active trial tests A. For A≥0 the melting solution TΓ=0, ΔM=ΔtA/L satisfies complementarity. For A<0 the coupled no-melt root TΓ=A/B<0 does. Thus the linear discrete problem has a unique admissible active branch, including a coincident zero-flux boundary case. This is a standard complementarity argument, **not a new mathematical-priority claim**.

Summing the finite-volume rows gives the required ledger:

ΔEdebris+ΔEice+LΔM=Δt qtop−Δt qdeep.

A negative deep-outward flux is heat entering from the fixed deep boundary and must retain its sign. Negative sensible storage is cold content relative to0°C, not negative mass. kg/m² is numerically mm water equivalent. Ice-interface retreat speed would be m/ρice; the current calculation does not advance that interface.

## Independent inspection of the implemented module

I inspected the new `model/cold_ice_column.py` without modifying it. It implements the full coupled branch, correct half-cell resistance, signed deep flux, sensible+latent ledger, nonpositive ice-temperature check and fixed-geometry scope. No inspected active-set, unit or ledger defect was found.

`test_cold_interface_basis.py` adds four independent tests:

-100 random one-step cases against an explicit two-cell Schur-complement oracle that does not use Thomas elimination, covering both branches;
-exact cold steady series resistance and zero melt;
-warm steady50W/m² heat supply converted to latent mass without an extra density/1000 factor;
-cooling from0°C produces sensible ice cold content and zero melt, with energy closure.

All four pass. These verify the declared numerical component; they do not validate glacier melt accuracy.

## Required scientific gates before interpreting a field prediction

1. **Initial ice cold content:** freeze the initial profile using independent evidence or an explicit assumption family. No evaluation melt may select it. A0°C initial profile is a temperate assumption, not a measured Pirámide profile.
2. **Finite ice domain:** a deep0°C boundary is an energy reservoir. Demonstrate insensitivity to its depth over the full interval or carry its influence as a model sensitivity. The diffusion penetration scale grows as√(αice t); a2m depth is not automatically semi-infinite over weeks.
3. **Retreat approximation:** this is not a moving-front Stefan solver. Fixed geometry neglects coordinate-motion advection and changing debris/sensor distances. Report cumulative retreat relative to geometry and a moving-interface Péclet scale vℓ/αice. If those are not small, the raw component is an approximation whose discrepancy needs testing, not a complete physical melt budget for the evolving glacier.
4. **Wet debris:** sensible-only debris cannot claim to account for freezing pore water, latent retention or moisture-dependent transport. Explicitly retain this limitation or implement and verify a separate water/enthalpy state before claiming it is solved.
5. **Gain interpretation:** coupling fixed known ice properties breaks the old normalized-heat linear scaling with debris C. Changing debris C while keeping α fixed changes k and relative effusivity, so a physical-C trial requires resolving the coupled system. A subsequent training-only output gain may remain an empirical bias correction, but the verified raw physical ledger must never be multiplied by that gain and called a conserved physical mass solution.
6. **Comparators:** use the same new interface physics and fixed assumptions for forward-fit and derivative-OLS thermal candidates; retain the same simple baselines and predeclared primary. This isolates parameter-inference differences from an unequal physics advantage.
7. **New freeze:** pin V2 code, assumptions, all candidate geometries/holds, gain rule and outcomes access order. Do not rewrite V1's failed eligibility result.

## Further verification requirements

Before a long field run: homogeneous-material composite-column equivalence in the entirely cold branch; spatial/time refinement against a smooth analytical diffusion mode; sudden cooling/rewarming active-set switches with stepwise/global ledgers; full-domain depth refinement; parameter/provenance refusal tests; and mass/energy consistency across any later empirical-gain reporting. A genuine moving-front Stefan benchmark would test a different solver and cannot certify this fixed-interface approximation merely by analogy.
