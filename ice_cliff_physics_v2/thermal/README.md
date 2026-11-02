# Conservative wet-debris thermal contacts and independent slab reference

This module fills a specific missing boundary condition: finite-water heat exchange with a debris-covered temperate ice patch and a parallel bare-ice patch. It does not establish field accuracy, historical novelty, a complete ice-cliff model, refreezing, collapse or free-surface CFD.

## Production-sized one-node kernel

`wet_debris.advance(ThermalState(Td), ThermalParameters(C,K,Gb,Gd), dt, Tin, F, atmospheric_power_W=Qa)` accepts extensive units:

- C: debris heat capacity, J/K.
- K: debris-to-ice conductance, W/K.
- Gb/Gd: water contact conductances for bare/debris patches, W/K.
- F = rho_water cp_water Q: advective water capacity, W/K.
- Qa: atmospheric power into the dry debris surface, W; signed.

All ice boundaries are at 0°C. C and conductances are nonnegative; time is positive. Cold substrates and trajectories crossing below 0°C are explicitly refused. C=0 is an algebraic steady limit. F=0 means the contact is dry and water exchange is zero. Atmospheric power on bare ice is a separate caller-owned port.

For F>0, G=Gb+Gd, alpha=Gd/G, phi=(1-exp(-G/F))/(G/F). The mean stream temperature along the contact is

    Tbar_water = phi Tin + (1-phi) alpha Td.

Consequently water-to-debris power is a Tin - b Td, where

    a = Gd phi,
    b = Gd [Gb/G + phi Gd/G].

Bare-ice power is Gb phi Tin + Gb(1-phi)alpha Td. The single temperature reservoir obeys

    C dTd/dt = Qa + a Tin - (b+K) Td.

The kernel integrates this linear ODE and all heat inventories analytically. Signed water-to-debris heat is retained: hot debris can warm the outgoing water and/or transfer stored heat to the bare contact. The closed budget is

    Qa dt + E_water,in - E_water,out
       = Delta(C Td) + E_bare_ice + E_covered_ice.

Interface heat is returned, not silently converted to a geometric volume or used twice. Callers with finite cold debt must not call this temperate boundary and then claim it resolved a cold-ice temperature field. Area changes require explicit remapping of enthalpy and incoming material.

The nonnegative-temperature solution gives a conditional passivity bound:

    E_ice <= C Td_initial + E_water,in + max(Qa,0) dt.

This is an elementary energy consequence of this defined closure, not a claimed new theorem. Bare/debris patches are lumped over a shared quasi-steady stream; no vertical mixing or local stream profile is resolved. At vanishing positive throughflow, stored debris heat may still redistribute to a bare patch through the shared-water closure. Setting F=0 declares that water absent; this is a dry/wet modelling switch, not a validated smooth contact transition.

## Why the resolved reference is necessary

Exact ODE integration and conservation do not validate a uniform-temperature debris approximation. `resolved_slab.py` independently solves one-dimensional debris diffusion using cell-centred finite volumes, implicit Euler and a direct tridiagonal solve. It imports no lumped solver routine. The top surface is algebraically coupled to the same declared finite-water boundary; the bottom ice is fixed at 0°C. Reported interface flux is the same discrete flux removed from the debris cells. Spatial discretisation is second order on the tested smooth solutions; time discretisation is first order.

Use either:

    SlabParameters(thickness_m, conductivity_W_m_K,
                   volumetric_capacity_J_m3_K, area_m2)

or

    ArealSlabParameters(capacity_J_m2_K, conductance_W_m2_K, area_m2).

The second form uses xi=x/d. Its C_A=rhoc*d and K_A=k/d fully specify temperature/heat evolution without inventing separate physical thickness or material conductivity. The diffusion time is C_A/K_A. It cannot identify d, k and rhoc separately.

State is `SlabState(tuple(cell_temperatures_C))`; all are nonnegative. The call is

    resolved_slab.advance(state, parameters, dt_s, Tin, F,
                         bare_water_conductance_W_K=Gb,
                         debris_water_conductance_W_K=Gd,
                         atmospheric_power_W=Qa).

It returns a dictionary with state, surface temperature, bare/covered ice heat, signed water-to-debris heat, storage change and separate discrete/stream-law residuals. In areal coordinates the energy of cell i is C_A*A/N*Ti. Changing A requires explicit remapping of every cell's extensive energy; it must not preserve J/m² on newly created area without an incoming-enthalpy term.

The reference is independently checked against steady Robin flux, a sine-series diffusion limit, hot-slab analytic heat release, space/time refinement and energy accounting. It still shares the prescribed water-contact law and fixed 0°C ice boundary. It does not verify that either closure describes a natural boulder deposit.

## Execute

    python -m unittest discover -s . -p 'test_*.py'
    python verify_thermal.py
    python run_slab_comparison.py

Tests and numerical calculations use the Python standard library. The optional comparison PNG/PDF uses Matplotlib. Results are synthetic implementation/closure tests, not field scores. Do not fit an effective C to these reference runs and call the fitted one-node approximation an independent prediction.

## Separate wet and dry normal columns

`column_partition.py` keeps distinct wet and dry normal-cell temperature profiles. Its `ColumnPartition(C_A, wet_area, dry_area, wet_profile, dry_profile)` stores each class's extensive enthalpy as C_A*A/N*sum(Ti). `remap(partition, total_area_new, wet_area_new, incoming_profile)` first imports any new area into the dry class, then transfers only material whose classification changes. The returned result reports imported energy, internal wetting/drying transfers and every normal cell's conservation residual. Empty classes have None profiles. Total-area contraction, negative/cold profiles and grid-size changes are refused.

Heat the wet class with its own finite-water contact and no dry atmospheric power; heat the dry class with F=0 and dry atmospheric power. This prevents instantaneous averaging of hot dry debris with cold wet debris when the classification is unchanged. Class conversion uses the donor's average normal profile. It does not resolve the height-wise correlation between temperatures and which particular strip the moving wet edge crosses. New area entering the dry pool first is an explicit material convention, not a universal geometric fact. There is no lateral conductive exchange between the two classes.

Ten partition tests include unchanged-split profile preservation, per-cell wetting conservation, all-wet/all-dry cycling, hot/cold incoming material, empty-domain creation, refused contraction and a coupled finite-water/diffusion/remapping energy budget across several boundary changes. These ten partition checks supplement twelve one-node and eight finite-volume checks.

The partition additionally carries optional `area_reference_m2`, preserving the exact total requested by geometry even when separately rounded wet/dry areas sum one ULP differently. Class areas must agree within a four-ULP representation check. This is not a contraction tolerance: a caller requesting a genuinely smaller canonical total, even by one ULP, is refused. Material/temperature updates using dataclasses.replace retain the reference.

## Classical analytic reference

`robin_reference.py` supplies an independent zero-initial constant-coefficient Robin-slab eigenfunction reference, including stable exact integrated-series identities. See `ROBIN_REFERENCE.md` for the derivation, generic parallel-contact mapping, APIs and limitations. Ten additional tests verify it against the PDE, boundary conditions and independent finite-volume convergence. It must not be reset to zero when benchmarking an inherited nonzero thermal state. Total current tests: 40.
