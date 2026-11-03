# Independent review of the additive thermal-memory component

Reviewed solver SHA256 `777e7c1799308bc9245fddb8ce30d04e44c6fe3f805a7d20b5017098659b4539`. No solver defect was found in the inspected half-cell coefficients, heat-capacity units, backward-Euler ledger or finite-scenario envelope wording. The reviewed solver is a fixed-thickness conductive sub-debris component; it is not a cliff model or a validated field forecast.

## Discretization and units

For n cell centres separated by Δz=d/n, interior conductance is k/Δz and each Dirichlet boundary conductance is **2k/Δz**. The top/bottom centre lies Δz/2 from its boundary. Storage per area per temperature is **Cvol Δz**, with Cvol in J m⁻³ K⁻¹. A mass-specific heat capacity in J kg⁻¹ K⁻¹ must first be multiplied by the appropriate bulk density; it cannot be passed directly. Porosity/moisture-dependent effective heat capacity requires independent physical inputs and is not automatically supplied here.

The backward-Euler row has storage/Δt plus its adjacent conductances on the diagonal and negative interior conductances off diagonal. In the single-cell limit both boundaries remain: the diagonal contains storage/Δt+4k/Δz. Independent review checks this edge case explicitly.

Summing rows cancels internal fluxes. Thus

Cvol Δz Σ(Tnew−Told)=Δt(qsurface−qbasal)

with both fluxes evaluated at the same implicit new state. This is exactly the recorded discrete ledger. Negative surface supply is physically permissible when stored heat exits upward; basal melt can exceed concurrent supplied heat by releasing initial sensible storage. Suppressing that sign or comparing basal heat to incoming heat without the endpoint storage would be incorrect. All accepted boundary/initial temperatures are nonnegative and the basal interface is fixed at melting, so refreezing/cold-ice conditions are intentionally unsupported.

## Independent numerical checks

`thermal_memory_independent_review.py` writes `THERMAL_MEMORY_INDEPENDENT_RESULTS_V3.json` and leaves agent-owned model files untouched.

- Single cell, stationary linear profile: basal flux53.333333W/m²; maximum deviation7.1×10⁻¹⁵W/m²; step ledger exactly zero.
- Homogeneous-boundary sine mode: exact discrete eigenvalue is4α sin²(π/2n)/Δz² and exact backward-Euler amplitude factor is(1+Δt λ)⁻ᴺ. Maximum discrepancy8.2×10⁻¹⁵°C. This checks the operator against an analytic eigenpair rather than a second implementation of Thomas elimination.
- Continuous sine-mode decay at fixed1s timestep: maximum temperature errors0.004415,0.001146,0.000313°C at8,16,32cells. Spatial refinement approaches the continuum solution; the remaining error includes backward-Euler time truncation. This is not presented as an isolated exact spatial convergence order.
- Eight-day periodic demos close their integrated sensible-plus-basal energy balances to1.6×10⁻¹³ relative or better.

## A physically meaningful synthetic demonstration

All forcing and parameters are declared synthetic, independent of every Arolla target. Conductivity0.96W/mK and volumetric heat capacity1.5MJ/m³K are dimensional illustrations, **not a provenance claim about a particular site**. The top boundary is an8°C-amplitude daily oscillation about10°C, held hourly; the initial profile is the linear10°C-to0°C profile. Both cases use40cells and60s stepping. Harmonics are evaluated after6days.

The analytic continuum transfer function from surface temperature to basal heat flux is

H(ω)=k λ/sinh(λd),  λ=√(iω Cvol/k).

|Debris thickness|Basal harmonic amplitude|Quasistatic amplitude k/d×surface amplitude|Phase lag modulo24h|Analytic transfer discrepancy|
|---|---:|---:|---:|---:|
|0.18m|39.58W/m²|42.54W/m²|2.281h|0.135%|
|0.60m|1.756W/m²|12.76W/m²|14.254h|1.113%|

This demonstrates strong attenuation and delayed ice heat delivery at greater thickness. It does **not** show an improvement in cumulative field melt. Harmonic phase is identifiable only modulo its forcing period; the14.254h entry must not be misreported as an advance of9.746h. The continuum comparison includes finite mesh/time error and the hourly-held forcing.

Over the eight days, final-minus-initial storage is−0.327MJ/m² at0.18m and−0.850MJ/m² at0.60m. Basal heat exceeds surface heat by these amounts, as required. Endpoint storage, not a new heat source, explains that difference.

## Scientific claim limits and integration safeguards

1. Boundary temperature must genuinely be debris-surface temperature, with an independently justified retrieval or measurement. Air temperature is not interchangeable. The solver's explicit quantity gate is useful but cannot authenticate a declaration's truth.
2. Independent surface temperature can test conduction/storage while bypassing a surface-energy closure. That is a different information setting from predicting surface temperature from meteorology. Compare equivalent inputs before claiming model superiority.
3. Fixed debris thickness and basal0°C assume local geometry/material conventions. Melt-driven debris settling, changing moisture, snow insulation and cold ice are not included.
4. `ensemble()` correctly labels a finite declared-scenario envelope. Taking extrema over finitely many scenarios does **not** bound every intermediate continuous parameter/state. No monotonic-in-heat-capacity or monotonic-in-conductivity claim should be invented to upgrade that envelope: time-varying boundaries and initial storage can reverse simple orderings.
5. Numerical verification establishes that the declared conductive model is solved consistently. It does not identify initial thermal storage, prove field parameter accuracy, or repair missing target/geometry metadata.
6. Harmonic thermal development is established physics, not a new theorem or environmental discovery. New field significance would require a target-independent prediction with independently constrained boundary, initial state and material properties, scored against genuinely new native ice-ablation observations.

Review updated for the added bulk-volumetric-capacity and temperate-basal-state provenance gates. The independent fixture explicitly declares both. All4 independent thermal review tests pass against this reviewed version; no numerical results changed. Earlier review outputs remain versioned.
