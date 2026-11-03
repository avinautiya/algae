# Kennicott temperature-assisted physical melt transfer protocol

## Fixed pre-target design

Primary `physical_layered_local` is a fixed-geometry heterogeneous-debris/cold-ice conduction model with interface complementarity. It infers two positive effective diffusivities from local temperature profiles, then calibrates ONE global volumetric debris capacity using the entire source-bound2020TRAIN mass-loss interval. No2023 melt target may be accessed before every2023 prediction is published. This is a conditional temperature-assisted hindcast, not an advance weather forecast, new turbulence law or guaranteed improvement.

All ten methods are retained:

1. `physical_layered_local` (primary).
2. `physical_homogeneous_local` (same temperature data, one positive diffusivity).
3. `physical_fixed2020_layered` (all sites use2020 layered coefficients).
4. `physical_weather_boundary_layered` (entire-period air-derived boundary, same local layered coefficients; temperature calibration still uses initial post-installation data).
5. `homogeneous_empirical_scale` (same homogeneous thermal method at fixed referenceC1.5MJ/m³K, then one2020 output gain; scaled output is explicitly not an energy-conserving physical melt volume).
6. `air_pdd`.
7. `buried_temp_exposure`.
8. `conventional_gradient`, the covered-temperature degree exposure divided by stake thickness minus measured upper boundary depth.
9. `depth_corrected_air_pdd`, the air dose divided by the same thickness.
10. `constant_rate`.

The five deterministic indices and empirical homogeneous comparator each receive the same single global2020TRAIN amplitude. No test-site melt calibration, local capacity fitting or model selection is allowed. A physical baseline calibration refusal is retained as unsupported, never silently dropped or counted as a numerical win; the empirical homogeneous comparator ensures the closest standard thermal method remains numerically comparable.

## Source correction and material inference

A source summary describes2020b as having two thermistors. The released file has a third22cm column identically zero. It is treated as an appended boundary, not an independent measurement. The primary instead uses genuine measured2020a temperatures to infer diffusivity, transferring them to2020b's forcing/stake column. Identical source coordinates, similar grain description and concurrent records support an explicit proximity/material assumption; they do not prove the same pit. The warm25cm2020a measurement is not forced to0C or treated as a confirmed ice interface.

At each material profile, use the deepest THREE genuine thermistors, six-hour consistent-P1 integrated balances, and the first14days after a five-day installation-reset exclusion. The exclusion applies only to parameter fitting and upper-boundary reconstruction; it does not truncate the observed whole-campaign melt target. Profiles are2020a11/18/25cm;2023a12/17/22cm;2023b16.5/24/31.5cm;2023c9.5/17/24.5cm. Parameters are static per site, not causal rolling adaptation. Source-derived temperature uncertainty0.2C supplies a necessary outer residual diagnostic, not confidence intervals.

NNLS positive-material and singular-ratio0.001 gates pass for all four profiles. All56 integrated balances at each profile satisfy the conditional sensor-error outer radius. This does not establish a single coherent sensor-error explanation, pure conduction or independently measured conductivity.

The physical homogeneous comparator constrains alpha1=alpha2 in the same integrated equations, without another search parameter. The fixed2020 comparator keeps the2020 coefficient pair and uses the same site-specific layer-interface geometry as the primary. This isolates coefficient adaptation while retaining temperature/geometry information available to every thermal method.

## Forward geometry and boundary

The computational top is the first genuine measured buried-temperature plane, not radiometric surface temperature. The bottom is the independently documented stake thickness:2020P_AWS_1=0.21m,P11=0.24m,P12=0.29m,P13=0.255m. The top planes are0.05,0.02,0.09,0.02m. Thus the computational thicknesses are0.16,0.22,0.20,0.235m. Profile thicknesses are kept as sensitivity metadata, not substituted silently for stakes.

Inferred alpha1 extends upward to the top plane; alpha2 extends below the middle thermistor to the stake ice boundary. These are explicit effective-material assumptions. For P12 the31.5cm profile sensor informs material inference only; it is never inserted inside the29cm forward column, nor imposed as a boundary there.

The boundary helper fits an affine buried-temperature/AIR2 relation on the same fixed14post-reset days; its blocked7→7day temperature residual diagnoses reconstruction error without melt targets. The initial missing36/48h and five-day disturbed period are supplied by that proxy, not zero-filled or removed from the seasonal target. All later measured top temperatures are used. Hourly boundary samples are interpolated linearly to300s backward-Euler endpoints. Native AIR2 is the primary shared weather source; AIR3 and±2blocked-RMSE reconstruction shifts are predeclared sensitivity scenes, not probability intervals or outcome-selected corrections.

The nominal visit-time operator is noon ADT(UTC−8), a declared assumption because the stake source supplies dates without hours. It is not a measured timestamp. Clock-window/initial-state/geometry sensitivity must be reported as conditional scenarios; the nominal point cannot certify coverage of all possible visit times. Some late-day clock scenes may exceed available weather and must refuse rather than extrapolate.

## Physical solver and calibration

The model uses16 debris cells,20 ice cells over2m, and300s steps. Debris conductivity is k_j=C*alpha_j. Internal conductances are harmonic half-cell resistances. Ice k2.1W/m/K and volumetric capacity917*2100J/m³/K are explicit model assumptions, distinct from the source's900kg/m³ density used to report field mass. The native mass-change target is already kg/m²; it is not density-converted again.

At each step, first solve a melting candidate at interfaceTb0C. If net interface heat is nonnegative, that exact heat becomes latent melt. Otherwise solve the coupled debris/ice sensible system with a harmonic interface resistance and zero melt, yielding Tb<=0. No refreezing water is invented. Initial debris is a linear top-to0 profile and ice initially0C; colder-ice scenes are separate. Deep ice is fixed at0C. The solver is fixed-interface, not a moving-front/advection model. Wet-debris phase change, moisture advection and collapse are absent and must remain scope limitations.

The ledger is surface energy = debris/ice storage change + latent energy + deep outgoing energy, using the SAME discrete fluxes as the equations. Posthoc scaling of this physical ledger is forbidden.

Before TRAIN release, probe C at0.8,1.5,3.5MJ/m³K. These are engineering limits, not a measured site uncertainty interval. The sampled2020 response increases for all three physical calibration families, but no continuous monotonicity theorem is asserted. After publication and TRAIN access, require the target inside the unchanged bracket and bisect to relative tolerance1e−5, at most24steps, checking each new response against its bracketing values. Outside-bracket or nonmonotone responses cause explicit refusal; no capacity widening or output-gain rescue of the primary is allowed.

Homogeneous and weather-boundary physical baselines each calibrate their own ONE global2020C using the same limits. Fixed2020 layered shares the primary's TRAIN response/C. Every2023 physical prediction keeps its calibrated C and physical energy ledger unchanged.

## Verification and feature-only outputs

Five core physical tests pass: exact uniform and heterogeneous steady thermal resistance/latent flux; signed cooling into cold ice without imaginary melt; cold initial ice reducing melt with ledger closure; and invalid-material/above-melt-ice refusals. Eight driver tests pass: bounded calibration, endpoint equality, nonmonotone/outside-bracket refusal, no density reconversion, P12 geometry and corrected indices.

An independent Fourier-series continuum step test separates conservation from numerical accuracy. At a sharp two-hour step,16cells/dt300s overpredict cumulative melt by about4.49%; dt150s reduces it to2.25%, with much smaller spatial variation16→32cells. This numerical limit is preserved in `TRANSIENT_VERIFICATION.json`; the excellent discrete ledger is not treated as physical accuracy. A declared refinement scene can bound effects on seasonal predictions. No time step is selected using field errors.

`FEATURE_ONLY_PLAN_FINAL.json` contains source bindings, temperature-only parameter estimates and all ten IDs. `PRETRAIN_CAPACITY_PROBES_FINAL.json` records target-free2020 capacity responses. `NOMINAL_REFERENCE_EXPOSURES_FINAL.json` holds all four physical reference runs atC1.5MJ and five deterministic index doses for all four stakes; the empirical homogeneous exposure is its corresponding physical reference total. None uses mass outcomes. Earlier plan versions are retained as pre-target metadata-design history, with the FINAL files authoritative.

## Release sequence

1. Publish code, protocol, source hashes, feature plan and physical capacity probes before2020TRAIN magnitudes.
2. Release ONLY the source-bound2020TRAIN target; perform the fixed one-parameter calibrations, retaining refusals.
3. Compute all three2023 stake predictions, all ten methods and every predeclared sensitivity scene without reading test magnitudes.
4. Publish final predictions and hashes.
5. Release exactly the three frozen2023 target rows, score same-operator seasonal errors and all predeclared comparisons. The three stakes share a glacier/year; they are not three independent climates. Report any failed superiority result unchanged.

## Final pre-TRAIN seasonal refinement

The interface is now integrated exactly within each cell using its harmonic layer-intersection resistance. This preserves the total series resistance at every grid size; the layer is not assigned merely by cell-center location. All four source-forced sites, primary and homogeneous methods, pass independent time/debris-mesh refinement at reference C1.5MJ: (16cells,300s) to (16cells,150s) to (32cells,150s), keeping20icecells fixed. Across all eight comparisons, the largest seasonal integral change is0.003814%, well below the declared1% engineering gate. This does not bound physical/model/source error or establish ice-mesh convergence. `SEASONAL_REFINEMENT_FINAL.json` contains every run. The300s primary is therefore retained before any TRAIN target is read.

Earlier center-assigned reference artifacts and pre-final plan files remain as pre-target design history. Their P11/P13 values are superseded by the harmonic-intersection FINAL exposures; no outcome influenced this correction.
