# Independent CFD verification and response-analysis tools

Python 3.9+ standard library. This package does not import, modify or launch the active OpenFOAM/reach solvers. It consumes exported data and supplies independently implemented classical references. No CFD installation or network is required.

Run from the repository root:

```sh
python3 -m unittest discover -s ice_cfd_checks/tests -v
python3 -m ice_cfd_checks reference --output /tmp/ice_reference_01
python3 -m ice_cfd_checks score-jump --input /tmp/ice_reference_01/jump_exact_example.json --output /tmp/jump_check_01.json
python3 -m ice_cfd_checks score-rest --input /tmp/ice_reference_01/hydrostatic_exact_example.csv --surface 1 --depth-scale 1 --output /tmp/rest_check_01.json
python3 -m ice_cfd_checks score-dam --input /tmp/ice_reference_01/dam_break_wet_exact_example.csv --left-depth 1 --right-depth 0.1 --time 0.3 --output /tmp/dam_check_01.json
python3 -m ice_cfd_checks score-thermal --input /tmp/ice_reference_01/thermal_exact_example.csv --output /tmp/heat_check_01.json
python3 -m ice_cfd_checks integrate-faces --input /tmp/ice_reference_01/face_flux_example.csv --output /tmp/face_check_01.json
python3 -m ice_cfd_checks convergence --input /tmp/ice_reference_01/convergence_exact_example.json --output /tmp/grid_check_01.json
python3 -m ice_cfd_checks linear-response --input /tmp/ice_reference_01/linear_response_manufactured.json --output /tmp/response_check_01.json
```

All output paths must be new. Generated examples are analytical/manufactured fixtures, **not observed CFD results**. Replace example measurements with correctly exported simulation data to score an actual run. The tool does not assign a universal pass threshold or a novelty verdict.

## What is ready

| Component | Function | Important limit |
|---|---|---|
| Stationary rectangular jump | Sequent depth, discharge, ideal momentum, head loss, mechanical power loss | Horizontal, hydrostatic remote sections, unit profile factors; no prediction of jump position or turbulence |
| Hydrostatic rest | Depth and wet-surface deviations, parasitic speed | Spatial samples alone do not establish a volume budget |
| Wet/dry dam break | Classical shallow-water exact states and finite-volume cell averages | Long-wave, flat-bed, inviscid reference; not an exact short-time gate-opening VOF solution |
| Stefan reference | Front, liquid temperature, interface flux and cumulative energy at Ste=0.001, 0.003, 0.013 and 0.2 | One-phase, equal-density conduction limit; not a density-jump or moving-mesh verification by itself |
| Face-flux integration | Separate inflow/outflow volume and donor-temperature enthalpy | Requires actual numerical liquid flux; air flux and arithmetic temperature are inappropriate substitutes |
| Thermal ledger | Absolute/relative residuals including advective, conductive, external and water-to-ice terms | Terms must be complete and disjoint; a closed ledger does not certify the physics |
| Three-grid diagnostic | Conditional observed order and Richardson estimate | Constant refinement ratio and monotone convergence; no automatic confidence interval |
| Finite-time response | Squared singular gain of C exp(A t) B | External physical operator and normalization required; small matrices <=32 only |

## Input schemas

* Jump JSON: `width_m`, `upstream: {depth_m, discharge_m3_s}`, `downstream: {depth_m, discharge_m3_s}`. Optional `gravity_m_s2`, `rho_kg_m3`.
* Rest CSV: `bed_m,depth_m,velocity_m_s`. Optional location columns are retained only in the source file; each sample receives equal weight.
* Dam-break CSV: `left_m,right_m,depth_mean_m,hu_mean_m2_s`, sorted nonoverlapping cells, with x=0 at the initial gate. Use cross-sectionally depth-integrated quantities for a VOF comparison. Cell means are compared to exact cell means, including cells cut by the shock.
* Face CSV: `area_m2,liquid_fraction,normal_velocity_m_s,liquid_volume_flux_m3_s,temperature_excess_K`. All normals point outward. Export the liquid flux from the actual numerical transport discretization and the consistent donor/face temperature. Reconstruction `alpha*U_n*area` is reported separately to reveal differences; it is never substituted automatically. Opposite-direction flows are separated before temperature averaging.
* Thermal CSV: `time_s,water_sensible_J,advective_in_J,advective_out_J,conductive_in_J,conductive_out_J,external_J,to_ice_J`. Flux columns contain time-integrated energy since a common origin. The score uses changes from the first supplied row, allowing a consistent restart origin. Store heat into water as positive input, and heat leaving water as positive output. `external_J` may be signed. Avoid overlapping conductive wall patches with `to_ice_J`; the latter is total water-to-ice heat, including conduction into cold ice when relevant, not merely latent heat.
* Convergence JSON: three objects `[{spacing_m, value}, ...]`, ordered coarse to fine. Specify a physical grid spacing, not total cell count. The compact utility deliberately rejects unequal refinement ratios.
* Linear-response JSON: `A_per_s`, `times_s`, optional `B` and `C`, plus nonempty `state_normalization`, `input_normalization`, `output_normalization`, `operator_provenance` strings. A acts on normalized perturbation coordinates. Defaults B=C=identity are appropriate only if those coordinates and their Euclidean norm have a physical justification.

Each score hashes its input. Missing or nonfinite data are rejected; earlier outputs are preserved.

## Relation to Claude's active milestone

The active task owns the CFD build, case setup, hydraulic checks, thermal solver and melting implementation. This package is offline support. There are no OpenFOAM dictionary edits, runtime hooks, environment changes or job controls here. Current runtime errors cannot be diagnosed from this package; read the real solver log in the active environment.

The reference equations are established science. `SCIENTIFIC_STUDY.md` proposes a separate, falsifiable physical study after the existing milestone passes. It does not authorize a new production campaign now.

Primary-source links, reading scope and numerical assumptions are in `SOURCES.md` and `HANDOFF.md`.
