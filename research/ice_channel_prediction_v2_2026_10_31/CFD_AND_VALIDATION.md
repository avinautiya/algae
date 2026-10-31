# CFD implementation and predictive evaluation specification

## A. What CFD means here

The intended reference solves spatially resolved water momentum and temperature, with a free surface and a thermally evolving ice boundary. A depth-averaged hydraulic solver is also a numerical fluid model, but it does not resolve the vertical recirculation or thermal boundary layer needed for the key comparison.

The model chain is: analytical checks → resolved laminar benchmark → selected turbulent CFD → conservative reach predictor → held-out prediction tests. The final predictor need not run expensive CFD for every forecast; any approximation it uses must be evaluated at the output level.

## B. Solver capability audit before implementation

Use a mature numerical framework. Preferred initial audit: an available, pinned OpenFOAM release. Document the exact executable/source revision, compiler, mesh treatment and case dictionaries. Inspect:

1. Air–water free-surface transport and flux consistency.
2. Water temperature/enthalpy transport using the same phase-consistent mass flux.
3. Ice conduction with distinct material properties, when cold ice is retained.
4. The Stefan boundary energy balance and phase-change mass transfer.
5. Geometry motion, remeshing and conservative transfer between meshes.
6. Turbulence and turbulent heat diffusivity; near-wall resolution/thermal treatment.

Do not assume that interFoam plus an arbitrary melting source combines correctly. Official OpenFOAM documentation lists free-surface, conjugate heat-transfer and solidification capabilities, but their existence does not verify the required combination. A liquid–vapour Stefan example is not interchangeable with ice melting. Basilisk is a possible alternative only after equivalent checks of the selected research modules. Choosing a solver by name is not a completed architecture.

For the first resolved case, hold ice geometry fixed and prescribe a melting-temperature boundary, solve water flow/temperature and extract heat flux. This isolates the hardest heat-transfer assumption before mesh motion is introduced. Then introduce moving ice and verify against a one-dimensional Stefan solution and a uniform-heat recession case.

## C. Boundary conditions and approximation decisions

* Match hydraulic boundary-condition count to incoming characteristics. Do not prescribe discharge, depth and downstream stage indiscriminately.
* Specify whether inflow temperature is uniform, fully developed or measured. If bulk temperature is supplied without a profile, include a documented entrance region/uncertainty.
* Use a genuine atmospheric free surface when jump shape and entrainment matter. A rigid lid removes relevant physics.
* Retain physical surface tension for the millimetre-scale reference unless a local scale analysis justifies removal. Surface tension at the free surface is different from ice-interface regularization.
* State initial and lower-boundary ice temperatures. Demonstrate independence from artificial bottom depth. Temperate ice is the primary initial model.
* Apply the phase-change normal velocity and mass flux from MATHEMATICAL_UPGRADES.md. Tangential ice material velocity and geometric wall velocity are different quantities.
* Include buoyancy only after evaluating its importance relative to shear; near freshwater's density maximum a constant expansion coefficient may be inadequate. Do not add a density model without testing the regime.
* Introduce atmospheric radiation/turbulent fluxes only with stated inputs, units, area and energy bookkeeping. In some short warm-water cases these may be negligible, but that must be a quantitative bound.

## D. Verification sequence

| Case | What it checks | Required evidence |
|---|---|---|
| Hydrostatic rest with varying bottom | Pressure/source balance | Velocity and free-surface error versus mesh |
| Laminar flow with known velocity profile | Momentum/viscosity and wall conditions | Profile error and observed convergence |
| Thermal advection/conduction on a fixed domain | Enthalpy transport and diffusivity | Pulse/steady solution, conservation, refinement |
| One-dimensional melting | Latent heat and interface motion | Stefan front and temperature error versus mesh/time step |
| Phase-change mass ledger | Density ratio and moving mesh | Ice mass loss equals added water mass; no double counting |
| Fixed-bed hydraulic jump | Free surface, momentum loss, statistics | Conjugate depths plus a suitable published flow benchmark |
| Warm flow over fixed ice-shaped geometry | Spatial heat transfer | Wall-flux and outlet-enthalpy balance, thermal-layer refinement |
| Evolving single step | Coupled geometry/heat feedback | Joint space/time/morphology-step convergence |

A correct conjugate-depth ratio alone does not validate near-wall turbulence or heat transfer. Run at least three systematically refined meshes on representative cases, reduce the time step separately, and vary morphology-update interval. Use Richardson/GCI-style estimates only when a coherent asymptotic trend supports them; otherwise report a resolution envelope. Time-averaged turbulent statistics need a sampling uncertainty based on correlation time, not the number of saved frames.

Numerical acceptance should be tied to the prediction tolerances and uncertainty budget. A default scaled conservation target can guide development, but rounding-level agreement is neither expected nor evidence of correct physics in the full multiphase case. Show residual definitions and denominators, particularly near zero flow.

## E. Model reduction and coupling

Start the reach predictor with the existing best thermal closure. Replace it only when resolved evidence identifies a deficiency. A candidate extension may transport a thermal-profile moment, use a spatial development length, or use geometry-conditioned heat transfer. Choose among them using development cases; compare with a flexible local empirical baseline. Do not require all three.

The heat-transfer closure must return both the exchanged power and its spatial distribution with a clear area convention. Use the same power in water cooling and interface melting. Evolve the ice shape and then recompute hydraulic geometry. Preserve input energy and remapping budgets. If only longitudinal geometry is evolved, width is an input and lateral retreat is not an output.

Interface controls should expose: physical parameters; applicability range; uncertainty representation; required state/history; and failure flags. Invalid cases must not silently receive clipped temperatures, artificial depth floors or arbitrary heat multipliers that alter the physics.

## F. Predictive evaluation

1. **Freeze development/test families.** Separate geometry shapes, forcing waveforms and thermal regimes. Keep nearby variations of the same trajectory in the same split.
2. **Fit baselines fairly.** Give local empirical models access to the same permitted training information. Include the incumbent physical model and report all tests, not the best case.
3. **Predict forward.** Inputs stop at the forecast origin except for explicitly supplied exogenous forcing. Later observed geometry is a target, not an input.
4. **Score end quantities.** Bed-lowering error in metres; integrated outlet heat in joules; ice loss in m³; feature displacement in metres; matched-accuracy runtime and memory.
5. **Test causality and attribution.** Hold geometry fixed, suppress shape feedback, and toggle only supported thermal terms. Keep energy supply matched. These ablations explain effects; they are not automatically the strongest baselines.
6. **Quantify uncertainty honestly.** Distinguish discretization, turbulence/closure spread, input uncertainty and sampling uncertainty. An ensemble envelope is not a calibrated 95% predictive interval.
7. **State scope.** Numerical-reference prediction, laboratory observation prediction and field prediction are three different accomplishments.

For observational tests construct an observation operator. Near-bed temperature, depth-mean temperature and flow-weighted temperature differ. Instrument accuracy stated as ±a is not automatically a Gaussian standard deviation. Derived Nusselt numbers share inputs with incision/temperature and cannot be counted as independent additional observations.

## G. Environmental application gate

Initially support clear fresh water, an open atmospheric channel and an ice boundary without significant sediment transport. Document geometry scale, Reynolds/Froude/Stefan definitions and horizon. Test temperate ice, then cold ice if pertinent. For highly aerated jumps, three-dimensional turbulence, lateral widening, debris, slush, roof closure or fractures, either supply a demonstrated model or flag the case outside scope.

Only couple to a lake-storage/drainage model after reach predictions are useful. Supply the lake hypsometry and inflow, predict outlet discharge with a physically compatible boundary and compare discharge-volume timing. If the effect is hydraulically isolated from the lake, do not manufacture a hazard benefit. The project can remain valuable as an incision/heat-export model without becoming a flood forecast.

## H. Resource plan

Do not start full glacier-scale DNS. Profile a small valid case: cells, physical timestep, simulated seconds per wall hour, peak memory, checkpoint size and disk growth. Estimate later cost from that measurement. Two-dimensional laminar tests are for verification; expensive three-dimensional cases are selected to resolve specific uncertainties.

Use a stable execution environment, per-case scratch, atomic checkpoints and resumable case metadata. Do not rely on an active chat to prevent worker reclamation. Pin physical parameters and never silently lower latent heat, Reynolds number or thermal resolution to fit a transient container. If resources only support the laminar scope, report that limitation instead of calling it environmental turbulent CFD.
