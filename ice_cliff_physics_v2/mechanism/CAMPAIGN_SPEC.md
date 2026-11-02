# Pre-outcome mechanism campaign specification

This specification precedes the new coupled mechanism campaign, not the independent analytic unit tests. Inputs below are engineering scenarios, not fitted field estimates. Record the exact parent input/configuration and code hash before running.

## Primary question

Does finite wet-debris heat storage, coupled to evolving ice geometry, change sustained post-control-clearing toe contact or matched-footprint face melt by a consequential amount relative to the same passive boundary with zero thermal storage? The question is not whether a warm initial rock can supply extra heat.

## Screen first

For each configured wet patch calculate tau=C/(K+b), NTU=G/F, Gamma=(K+d)/(K+b) and Gamma C DeltaT/(rho_i L A). Record the declared DeltaT bound and ice footprint. If the maximum energetic displacement is below numerical/model uncertainty and no independently demonstrated threshold amplification exists, stop after the screen rather than launch a large campaign. A boundary flip inside the numerical bracket is unresolved, not a discovery.

## Frozen material/hydraulic tuples

Use all 12 tuples in campaign_cases.json: thickness d in {0.02,0.10} m; conductivity k in {0.5,2.0} W/(m K); supplied inflow Q in {0.005,0.02,0.07} m3/s. Volumetric heat capacity is fixed at1.5e6 J/(m3 K), so per-area capacity is Cv*d and conductance k/d. These covary physically; do not independently choose a huge capacity and tiny resistance simply to make a dramatic result.

Use one parent geometry fixed before scores. Prefer its existing admissible development geometry. Set explicit lower-face cover fraction0.05, with upstream hydraulic control and attached coating identified as different objects. Do not flood a3m face merely to activate upper-face debris. Keep all cases where the wet kernel stays inactive; label them unexcited controls. If hydraulic admissibility rejects a tuple, retain the refusal with reason and no class score.

Primary transfer closure: standard DB0.4 with its declared domain checks; no test-fitted multiplier. A second closure may be a fully reported sensitivity, not a selected winner. Run at most24 paired cases (12 per closure) before examining whether more compute is justified.

## Paired runs and inputs

Each tuple compares finite C with the C→0 equilibrium boundary, keeping K, geometry, water route, cover mask and model family identical. Start both with zero extensive debris energy, zero surface cold debt and temperate fresh ice. Same prescribed Q and water superheat0.5K. Net dry-face atmospheric flux is180W/m2 for the same eight daytime hours and0otherwise; floor/bypass atmospheric flux0. There is no artificial warm initial reservoir. The upstream control is removed at24h; attached face debris is not removed. Horizon72h (48h after removal). All forcing and intervention boundaries are hit exactly.

The C→0 run is a mechanism ablation, not a claim that professional debris models ignore storage. Also retain a geometry-state baseline and an infinite-water-temperature-reservoir limiting check where applicable. A model that already carries these physical states can reproduce the effect and must not be excluded from comparisons.

## Mandatory falsifiers

1. Fixed geometry and fixed conductances: the supplied analytic identity must hold, including after two equal-total forcing orders and a common relaxation tail. No lasting thermal-only history effect is allowed after terminal storage matches.
2. No attached debris: capacity changes must have no effect.
3. Zero air and zero initial storage, fixed water forcing: the finite-storage response cannot manufacture more local cumulative phase energy than the complete available budget.
4. Same ice geometry at restart: distinguish remaining thermal state from pre-existing substrate change. Any diagnostic thermal reset must book its energy exchange explicitly.
5. If a total deposited-heat replay is used to distinguish timing/allocation from energy quantity, label it a non-predictive diagnostic and preserve the prescribed heat packets exactly.

## Outputs and acceptance

Report time series of toe flow/fraction, wet bare/debris areas, face/floor/bypass melt on identical footprints, outgoing water heat, debris enthalpy, ice energy and all conservation residuals. Distinguish maintained, lost, genuine reconnection and failed reconnection. Shared pre-removal trajectories are one history.

A candidate consequential reduced-model result needs >=0.10 absolute post-clearing toe-share change OR >=10% matched-footprint face-melt difference, persisting in a predeclared late6h window after the thermal-state mismatch is <1% of the phase-energy difference. It must occur in at least two adjacent predeclared tuples, not one finely tuned case. Refinement error must be <20% of the claimed difference; halve coupling steps twice on the strongest and weakest qualifying cases, with failures retained.

These are explicit engineering decision thresholds, not statistical significance. No confidence interval is inferred from12 dependent parameter tuples. A pass supports a mechanism within the declared reduced closure only. Independent hydraulics/heat/geometry comparison, field input identification and prior-art checking are still required before environmental novelty or predictive-superiority claims.

A fail is useful: classify the change as implementation robustness rather than manufacture a headline by selecting warmer initial debris, larger heat capacity or a different observation threshold.
