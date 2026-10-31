# Decisions on the requested extensions

These additions are not all prerequisites for a useful thermal-channel model. Each adds equations, parameters and verification obligations. This contribution implements lake storage/heat coupling and its executable demonstration. It does not contain placeholders advertised as completed physics.

| Extension | Decision | Evidence required before implementation/claim |
|---|---|---|
| Lake storage and thermal channel coupling | Implemented here at the interface and reduced demonstration level; CFD adapter is next. | Actual coupled fluxes, channel storage and boundary-condition tests before claiming CFD integration. |
| Lake-hazard coupling | Later, if outlet evolution materially changes release. | Demonstrated hydraulic control plus downstream routing, topography, uncertainty and an explicitly evaluated hazard quantity. A lake hydrograph alone is not an inundation or risk map. |
| PINN / physics-informed surrogate | Optional later comparison after trustworthy CFD cases exist. | Compare total training+inference cost and held-out prediction error against the conservative reach model, an interpretable closure and an unconstrained learned baseline. No training campaign is launched here. |
| Hydrofracture | Separate research extension if the target includes drainage through cracks or to the glacier bed. | Fracture geometry, toughness, stress state, viscous/elastic deformation, water-pressure coupling and independent fracture benchmarks. Do not emulate it with a threshold outlet flag. |
| Debris transport | Add only if the selected reach contains mobile sediment that materially affects heat or hydraulics. | Grain-size/concentration inputs, entrainment/deposition closure, sediment mass conservation and a coupled benchmark. A roughness multiplier is not debris transport. |
| Full glacier mechanics | Couple an established ice-flow solver if deformation changes geometry on the prediction horizon. | Compare deformation and incision timescales using justified parameters; constrain ice geometry/rheology and test conservative transfer. Do not rebuild a glacier-scale Stokes solver merely to enlarge the project. |

## PINN experiment to keep available

Keep the finite-volume or CFD model as the reference and conservation owner. If speed or missing heat-transfer structure becomes a bottleneck, test a physics-informed thermal surrogate on entire held-out geometry/forcing families. Identify whether the candidate is a PDE-residual PINN, a neural operator or a learned closure; those are not interchangeable labels.

Inputs must be available at prediction time. Do not include future geometry or reference heat flux as hidden inference inputs. Predict either the needed thermal state/flux or a bounded correction to an interpretable closure; apply water-to-ice heat exchange once in the conservative solver. Report boundary residuals, global energy error, wall-flux error, propagated incision error, OOD failures and total compute. An energy penalty alone does not guarantee conservation or physical validity.

Only train after the CFD generation pipeline, split manifest, baseline and scoring code exist. A sensible initial comparison is: equilibrium physical closure; training-calibrated physical closure; conservative thermal-state model; unconstrained ML; physics-informed ML. Keep fitting data and budgets comparable. Retain the neural model only if it earns a measurable accuracy/cost or transfer advantage.

## Primary sources checked for scope

* [Kingslake et al., Modelling channelized surface drainage of supraglacial lakes](https://www.cambridge.org/core/journals/journal-of-glaciology/article/modelling-channelized-surface-drainage-of-supraglacial-lakes/E5154217B4D0307A242FE37C958F228E): lake storage, critical controls and incision-driven drainage are established. The present demonstration does not claim priority for them. Methods/equations and assumptions were inspected.
* [Hageman et al. 2024](https://doi.org/10.5194/tc-18-3991-2024): a coupled thermo-hydro-mechanical fracture model shows that viscous ice deformation can matter during rapid lake drainage. Abstract and bibliographic record inspected; not a full equation-level implementation audit. This is why a crack-opening shortcut is not being added.
* [Hill et al. 2025, GlaDS-GP](https://gmd.copernicus.org/articles/18/4045/2025/gmd-18-4045-2025.html): hydrology emulation already has substantial prior art and application-dependent tradeoffs. Discussion of emulator comparisons and physical constraints inspected. This is a different, subglacial setting, not direct validation of the proposed surface-channel surrogate.

No novelty claim follows simply from combining these components. The project's candidate novelty still needs a supported thermal/hydraulic mechanism and improved predictions on independent cases.
