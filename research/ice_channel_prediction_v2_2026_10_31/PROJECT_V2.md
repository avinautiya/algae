# A predictive CFD project for thermal incision in glacier channels

Revision 2, 31 October 2026. This specification supersedes the emphasis and implementation order of the earlier ice-step proposal. The earlier mathematical checks remain valid within their stated assumptions. No new full CFD simulation or field forecast was produced in this revision.

## Decision

Build a **predictive model of thermal incision in open glacier meltwater channels**, with a hierarchy connecting resolved CFD to a conservative, faster reach model. Investigate developed ice-step stability as the physical explanation for persistent melt concentration. Do not make cold-ice memory, a unique selected wavelength, or a positive improvement over every baseline prerequisites for the model to exist.

Working title: **Predictive Modeling of Thermal Incision and Hydraulic Instability in Glacier Meltwater Channels**.

Core question: when evolving ice steps and hydraulic jumps invalidate a local equilibrium heat-transfer approximation, can a conservative model informed by resolved water–ice calculations predict channel lowering, heat export and retreat more accurately on unseen geometries and forcing histories?

The project has a software deliverable even if the proposed mechanism is weak: a verified model and an honest validity map. A substantial research contribution additionally requires a new supported result or a useful improvement. A functioning program alone is not scientific novelty.

## 1. What changed after the second literature review

**Prior art is broader than a first reading suggested.** Jarosch & Gudmundsson (2012) already couple channel evolution to thermal erosion and ice mechanics, although their water-flow treatment and melt distribution are parameterized. Explicit turbulent melt morphodynamics also exist. A 2025 JFM study extends this work to melting timescales and parameter scaling. Thus, adding Navier–Stokes or an evolving ice boundary is not a sufficient originality claim. See REFERENCES_AND_DATA.md for primary links and reading depth.

**Prediction is now the primary deliverable.** The earlier proposal could have ended with a stability diagram without a usable forecast interface. This revision requires a forward model driven only by initial states and boundary forcing available at the forecast origin.

**Temperate ice is the base case.** Cold-ice response is a nested physical extension, retained only where relevant. That avoids making environmental usefulness depend on a cold-gradient mechanism that disappears at the melting point.

**CFD and the operational model have different jobs.** Resolved calculations supply heat-transfer evidence and test assumptions. The fast model performs repeated reach predictions. Neither is automatically ground truth.

**The strongest novelty target is a connected result:** a thermal explanation of step persistence or merging, a heat-transfer representation consistent with that explanation, and improved held-out predictions of incision/heat export. Do not publish an elaborate method without testing the final quantities.

## 2. Final prediction contract

The initial scoped product accepts the inputs below and returns physical predictions in SI units. Its short-horizon development cases should span hours to a day where ice creep, debris and roof closure are demonstrably secondary; this is a scope choice, not an assertion that those processes are always negligible.

| Inputs available at forecast origin | Required outputs |
|---|---|
| Initial longitudinal bed and cross-sectional geometry, including uncertainty and spatial resolution | Bed lowering b(x,t)−b(x,0), and cross-sectional change only if the model actually evolves width |
| Initial water depth/velocity and temperature state, or a documented spin-up | Water stage, velocity and outlet temperature/enthalpy time series |
| Either prescribed inlet discharge plus compatible outlet condition, OR lake stage/storage coupling with the appropriate hydraulic boundary conditions | Outlet discharge in transient runs; avoid claiming independent discharge prediction when it is imposed or nearly fixed by conservation |
| Inlet temperature history and its observation definition: point sensor versus bulk/flux-weighted temperature | Integrated melted ice volume and a closed mass/energy ledger |
| Ice temperature profile or temperate-ice declaration | Locations and rates of concentrated erosion; step retreat/migration if features are resolved |
| Atmospheric heat terms where important, friction/roughness inputs and their provenance | Uncertainty/sensitivity envelope, resolution diagnostics and explicit out-of-domain flags |

The initial product is a **reach model**, not a regional climate model, glacier mass-balance model, or GLOF warning service. A useful integration interface can later pass corrected channel geometry, conductance and heat export to a drainage network. Regional usefulness must be measured after coupling; it cannot be inferred from a local animation.

## 3. Precise novelty boundary

Established components: incompressible flow equations; heat advection/diffusion; Stefan melting; moving free boundaries; hydraulic jumps; spectral stability; reduced-order models; thermal memory in general.

Unresolved candidate contribution: a quantitatively tested stability/heat-transfer relation for finite-amplitude, free-surface ice-channel steps, used to predict incision on finite open reaches without prescribing the final step spacing or refitting each case.

Three gates must all be addressed before a strong headline:

1. **Distinct physics:** resolved calculations reveal a thermal–hydraulic effect not already captured by the closest thermal model, beyond discretization and closure uncertainty.
2. **Mathematical explanation:** a derived criterion or reduced operator explains that effect and predicts independent transitions or rates. Merely using an existing spectral method does not pass this gate.
3. **Prediction:** the resulting reach model improves important held-out quantities against both established and fairly calibrated empirical baselines, at useful cost.

If only the third gate succeeds, a useful engineering/modeling paper may still result. If only the second succeeds, a theoretical paper may result. Do not call either a proved new physical law. No search can guarantee absence of all prior work; the specific mechanism remains unconfirmed until the calculations and final citation audit.

## 4. CFD hierarchy

**Level 0 — exact and manufactured checks.** Stefan conduction, advection of a thermal pulse, hydrostatic rest, conjugate-depth jump relations, phase-change mass balance, and the prior packet's thermal admittance/jump derivatives. These diagnose implementation, not environmental skill.

**Level 1 — resolved laminar water–ice CFD.** Solve the incompressible Navier–Stokes and temperature equations in a vertical plane with a real free surface, heat conduction in ice where retained, and a Stefan-driven moving boundary. Use it for a small set of controlled, low-Reynolds-number geometries. Resolve the thermal boundary layer. Include surface tension where its local scale demands it. Two-dimensional calculations cannot substantiate three-dimensional turbulence claims.

**Level 2 — turbulent checks.** Use three-dimensional URANS with declared turbulence/thermal closures for screening, and selected LES where resources and resolution permit. Assess mean wall heat flux, free-surface statistics and jump/reattachment location, not just the mean water profile. RANS agreement with itself is not validation; LES is not DNS, and two-dimensional unsteady flow is not resolved three-dimensional turbulence. Compare turbulence/thermal closures and discretizations on a small independent subset.

**Level 3 — conservative reach predictor.** A depth-averaged flow/enthalpy model updates channel geometry through the same computed interface energy exchange. Use the simplest heat-transfer representation that survives Level 1/2 tests. Add a transported thermal state or history dependence only after matched-state tests show a material need. Reuse the existing network's infrastructure after a code audit, not its labels as evidence of physical capability.

**Implementation preference:** use an established CFD code rather than writing a full multiphase turbulent solver from scratch. OpenFOAM supplies relevant flow and heat-transfer building blocks; it is not verified here as a turnkey free-surface ice-melting solver. Standard gas–liquid phase-change examples are not solid–liquid Stefan benchmarks. The next milestone must identify a pinned, working solver path and demonstrate the required coupling. Basilisk has relevant research modules, but their composition, verification and maintenance also require audit. See CFD_AND_VALIDATION.md.

## 5. How to make the expensive calculations count

Start with prescribed geometries to isolate the heat-transfer mechanism. Use a smooth reach, a single step/jump and a short train with controlled spacing; keep input water/heat budgets comparable. Separate geometry effects from changes in bulk discharge or supplied heat. Avoid tuning inlet turbulence individually to match each target.

Then allow melting and compare predicted geometry against the frozen-geometry heat-flux tendency. Track whether the flow follows the evolving bed quasi-steadily or requires unsteady coupling. The mean water transit time, jump fluctuation time, thermal relaxation time, morphology time and forcing time must be estimated separately.

The moving-boundary calculation must preserve physical latent heat. If morphology is advanced using a large outer time step while the flow is sampled on a smaller step, establish convergence in that outer step and sampling window. Changing latent heat to make melting faster changes a dimensionless parameter; it is not merely a faster clock.

Do not assume infinite periodic steps model an open channel. The finite reach receives heat at its inlet and exports remaining heat at its outlet. A periodic thermal calculation must state its heater/forcing. The jump stability family is an analysis tool; the final forward model must evolve from initial geometry without being assigned the desired endpoint.

## 6. Evidence before elaborate heat-transfer models

MATHEMATICAL_UPGRADES.md contains a constructed pair of temperature profiles with identical bulk hydraulic/thermal descriptors and different instantaneous ice heat flux. This illustrates why an instantaneous bulk formula is not universally exact. It does not establish that such differences persist in glacier conditions.

To establish practical relevance, generate matched histories in resolved CFD. At a comparison section, match current depth, flux, geometry and bulk temperature within a predeclared tolerance, but vary physically attainable upstream thermal development. Compare wall heat flux and subsequent incision, with error bars from mesh/temporal statistics. Repeat after several residence/relaxation times. If differences disappear or are below uncertainty, use a local closure and abandon the memory claim.

If a retained thermal state is needed, derive it from energy moments or projection of the resolved transport equation. Enforce exact exchange cancellation between water and ice, and verify limiting cases. A neural network is optional and must beat an interpretable state model on entire held-out cases. It is not the source of novelty.

## 7. Stability as explanation and prediction

Use the earlier moving-jump and thermal-response foundation to investigate the stability of developed profiles. Include perturbations spanning several steps, not only one period. Compare the maximum-admissible-spacing prescription with a dynamically stable band; do not require a unique selected wavelength.

Test predicted growth/decay and migration in independently discretized time evolution. Then assess whether the predicted stability regime explains finite-reach melt patterns and their persistence. An eigenvalue change with no observable consequence is not sufficient evidence for environmental impact.

Complete the bulk heat-transfer problem before adding full cold-ice memory, finite-amplitude spectral continuation, lateral erosion and lake coupling simultaneously. These are successive levels of evidence, not mandatory ingredients in the first solver.

## 8. Benchmarks and fair comparisons

Compare:

* conventional equilibrium heat-transfer correlations, with their stated applicability;
* a training-calibrated empirical heat-transfer law using the same permitted inputs;
* the present thermal reach model with fixed/local closure;
* the proposed closure/coupled model;
* the resolved reference on a manageable subset.

For step-spacing questions, include the closest published thermal profile model separately. Do not force models designed for other geometries to make predictions outside their scope and count that as a win.

Hold out entire geometry/forcing families, not snapshots from the same simulation. Split development and test parameter spaces before production. Keep forcing, initial-state information and parameter fitting budgets comparable. A published empirical fit using an entire field event is a contextual reference, not an independent held-out baseline for that same event; refit using training data for the fair test.

Primary endpoint: time-dependent bed lowering over the forecast reach, reported as MAE/RMSE in metres and relative to a fixed predeclared physical scale. Co-primary endpoint: outlet sensible-heat export integrated over time, in joules. Secondary endpoints: total ice volume lost, migration distance, concentration of erosion and runtime at matched accuracy. Avoid dividing by nearly zero melt. Report biased predictions even if an aggregate score improves.

Pre-register the smallest scientifically meaningful improvement after quantifying reference and observational uncertainty, before looking at test outcomes. Do not choose a universal 10% or 20% threshold without that justification. Confidence intervals are by independent forcing/geometry case; correlated snapshots do not increase the sample size.

## 9. Environmental impact and its causal limit

Better incision predictions could improve how drainage models evolve channel geometry and transmit warm water. Whether this materially changes lake release or network routing requires a separate coupled experiment.

For shallow-water dynamics the local characteristic speeds are u±sqrt(g h). An entirely supercritical reach does not transmit small hydraulic disturbances upstream. A downstream melt hotspot therefore cannot automatically change an upstream lake's release. Track actual hydraulic connectivity and control migration before claiming hazard impact. This characteristic argument does not exclude later large-amplitude or geometry-mediated changes.

The temperature and volume of supplied water also impose an energy ceiling: reshaping a channel redistributes supplied heat; it does not create unlimited additional melting. Report retained/exported heat alongside melt. If the net effect is mainly spatial redistribution, describe that accurately and explain whether the location matters to a hydraulic control.

For Himalayan use, begin with clear-water, ice-bounded reaches meeting the model assumptions. Debris transport, slush, pressurization, hydrofracture, roof closure and major ice deformation need separate processes or exclusion flags. Satellite imagery may eventually constrain geometry/forcing at adequate resolution; it cannot supply unresolved boundary-layer state or prove centimetre-scale incision accuracy.

## 10. Data and attainable validation

The public ETH collection for Ogier et al. provides a potential field test. Its metadata lists a roughly 322 MB archive and README; downloads were rate-limited in this session. The article establishes substantial boundary/measurement uncertainty. REFERENCES_AND_DATA.md records the exact observation issues. Do not advertise a completed data audit or calibrated model yet.

Use published laboratory observations for the correct regime and independent numerical checks for mechanisms. Field comparisons require an observation operator and initial/boundary uncertainty, not direct comparison between bulk model temperature and a near-bed sensor. Geometry at later times cannot be used as an input when it is the prediction target.

The first deliverable can be a working computational predictor without collecting new observations. Its demonstrated validity will be the regime actually tested. A claim of improved real-world forecast skill requires independent observations sufficient to evaluate that claim.

## 11. Completion criteria

A successful final software package must accept a documented case, run from initial conditions, produce the promised time series/maps, conserve mass and energy within a verified numerical tolerance, reproduce a held-out benchmark, and refuse or flag unsupported regimes. It must expose parameters and uncertainty rather than burying them in source code.

A successful scientific paper must additionally establish something new and significant: a mechanism, criterion, transferable closure or predictive improvement. Neither Nature publication nor STS selection can be assured in advance. This revised scope gives a concrete way to earn a strong result rather than increasing complexity until the project sounds impressive.
