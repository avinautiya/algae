# Research decision and evidence plan

31 October 2026. Recommendation: a conditional pivot to **Thermal Stability and Melt Localization in Supraglacial Ice Channels**. More descriptive working title: **Dynamic Stability of Thermal Cyclic Steps in Supraglacial Meltwater Channels**. These titles describe the research question; neither claims a discovery already made.

## The contribution worth pursuing

Predict which water-carved ice steps persist, merge or disappear, and where the resulting flow concentrates melting. Replace an assumed step spacing with a tested dynamical stability criterion. Determine whether heat advection and cold-ice response change that criterion in a physically important regime.

This is substantially more fundamental than another forcing-averaging algorithm. Its outputs would be a stability boundary, a range of persistent spacings, migration/growth rates and spatial melt predictions, with a defined validity range. Its success is not established. A journal or STS outcome cannot be guaranteed, and there is no honest complexity threshold that guarantees originality.

## What the literature actually supports

[Wu & Izumi (2023)](https://doi.org/10.3390/geosciences13050128) construct thermal ice-step profiles and adopt a maximum-wavelength choice. [Wu's thesis](https://doi.org/10.14943/doctoral.k15621), printed pp. 33–36, makes the associated stability assumption explicit. That provides a specific testable assumption, rather than a claim that nobody models ice steps.

[Yokokawa et al. (2016)](https://doi.org/10.1002/2015JF003736) analyze flat-interface instability and reproduce ice steps experimentally. That is distinct from stability of an already developed, jump-bounded step train.

Crucially, [Balmforth & Vakil (2012)](https://doi.org/10.1017/jfm.2011.555) already analyze nonlinear cyclic-step stability and secondary instabilities for erodible beds. Consequently, neither step coarsening nor a traveling-wave eigenvalue calculation is a new general idea. [Johnson et al.](https://arxiv.org/abs/1803.03484) already supply a spectral framework for discontinuous hydraulic wave trains. The possible contribution lies in a demonstrably different thermal ice result, not renaming these methods.

No directly matching coupled thermal, finite-amplitude ice-step stability study was located in the targeted searches recorded in SOURCES.md. That is provisional evidence of a gap, not proof of priority. A source could exist under another terminology or outside the indexed material. Do not write “first” before a final forward/backward citation audit and mentor review.

## Directions screened and why this one is strongest

| Direction | What could be useful | Main originality/feasibility obstacle | Decision |
|---|---|---|---|
| Another event-aware averaging or outlet-screen method | Cheap network predictions | The user's frozen experiment already lost to simple forcing blocks | Stop expanding it |
| Ice scallops and turbulent melt patterns | Visually compelling morphology and heat transfer | Extensive existing laboratory and numerical work, including 2024 turbulent-stream morphodynamics | Not a novelty claim by itself |
| Oscillatory forcing changes melting | Phase-dependent heat transfer | Direct oscillatory melting studies exist in 2025 and 2026 | Reject generic version |
| Add cold-ice conduction or refreezing to drainage | Wider physical scope | Established physics; addition alone is insufficient | Use only when an effect is identifiable |
| A new neural surrogate/PINN for the same simulator | Speed or inversion | Does not repair missing physics automatically; prior project lacks a demonstrated surrogate advantage | No network in the core project |
| Source-aware transient conduit heat transfer | Better heat partition | Source-dependent coefficients and time-scale arguments already exist in Sommers & Rajaram (2020) | Backup only if a real separation-of-scales failure is established |
| Heat-budget bounds on critical outlet discharge | Provable conditional safety statements | Prior packet's identities use standard hydraulic/control mathematics; full environmental correction remains open | Retain as a separate theoretical tool, not inflate into a new law |
| Stability and thermal selection of developed ice steps | Predict persistent geometry and melt concentration without imposed spacing | Needs coupled eigenproblem and independent checks; sediment analog is strong prior art | Best scoped candidate from this review |

This is a focused comparison across several mechanisms, not an exhaustive search of all possible ice research.

## The scientific claim ladder

1. **Already completed here:** conservation-based formulation; flat-slab thermal response and passivity derivation; a diagnostic coupled dispersion matrix; moving-jump perturbation conditions; verification code for selected subproblems.
2. **First publishable candidate result:** establish or refute the maximum-spacing stability assumption for the thermal ice model over a stated parameter domain. A converged counterexample would matter; one noisy case would not.
3. **Stronger physics result:** identify a mechanism by which thermal coupling changes persistent spacing or merging, show it survives reasonable closure changes, and explain it using an interpretable dimensionless criterion.
4. **Environmental contribution:** demonstrate a change in melt distribution, retreat or hydraulic resistance in finite, physically consistent open reaches. Include uncertainty from turbulence and boundary conditions.
5. **Not presently supported:** improved whole-glacier mass balance, Himalayan forecasts, operational agency adoption, or superiority to all glacier models.

Do not prescribe that levels 2–4 must be positive. If the existing maximum-spacing assumption is correct and the thermal extension changes nothing measurable, record that result and stop the proposed headline.

## Environmental conditions are a test, not a label

Three deliberately separate tiers are needed.

**Tier A: published laboratory benchmark.** Use the published inputs and measurements, with source/uncertainty provenance. A concrete anchor is Yokokawa's CSIM120910B: width 0.02 m, initial depth 0.00165 m, discharge 1.329e−5 m³/s, water temperature 0.7 °C, slope 5 degrees. These imply U≈0.403 m/s and Fr≈3.17. Its reported ice-probe temperature is not automatically a lower-boundary temperature. The experiment's heat and wall conditions must be reproduced, not inferred from one number. Most runs are laminar/transitional; keep them separate from turbulent applications.

**Tier B: physically constrained, synthetic open reaches.** Build cases from depth, discharge, width, slope, roughness, inlet heat and ice temperature that satisfy their governing relationships. Do not independently sample all of them as though they were unconstrained. Include temperate and cold ice; steady and time-dependent forcing; inlet temperature depletion; realistic heat budgets; and wet/dry limits where supported. “Environmentally plausible scenario” is the correct label until linked to actual measurements. Do not call a synthetic grid a Himalayan dataset.

**Tier C: observational transfer.** [Kamintzis et al.](https://doi.org/10.1080/04353676.2018.1549945) document discharge-dependent relocation of erosion in a real supraglacial knickpoint. This motivates the output, but an isolated knickpoint is not a periodic step train and its published description is not a complete CFD boundary dataset. Himalayan drainage observations at [Lirung Glacier](https://doi.org/10.3389/feart.2017.00069) establish drainage-network relevance, not validation of our mechanism. No claim of matched Himalayan validation is justified by these sources alone.

The research can be conducted without collecting new field data: mathematical results and verified computational experiments are legitimate. They establish model behavior under declared assumptions. They cannot, by computation alone, establish unknown environmental heat-transfer coefficients or observational forecast accuracy.

## The thermal-timescale calculation already done

The exact flat-slab response in MATHEMATICS.md shows when steady conduction is inadequate. As an **illustration**, using kappa_i=1.2e−6 m²/s, an 8 cm slab driven at a six-hour period has a response magnitude 1.169 times its static value and an admittance phase of 25.8 degrees. At 24 hours those become 1.012 and 7.34 degrees. A 20 cm slab at six hours has a magnitude ratio 3.091 and phase 46.3 degrees.

These are calculated thermal admittances, not melt enhancements or instability results. Slab depth is a specified thermal-domain parameter, not arbitrary glacier thickness. The actual basal boundary and initial temperature profile determine whether this idealization applies. The zero-cold-gradient limit must be included; the proposed cold-memory feedback vanishes there.

## What would make the eventual result strong

* A stable/unstable spacing boundary reproduced by two numerical formulations, with a clear mechanism and a quantified numerical uncertainty.
* A prediction made before each independent nonlinear test, including cases near the boundary where wrong predictions are informative.
* A meaningful difference from the closest thermal baseline, not only from an obviously deficient no-melt model.
* An environmental output: integrated melt, downstream heat export, migration speed, or fraction of melting concentrated in a fixed fraction of the reach. Preserve the input energy budget when comparing cases.
* A result that survives at least two justified heat-transfer treatments, or a clearly documented explanation of which closure uncertainty prevents prediction.
* Reusable code with physical validity checks, dimensional inputs/outputs, state/checkpoint provenance and failure records. Operational usefulness begins here, not with an agency logo.

## Specific failure modes to prevent

1. A periodic warm-water boundary secretly supplies unlimited heat.
2. A fixed-temperature ice bottom stays artificially close while the modeled bed recedes substantially.
3. Jump locations are frozen during linearization, changing the stability calculation.
4. Grid viscosity or a chosen momentum regularization sets the reported spacing.
5. A spatial Floquet/Bloch calculation is described as a new mathematical technique.
6. An under-resolved hydraulic jump is interpreted as a heat-transfer experiment.
7. More localized melting is called more total melting without measuring the energy partition.
8. Laminar flume agreement is called turbulent field validation.
9. A weak local eigenmode is advertised as a robust nonlinear environmental effect.
10. The same equations implemented twice are called independent evidence for the heat-transfer closure.

## Scope and project management

The existing network is useful infrastructure, but its quasi-steady discharge treatment cannot provide the missing jump dynamics. Reuse its units, provenance, resource scheduler and conservation tests. Audit whether its independent HLL reference actually evolves streamwise bed morphology, supports moving hydraulic jumps and resolves water/ice heat transport before reusing numerical kernels. The current remote repository was not fetched or modified in this task; reuse recommendations are based on the status supplied by the user.

The package does not authorize an immediate large campaign. First establish one reproducible profile family and one verified perturbation calculation. That is substantial scientific implementation work, not a broad parameter sweep. Only then measure runtime and freeze a production design. A large number of runs cannot rescue an inconsistent base state.

The best next outcome is a precise new answer about the dynamics of melting ice. No responsible assessment can promise Nature publication or STS top 40 before that answer exists.
