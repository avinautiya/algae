# Physics-constrained probabilistic retrieval and effect interpretation

This is a concrete proposed next model, **not a newly trained or validated model on this branch**. Implementing a larger neural network on 64 field rows would not establish operational capability. The executed ML baselines remain necessary falsification/comparison tools.

## Prior art and documented limitations

[SNICAR-ADv4 (2022), sections 3.3.2 and 4](https://tc.copernicus.org/articles/16/1197/2022/) reports visible-albedo disagreement for highly colonized ice, loosely constrained surface properties in several comparisons, and remaining work on wet/weathered ice and clustered impurities. These are documented limitations of that study/configuration, not proof that current successors retain every limitation.

[GA_BLOOM (2026)](https://www.nature.com/articles/s43247-026-03758-8) identifies poorly constrained loss/phenological parameters and challenges comparing point samples with its 10 km output. Improved optics alone will not repair those ecological uncertainties.

[NASA HLS](https://hls.gsfc.nasa.gov/algorithms/) illustrates the preprocessing required for intersensor reflectance consistency; [NASA Landsat validation](https://science.nasa.gov/mission/landsat/calibration-validation/) emphasizes independent higher-level-product evaluation. Neither endorses this project or requires a particular ML architecture. [JPL coupled atmosphere/snow retrieval work](https://www.jpl.nasa.gov/site/research/bohn/) is relevant prior art. Joint retrieval and uncertainty are not first-ever contributions here.

## Proposed generative model

For pixel p and date t infer a low-dimensional latent state:
`z = [algae abundance, pigment loading, dust loading, ice optical radius, liquid-water state]`.
Only include a component if the selected sensors/prior information can identify it. Do not fit all components freely to four bands. Freeze unsupported species fraction or return prior-dominated flags.

Shared parameters `theta` represent joint pigment/cell calibration, not independently resampled per-pixel values. Ice-density/crust and sensor/geometry nuisance parameters must have measured or independently justified priors.

The observation operator is:
`observed bands = SRF[geometry/atmosphere operator(RT(z, theta, illumination))] + discrepancy + measurement error`.

For BOA reflectance use the corresponding directional observation operator. For radiometer targets use hemispherical albedo/flux. Do not make a learned scalar transform silently equate the two. Explicitly represent sensor covariance and geometry uncertainty.

Begin with a small hierarchical Bayesian model plus a regularized low-rank spectral discrepancy or residual GP. Learn discrepancy on TRAINING observations only, with shrinkage and counterfactual tests. Compare against zero-discrepancy physics and ordinary ML. Check whether discrepancy absorbs algae signatures; if it does, attribution is not identifiable. A flexible discrepancy is not a license to label residual darkening as biological.

If acceleration requires a neural emulator, train it on verified radiative-transfer simulations over a physically supported state domain, not pseudo-labeled satellite pixels. Separate simulator-emulator error tests from empirical model validation. Enforce albedo bounds, retain signed scattering effects, and require the same zero-algae/dust counterfactual. Active sampling should target emulator uncertainty/error, not favorable empirical outcomes. No claim of solving a governing PDE via PINN is justified unless that actual formulation is implemented and tested.

## Training and independent evaluation sequence

1. Audit new PROMBIO records, match satellite/AWS support and freeze campaign/site-year roles before fitting. Keep new holdout targets out of every calibration, discrepancy and prior-selection step.
2. Quantify identifiability with simulated parameter-recovery tests, sensitivity/Jacobian conditioning and posterior-prior comparison. Simulation success does not establish real-data success.
3. Establish observation-operator correctness on ice/dust controls before enabling pigment interpretation.
4. Fit the small hierarchical model on the training campaigns. Inspect independent chains, effective sample sizes, predictive residuals, calibration and sharpness. Carry joint theta draws through retrieval and effects.
5. Evaluate untouched albedo/absorbed-SW observations against host-native/established algae optics and empirical satellite baselines on identical supports. New abundance labels alone do not establish radiation improvement.
6. Advance to a validated host surface-energy balance only after radiation validation. Independent observed ablation is a separate endpoint, not `RF / latent_heat` renamed as melt.

## Candidate novel product: qualified effect rankings

Combine biological evidence and paired radiative effects, rather than publishing a forced abundance leaderboard. Maintain two separate outputs:

- Biological interpretation: supported/ambiguous algae-compatible signal, with dust/water alternatives and prior dominance. These are evidence categories, not invented species labels.
- Physical effects: median/interval signed absorbed-SW or integrated-energy perturbations from matched algae-present/absent counterfactuals, and probabilities of practically meaningful pairwise differences.

`effects.py` implements the latter using aligned joint draws, explicit quantities/units, physical margins, and mandatory QC/domain/identifiability/forward-validation gates. It returns unresolved or unsupported comparisons instead of forcing a total order. It preserves common calibration uncertainty during regional area aggregation. Probabilities remain conditional on the model/ensemble, not externally calibrated causal probabilities. The code does not establish that these ideas are unique in the literature.

Do not rank species causally from an unidentified community fraction. Do not sum pigment-specific counterfactual effects as if nonlinear absorption had a unique additive attribution. Distinguish W m-2 intensity from W regional power and J accumulated energy. No melt ranking is allowed by this component.

## Acceptance and handoff

Run `python3 -m unittest discover -s agency_design/tests -v` and `python3 -m agency_design.prombio`.

This branch adds only research/design, original public data and tested interpretation/ingestion components. It does not replace the current retrieval or train a new operational model. The next implementer must resolve support matching and observation geometry, choose a qualified training/holdout inventory, and supply validated joint effect draws. Never enable `forward_validation_pass` merely because numerical tests pass.
