# Melt-goal transferability audit

## Decision

The next defensible contribution is a **calibration-transferability certificate for phase-aware debris-to-ice heat models**, rather than another temperature-MSE fit or an assertion that moisture necessarily increases melt. Implement a source-bounded family of conservative models, retain models on TRAIN evidence only, and quantify the remaining diameter of held-out melt predictions before opening outcomes. A small diameter licenses a conditional point estimate; a large diameter supplies explicit counterexample models and an abstention. This is an implementable research hypothesis, not established historical originality, field superiority, or a Nature-level discovery.

## Closest primary comparisons

- [Laha et al.2022](https://www.cambridge.org/core/journals/journal-of-glaciology/article/estimation-of-ice-ablation-on-a-debriscovered-glacier-from-vertical-debristemperature-profiles/3C9505EF0F1D8E4EAF164FB796F12BF2) already infer homogeneous/two-layer thermal properties from debris temperatures using derivative and Bayesian approaches. New layering or Bayesian profile fitting cannot be the contribution.
- [Melo-Velasco et al.2025](https://doi.org/10.1029/2025JF008360) directly compare field-derived and melt-optimized properties. Their discrepancy is prior evidence, not a discovery of this project. Their field arrangement separates buried instruments from weather towers; the Tower3 tower terminates7December while buried observations continue.
- [Spantini et al., Goal-Oriented Optimal Approximations of Bayesian Linear Inverse Problems](https://kiwi.oden.utexas.edu/papers/goal-oriented-bayesian-inverse-problem.pdf) explicitly give a transient heat-equation inverse example. General goal-oriented thermal inversion is established. Their formulation estimates a distribution of desired quantities in a linear-Gaussian setting; a proposed distinction here would be deterministic bounded source errors, nonlinear melt/no-melt phase constraints, and calibration-transfer ratios, not the basic goal-oriented principle.
- [Renard et al.2010](https://doi.org/10.1029/2009WR008328) already explain how input/structural errors and nonidentifiability damage hydrological prediction. Equifinality and calibration failure are not novel concepts.
- [Giese et al.2020](https://tc.copernicus.org/articles/14/1555/2020/) account for moisture effects and show why conductivity alone is insufficient to infer net melt response. A universal positive moisture coefficient is physically unjustified.

This targeted review excludes broad originality claims. It does not prove that the precise proposed certificate has never appeared in inverse heat, robust control, or set-membership literature. Before a priority claim, compare full-text robust inverse/goal-oriented calibration papers, including the2025 article [Goal-oriented calibration of a model and its modeling errors](https://www.sciencedirect.com/science/article/abs/pii/S0045782525007583).

## Exact model-class result

Let theta denote thermal response shape, geometry, initial state and fixed phase rules. Let D_T(theta)>0 and D_V(theta)>=0 be integrated modeled TRAIN and future ice melt doses. A single local gain calibrated to positive observed TRAIN melt Y gives

    prediction_V(theta) = Y D_V(theta)/D_T(theta).

A common multiplicative scaling of both doses cancels exactly. It does not remove uncertainty in response shape, storage memory, geometry, or phase activation. The relevant identifiable quantity is therefore r(theta)=D_V(theta)/D_T(theta), rather than parameter uniqueness or training temperature error alone.

For an explicitly finite retained family F, define L=Y min_F r and U=Y max_F r. For every point p, max(|p-L|,|p-U|)>=(U-L)/2, attained by p=(L+U)/2. Thus radius>(declared acceptable error) proves that no point predictor has that uniform error guarantee over F. This is elementary minimax interval geometry, not a new historical mathematical theorem. If the true process is not in F, no coverage statement follows. An empty F means the model/error assumptions are falsified, not zero uncertainty. A finite grid is not a certified enclosure of a continuous parameter set.

For a linear fixed-state branch with response vector b_theta to future forcing and TRAIN dose a_theta, calibrated future predictions agree for every admissible forcing vector u iff b_theta/a_theta is identical across retained theta (on the linear span of admissible forcing). Proof: equality of all linear functionals implies equality of their coefficient vectors on that span; reverse implication is immediate. Changing phase branches invalidates a global linear-kernel argument and requires the nonlinear dose calculation.

## Executable counterexample and certificate

`transfer_certificate.py` implements the finite-family interval, minimax point, abstention and extremal witness IDs, using no glacier target data. Four tests pass. A passive first-order thermal memory model y'=(u-y)/tau starts at0, receives unit TRAIN forcing for1time unit and zero future forcing for1unit. Each tau is calibrated to TRAIN melt100. Tau .1/1/10 yield future11.110/108.616/187.206. All outputs are nonnegative, and exact TRAIN calibration does not ensure transferable predictions. This counterexample is about memory and one scalar calibration; it is not a claim that all three models fit the actual thermistors.

## Concrete next algorithm

1. Freeze geometry, source uncertainty, admissible closures and model family before any new outcomes. Include constant/layered conductive models, source-consistent ice cold content, and explicitly bounded initial states. Do not claim missing wet-pore phase or advection is covered.
2. Fit/retain thermal response shapes only using TRAIN thermistors and declared depth/temperature errors. Keep both accepted and rejected families. Do not choose the member with best evaluation melt.
3. Use the same one positive TRAIN melt gain for each model and each strong empirical comparator; report gain interpretation as empirical when it breaks physical scaling.
4. Before target release, compute the prediction diameter, witness models, energy/phase diagnostics, numerical sensitivity, and source-coverage gates for each period. Freeze a tolerance with scientific units, not one selected for favorable coverage.
5. Evaluate accuracy AND selective coverage: all-period error, non-abstained-period error, fraction abstained, and whether retained intervals cover observations. Never report selective error without coverage, or call a finite-family envelope a statistical confidence interval.
6. The falsifiable scientific question is whether TRAIN temperature agreement plus one local calibration controls the future melt functional across thickness/moisture/forcing changes, and which independent measurements shrink its diameter. A negative answer can be rigorous, but does not establish a superior point predictor.

## Tower3 eligibility without outcome release

Feature-only audit confirms native Tower3 weather ends2022-12-06T23:50. Native debris temperatures continue2023-03-06T13:00; moisture continues2023-03-02T13:00. Therefore Dec2–Dec16 TRAIN plus later14-day evaluation cannot be a strict-local-weather replication. Initial all-depth warmup and the top probe's buried location remain setup/geometry gates. No Tower3 camera values were read.

A distinct reconstructed-forcing contract is feasible: select standing Tower1 as forcing source before outcomes; fit Tower3/Tower1 air relation on their early overlap; freeze all models to the same reconstructed series and transfer-error scenarios. The reproducible chronological half-split audit gives MAE1.040C, bias−.902C, RMSE1.330C on the later overlap. This error is large relative to source instrument accuracy and its persistence after tower fall is untested. It must be propagated, not treated as exact local measurements. Do not optimize spatial relation using melt. A Tower3 test remains same-glacier/same-campaign spatial transfer, not an independent climate or Himalayan validation.

## Alternative public targets if the new contract is rejected

[Lirung2013–2014](https://zenodo.org/records/4446387) provides a52.4kB paired debris-temperature/stake archive from three expeditions. Metadata says two temperature sites and three seasons. Co-location, stake intervals, source clock, thermal coverage and independent parameter leakage still need auditing before eligibility; no target values were opened here. [Kennicott2016–2023](https://doi.org/10.1002/gdj3.70032) provides richer paired profiles and ablation, but published common heat capacity incorporates melt-colocated2023 profiles; reusing it on a2023 holdout would leak target information. Derive training-only capacity or use independently justified physical assumptions. Neither dataset is automatically a clean test simply because a download exists.

## Rolling weak-form candidate audit

[Evatt et al.2015](https://nora.nerc.ac.uk/id/eprint/512051/1/s2.pdf), discussion, explicitly identify estimating time- and depth-varying debris conductivity as desirable future work. That is a rationale, not evidence their paper implements the proposed inverse method. Online thermal inverse estimation is already established outside glaciers: [2018 EKF–RLSE heat-flux reconstruction](https://doi.org/10.1016/j.infrared.2018.10.023) and [2026 thermal-network parameter updates](https://doi.org/10.3390/en19132999) defeat broad claims to first adaptive thermal estimation. These two comparisons were inspected at primary abstract/formulation level, not reproduced.

A narrow candidate claim could be: causal rolling integrated-energy inference of positive layer conductances, with one TRAIN-only local melt scale and a pre-release calibration-transfer certificate on source-bounded geometry. No reviewed source here establishes that exact combination on glacier melt, but this limited review cannot certify first publication. Its environmental value depends on untouched stake outcomes and on demonstrated selective reliability, not algebraic complexity.

Important implementation checks: (i) a past7-day window must end before the predicted interval; centered/future windows are smoothing, not forecasting; (ii) using temperatures during an evaluation interval to update conductances or infer its heat flux is an observation-conditioned retrieval/hindcast even when the inverse is rolling; (iii) fixed positive heat capacity and positive conductances preserve the thermal network's passive structure, but fitted conductance can absorb omitted advection, latent moisture transport or sensor errors and is not necessarily material conductivity; (iv) variable heat capacity requires differentiating C*T or a justified enthalpy, rather than quietly using C*Tdot; (v) no-contact/cold-interface states require their own phase/geometry contract; (vi) select regularization using TRAIN thermal forecast, never Tower3 melt or the model having the smallest eventual melt error.
