# Heterogeneous inverse-flux gate

Research audit, 3 November 2026. This changes no frozen model or field score. Tower 2 camera outcomes were not opened. The earlier Tower 1 comparisons remain failures of all-baseline superiority.

## Decision

Proceed with a bounded DEVELOPMENT implementation of the positive, time-integrated layered conduction inverse, because it addresses a concrete failure of the homogeneous temperature fit. Do not describe its linear algebra, positive least squares, layered conductivity, or inverse heat flux as a new mathematical discovery. A possible methodological contribution is a source-audited **basal-flux feasible-set and abstention procedure**, tested without melt calibration against a genuinely withheld station. Its originality and field benefit are not established by this audit. An implementation alone does not establish either.

The root's inferred low middle-layer conductance followed by higher lower-layer conductance is a useful diagnostic, not independent evidence: the illustrative steady conductivities used the already-known TRAIN melt rate. A new fit must use temperature/geometry/C inputs only, with all choices frozen before Tower 2 camera scoring. Previously inspected Tower 1 held-out periods are development evidence now, not new tests.

## Closest primary work

1. **Nakawo and Young (1981; 1982)**, DOI [10.3189/172756481794352432](https://doi.org/10.3189/172756481794352432), and [1982 article](https://www.cambridge.org/core/journals/journal-of-glaciology/article/estimate-of-glacier-ablation-under-a-debris-layer-from-surface-temperature-and-meteorological-variables/3F6103F2F2EE1DFE63BCF9503731A1E1). The 1981 work explicitly combines dry and wet layers as series thermal resistance. The 1982 approach infers effective resistance from surface meteorology/temperature and then predicts melt. Neither permits a claim that series resistance, layering, or temperature-derived melt is new.
2. **Brock et al. (2010)**, DOI [10.1029/2009JD013224](https://doi.org/10.1029/2009JD013224), estimates effective conductivity from coincident mean melt flux, debris thickness and temperature contrast. This uses melt in the conductivity determination. A temperature-only test is materially different but must not use a later fitted melt gain and call that absolute temperature-only calibration.
3. **Conway and Rasmussen (2000)** as implemented and analyzed in **Beck and Nicholson (2025)**, DOI [10.5194/tc-19-2715-2025](https://tc.copernicus.org/articles/19/2715/2025/), estimates diffusivity from vertical temperature arrays and obtains conductivity using an external volumetric heat capacity estimate. Beck–Nicholson explicitly examines interacting spatial/temporal truncation, sensor precision and position errors, and their impact on inferred melt. Therefore generic sensor-error correction or avoiding noisy finite differences is not a new research gap.
4. **Laha et al. (2022)**, DOI [10.1017/jog.2022.35](https://doi.org/10.1017/jog.2022.35), already includes homogeneous regression CRh, layered CRi, and forward Bayesian homogeneous/layered alternatives. CRi generalizes to a two-layer, three-parameter regression. The full primary article was read in the earlier audit, but its current Cambridge endpoint returns HTTP429; exact typography of equations should be checked from the retained paper before an equation-by-equation replication claim. Layered inference and Bayesian forward inversion are direct neighboring methods.
5. **Melo-Velasco et al. (2025)**, DOI [10.1029/2025JF008360](https://doi.org/10.1029/2025JF008360), sections 2.3.1.1–4 and Table 5, compares NYB/CRh/CRi/Opt at this very field site. It already reports disagreement between conductivity methods and measured melt. Opt uses melt observations and is not a temperature-only baseline. The source also establishes that the stake and profile were installed in the same pit; a 3 m stake is not a measured 3 m horizontal mismatch.
6. **Petersen et al. (2022)**, DOI [10.1029/2021JF006520](https://doi.org/10.1029/2021JF006520), shows depth-dependent apparent diffusivity and nonconductive heat transport in glacier debris. A fitted heterogeneous k can absorb convection or phase processes, and its positivity does not establish that it is a material conductivity profile.
7. **Gutman and Ha (2009)**, DOI [10.1016/j.matcom.2008.12.002](https://doi.org/10.1016/j.matcom.2008.12.002), explicitly studies identifiability of piecewise-constant 1-D conductivities from finitely many internal temperature measurements and boundary excitation. Conductivity identifiability with sparse sensors is established mathematical territory.
8. **Battaglia (2024)**, DOI [10.1016/j.ijthermalsci.2023.108840](https://doi.org/10.1016/j.ijthermalsci.2023.108840), develops thermal system identification using integration rather than differentiation to avoid a signal prefiltering step. Integration-based noise mitigation is not new in inverse heat conduction.
9. **Djatouti et al. (2017)**, [original conference paper](https://publications.ibpsa.org/proceedings/bs/2017/papers/BS2017_513.pdf), explicitly updates a thermal model for heat-flux quantities rather than global temperature error. **Spantini et al. (2017)**, DOI [10.1137/16M1082123](https://doi.org/10.1137/16M1082123), supplies goal-oriented Bayesian inverse theory. A flux-oriented objective is appropriate, not historically new by itself.

This is a targeted audit, not an exhaustive proof of absence of an identical glacier implementation. It supports a conditional implementation claim, not priority language.

## Exact integrated model

For nodes z_0<...<z_N, per-area P1 finite elements with known volumetric C give consistent element mass

    M_e = C_e h_e/6 [[2,1],[1,2]],   G_e=k_e/h_e >=0.

For an interior node i over [a,b], define S_j=integral_a^b T_j dt. Then

    (M [T(b)-T(a)])_i
      + G_(i-1)(S_i-S_(i-1)) + G_i(S_i-S_(i+1)) = d_i.

Here d_i is integrated omitted heating/model-discretization discrepancy, in J/m². Every term is energy per area. This is linear in G only when C, geometry, sensor-to-node mapping and observed interpolants are fixed. The consistent mass matrix includes neighbor endpoint changes; replacing it by a lumped diagonal is a different declared discretization, not algebraically the same method. A sum of equations can yield a total heat ledger; do not drop boundary fluxes while claiming conservation.

With fixed observed temperatures, nonnegative least squares estimates passive G. A small regression residual does not establish forward prediction, since both regressors and responses contain measurement error. Validate using temperature-only development splits and forward integration, rather than selecting the model by the inspected camera targets.

### Fundamental scale obstruction

If C is not independently fixed, simultaneous C→sC and G→sG leave all temperatures and integrated temperature equations unchanged but scale all heat/melt predictions by s. Regularization and positivity cannot resolve this. The independently sourced C proxy breaks scale only by assumption. A free output melt gain subsequently restores empirical melt calibration; it must not be described as inferred physical C or temperature-only absolute melt prediction.

### Safe conditional interval relaxation

Let every nodal temperature error satisfy |e_j(t)|<=delta_j throughout the interval. Let |d_i|<=eta_i represent an independently justified discretization/model discrepancy. For measured residual r_i(G), a necessary condition for an admissible true model is

    |r_i(G)| <= 2 sum_j |M_ij| delta_j
        + duration * [ G_(i-1)(delta_i+delta_(i-1))
                     + G_i(delta_i+delta_(i+1)) ] + eta_i.

The factor 2 in the endpoint term bounds two endpoint errors. This inequality is linear in nonnegative G after splitting the absolute value, so all windows produce a convex polyhedral outer set. It preserves any physically admissible model under the declared assumptions; it does not enforce one coherent error trajectory across overlapping windows, so it can be wider than the exact feasible set. It does not create independent information from overlapping windows. A time-constant sensor bias cancels from its own endpoint difference but persists in neighboring temperature gradients; explicit shared offsets may tighten the set.

For a known continuously melting basal node T_N=0, Q=G_(N-1) integral T_(N-1)dt. The observed integral has a bounded error duration*delta_(N-1). Since G>=0, linear objectives with lower/upper integral coefficients bound Q over the polyhedron. If the integral lower bound is negative, do not call its lower bound positive melt. If subzero basal ice, refreezing, or phase switches occur, this formula is not the total phase-change operator; the cold-ice complementarity model and its state/storage must be treated explicitly.

Unbounded LP means the observations do not constrain that flux under the declared model. Empty polyhedron means the model, geometry, error bounds or discrepancy bounds are incompatible. Neither outcome licenses retuning bounds until a desired field result appears. An NNLS point can exist even where a physically credible interval gate fails.

## What would be a useful contribution

A predeclared analysis can answer: **Can available thermistor profiles constrain the basal energy flux more tightly than they constrain internal conductivities, and can a physically passive heterogeneous model transfer its flux prediction to an uninspected station without melt calibration?** A non-identifiable k profile can still have a narrow Q interval; conversely a tiny temperature residual can leave Q unbounded. Demonstrating that distinction in these data with honest error bounds would be a precise result, though general inverse theory already recognizes it.

Necessary gates:

- Fix training periods, sensor geometry cases and independently sourced C/porosity/moisture cases before Tower 2 outcomes.
- Retain CRh, source-equivalent CRi where executable, homogeneous forward, constant rate and air-temperature baseline; report any baseline requiring melt calibration separately.
- Report both unscaled physical melt and any one-TRAIN-gain empirical prediction. Do not compare calibrated candidate with uncalibrated competitors without a separate fair-calibration analysis.
- Diagnose singular values/nullspace for the conductance matrix and the flux functional, not merely optimizer convergence.
- Bound or label unresolved spatial interpolation, sensor offsets, wet phase, depth uncertainty and moving geometry. A finite case envelope is not a continuous physical-coverage theorem.
- Use consistent mass/source budgets, independent synthetic layered solutions, and an out-of-fit forward-temperature check.
- Freeze all fitting and selection before any new camera reading; score all available eligible intervals with identical camera operators and source-missingness rules. One station supplies serial intervals, not independent glaciers.
- A physically based cold/nonconductive exclusion should be set from temperature/source metadata before camera access, not after seeing melt errors.

No result here establishes improved melt prediction, a new mathematical theorem, Nature suitability or award success. The implementable next test is materially better motivated than another unrestricted parameter search, and retains a genuine prospective validation opportunity if Tower 2 remains sealed.
