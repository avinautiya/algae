# Results of the proof and executable demonstration

31 October 2026. All incision numbers below are from the provided synthetic reduced-model scenario. They are not field measurements, newly resolved CFD results or an accuracy improvement against a glacier dataset.

## Mathematical result

The proof establishes an exact map with meltwater mass addition, constructive sharp regional bounds, forcing-uniform attaining profiles and an optimal scalar worst-case forecast. It also identifies what exact prefix thermal information removes regional ambiguity. The statements, assumptions and proofs are in THEOREMS.md; mathematical priority is assessed separately in PRIOR_ART.md.

## Executable forecast

For the synthetic one-hour case, admissible profiles sharing the same total conductance and outlet predictions yield declared-region mean lowering between **0.52994 mm and 4.88185 mm**, a factor of **9.21** between the extremes. The feasible nominal profile gives 2.68073 mm. The scalar minimax estimate is **2.70590 mm**, with conditional worst-case absolute error **2.17595 mm**.

The interval width is a deterministic information limit for the declared class. It is not a confidence interval, a measured glacier uncertainty or an improvement factor. The two explicit witness profiles attain its endpoints within this model. Their realizability by a particular full fluid/ice geometry remains a separate physical question.

A declared **2 mm** regional-incision threshold can be reached as early as **1518.70 s** (25.31 min) by one attaining profile. Another admissible profile does not reach it within the one-hour horizon. The nominal forecast crosses at 2375.61 s (39.59 min). These are model threshold statements, not flood or hydraulic-failure predictions.

The nominal run melts 0.05362 m³ of ice. Its maximum water residual is approximately 1.48e-14 m³ and thermal-energy residual 2.24e-8 J, compared with 360 m³ water supplied and 3.54228e8 J inlet sensible energy. These are conservative accounting checks, not environmental accuracy measures.

## Independent numerical validation

* **22 unit/integration tests passed**, including actual CLI execution and output preservation.
* Exhaustive small-profile checks match all contiguous-region extrema, including exact internal information and interval total bounds.
* **500 continuous cases** were checked by independently enumerating intersections of three polytope faces. The largest discrepancy, normalized by inlet power, was 1.33e-17.
* An independent RK4 integration of the coupled mass/enthalpy equations converges at approximately fourth order toward the map.
* Tests cover zero thermal input, the small-Stefan limit, calibration inversion, forcing-uniform witnesses, scalar minimax risk and threshold classification.

Testing corroborates the implementation; the proof supplies the general conditional result.

## Paper figures

Six reproducible method figures use this synthetic case and include underlying data, captions and source hashes. A seventh figure replots completed OpenFOAM verification metrics from repository commit 9a52524. The stored timestep sensitivity is tiny relative to spatial error but does not support the literal earlier claim that every error metric changed by less than 1e-13; the source-based correction is recorded beside that figure.

All seven PNGs were visually inspected. A crowded tick layout was corrected before delivery. Vector PDF versions are included.

## Remaining scientific gate

The physically relevant next use is to derive justified conductance uncertainty from verified CFD/calibration and test whether the assumptions and interval are useful for real incision forecasts. A proved information limit is a mathematical contribution; it does not itself demonstrate a new buoyancy/instability mechanism or a Nature-level physical discovery. The active CFD work continues independently.
