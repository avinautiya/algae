# Reserved temporal test of fixed-midpoint heterogeneous diffusion

This is a new hypothesis informed by the completed days14–28 study, not a claim of an independently conceived model. The initial benchmark cropped all arrays at28days and never scored later interior temperatures. Source availability/depth metadata were inspected, and the separate ablation study used later boundary forcing, but no days28–42 interior-temperature accuracy has been computed by this team before this freeze.

Include every previously eligible profile with complete native5minute observations through42elapsed days: KDT_2021a,b,d; KDT_2022a,b,c; KDT_2023a,b,c. Exclude2020a (41.517days),2021c (36days),2020b (constructed bottom). This leaves9profiles/23interior series. No temperature-value eligibility filter.

Train on hourly samples elapsed[1day,14days). Initialize at first observation once; no state updates using interior data thereafter. Evaluate hourly samples[28days,42days], inclusive, equal probe weights within each profile, equal profile weights in aggregate. Measured endpoints are contemporaneous forcing, not meteorological forecasts.

New model: common constant volumetric capacity, two conductivities with interface fixed at geometric midpoint between measured boundary probes. Fit independent top/bottom diffusivity from the Cartesian9x9logspace grid1e-8..1e-5m²/s. Finite-volume40cells, harmonic interface conductance, continuous heat flux, exact modal integration of linearly interpolated boundary input. Choose minimum training mean-squared interior error; deterministic first-grid tie-break. No later fitting, interface selection, priors or parameter adjustment. The coarse81point grid is deliberate and may limit accuracy.

Comparators: preserve original early-training homogeneous diffusion alpha, lag-filter taus, and ridge coefficients without refitting. Quasi-steady interpolation also retained. Include an additional algorithmic control: two-layer numerical solver with both alpha values equal to the ORIGINAL homogeneous fitted alpha. This control is not trained or selected and distinguishes discretization from heterogeneous material effects.

Primary result: every model's equal-profile MAE and RMSE and pairwise profile wins; no superiority claim unless lower aggregate MAE and wins at least8/9 versus every original comparator. This decision is a descriptive engineering gate, not a statistical significance test. Report all failure cases and all81training objectives. A fit on a grid endpoint must be flagged, not expanded using test results.

Numerical gate:40versus80cells for chosen two-layer pair and homogeneous FV control, maximum difference on hourly test predictions<0.02°C. If any fails, preserve all outputs, label numerical gate fail, and do not claim resolved superiority; a separate versioned refinement would be required. Equal-alpha analytic and independent steady series-resistance backend tests must pass before run.

No field melt, hydraulic-contact or cliff-lifetime validation follows. Temperatures identify conductivity/capacity ratios, not absolute heat flux. Fixed common capacity and midpoint may be inadequate. Profiles share one glacier/environment and are not independent climates. This temporal test is a reserved continuation of already inspected profiles, not cross-site generalization or a novelty-priority proof.
