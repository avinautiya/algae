# Next implementation milestone

Read RESEARCH_DECISION.md, MATHEMATICS.md and SOURCES.md first. This is a proposed research direction with verified subproblems, not a positive discovery that your simulations must reproduce. Preserve all earlier frozen results and negative findings.

## Objective

Establish whether a thermal ice-step profile can be reproduced and its perturbation stability calculated reliably. The longer-term hypothesis is that thermal coupling changes dynamically persistent step spacing and melt localization. Do not launch a broad production campaign in this milestone.

## 1. Protect and inventory the existing work

Read repo instructions; inventory running processes, memory/disk limits, branches and checkpoints. Use an isolated branch/worktree if needed to protect active jobs. Do not overwrite frozen studies. Locate the current HLL reference and verify its actual state variables, geometry update, energy convention, shock handling and thermal closures. Report exact file/function references. A solver that updates only a few incision parameters is not automatically a distributed morphodynamic reference.

Copy this packet into a clearly labeled research directory. Run verify_foundations.py unchanged and record its output. Keep the supplied results as baseline evidence and version any correction rather than silently replacing it.

## 2. Close the narrowly identified prior-art questions

Read Wu & Izumi 2023 and Wu's thesis formulation, solution and wavelength-choice sections. Read Balmforth & Vakil 2012 §§2 and 5 and Johnson et al. §§1.2 and 4. Inspect original ice-interface stability work and forward citations to the closest ice papers. Update an equation-level comparison: state variables, quasi-steady approximation, ice flux model, jump conditions, treatment of wavelength, actual stability calculation and validation.

If a matching thermal finite-amplitude stability analysis already exists, report the overlap and stop the novelty claim. Do not invent another layer merely to retain the word novel.

## 3. Build a reproducible baseline profile family

Implement the published Wu/Izumi model faithfully before extending it. Separate reproduction parameters from new-model parameters. Export each profile, full parameter set, units, residuals, critical-point treatment, jump states and the wavelength-family parameter.

Reproduce at least one published plotted profile and a small branch of admissible wavelengths. Digitized graphical comparisons must carry digitization uncertainty; do not claim roundoff-level reproduction from a printed figure. Continue the family using the previous converged profile; save each solution atomically. Do not choose the largest wavelength as the new model's answer.

For every profile, report the integrated water, heat and bed-recession residuals and explicitly account for any imposed atmospheric heat source. Do not silently reset water temperature at a periodic boundary. If a frozen-time cold boundary is used, compute the ratio of expected bed motion to thermal-domain depth during the intended analysis period.

Initial implementation targets, to freeze before new stability results: scaled smooth-equation residual <1e−6 and scaled hydraulic jump residual <1e−8, with convergence under doubling resolution. These are engineering verification thresholds, not physical error guarantees. If unattainable, diagnose/justify a versioned revision before running evaluations.

## 4. Verify the stability machinery in increasing difficulty

First reproduce the flat-slab admittance, its limits, and the conduction-only stable case from MATHEMATICS.md. Implement the flat-base coupled dispersion matrix and independently check its linearization. Do not claim these established mechanisms as discoveries.

Then linearize the reproducible traveling profile. Include jump displacement as an unknown. Check mass and momentum Rankine–Hugoniot conditions and their derivatives. Keep actual interface geometry continuous and include its perturbation matching. If thermal jump heating is neglected, quantify the jump Eckert ratio. Fix translation with a phase condition; verify the expected neutral mode where symmetry applies.

Use a direct shock-fitting generalized eigenproblem first if practical. An Evans determinant is optional; do not implement it solely to sound advanced. Include Bloch phases for perturbations spanning multiple steps. A reasonable development grid is beta in {0, pi/4, pi/2, 3pi/4, pi}, refined near any indicated transition. It is not sufficient for the final stability map.

Before any cold-memory extension, establish whether the underlying published thermal model's profile is stable or unstable. The full cold-ice operator on a curved boundary must be derived and discretized; the flat-slab Fourier multiplier supplied in this packet is only an exact benchmark. It cannot be substituted as the curved operator without an explicit approximation/error study.

If no stationary/relative-equilibrium thermal base is physically consistent, report that and formulate tangent evolution around a transient base. Do not force an eigenanalysis of a rapidly changing state.

## 5. Independently test a few perturbation predictions

Select at most three converged profiles spanning the admissible family. For each, pre-record the leading growth/decay prediction and expected jump displacement. Test with an independently discretized nonlinear solver and at least three perturbation amplitudes, checking that the linear regime is observed. Compare physical time, not iteration count.

Use successive spatial/time refinements and report numerical uncertainty of the measured rate. A mode counts as supported only if its sign is resolved beyond estimated numerical error and the predicted/measured rates are compatible with that uncertainty. Near-zero unresolved modes stay unresolved. No arbitrary positive real part becomes a discovery.

Measure runtime and memory on these cases. Set resource-aware per-case limits and durable checkpoints before planning the larger campaign. Do not run long silent jobs or commit files that are being rewritten mid-run; write to staging and atomically publish complete case records.

## 6. Deliver a decision, not a parameter sweep

Required files:

* equation/source comparison with exact code references;
* reproducible profile family and conservation/critical-point residuals;
* verified eigenproblem and independently measured perturbation rates, if reached;
* failure log and resource/runtime measurements;
* figure of profiles, admissible spacing family and resolved/unresolved stability results;
* verdict: reproduced / unresolved / contradicted, with exact remaining blockers.

Stop this milestone once those results exist or a concrete blocker is diagnosed. Do not start hundreds of environmental simulations, add a PINN, or promise a publication.

## Later campaign, only after this milestone works

Compare (A) maximum-spacing prescription, (B) stability in the existing thermal closure, (C) unsteady water heat transport and (D) cold-ice response. Disentangle changes to the base profile from changes to the perturbation operator. Use shared-base diagnostics where mathematically consistent, and separately report the fully re-equilibrated comparison.

For a stable-spacing candidate, test multiple-step evolution from several initial conditions; report a band if final spacing is not unique. Open finite reaches must test inlet heat depletion, realistic forcing and outlet conditions. Include temperate ice as a null case and distinguish laminar from turbulent closures. Compare against an independently resolved water/ice calculation on selected cases; the existing depth-averaged HLL model is only a model-class cross-check, not resolved turbulence or an independent closure validation.

Record integrated melt, melt-location concentration, migration, downstream heat export and hydraulic resistance. Compare with closure uncertainty and discretization uncertainty. Freeze the eventual test design only after development is complete. If differences are negligible or closure-dependent, state that instead of enlarging the model until something changes.

The eventual computational contribution may be a new thermal stability result and a usable reach-scale closure. It is not yet a validated glacier-wide forecast model. Preserve that distinction in every report and figure caption.
