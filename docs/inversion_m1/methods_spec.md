# Methods spec: smallest chemistry-aware extension (specification only; not yet implemented or run)

Baselines (existing methods, fixed optics): **A** = our BioSNICAR bounded least-squares inversion (`inversion_m1/baseline_fixed_optics.py`, config frozen in `config_baseline_A.json`);
**E** = the published Dark Zone target density re-implemented in NumPy (`inversion_m1/reproduce_darkzone.py`, MAP not NUTS).
The extension **B** modifies A only. E is a reference for what an existing published method achieves; B is not built on E because E's emulator cannot accept altered optics.

## Parameter classes

| Class | Parameters | Rule |
|---|---|---|
| Fixed (never fitted) | RT solver, ice refractive index, column geometry/densities, illumination rule, end-member absorption spectra (shape), SSA-radius relation, observation operator | identical to A |
| Measured (enter as data, with their own uncertainty) | sample-level absorption spectra where they exist (e.g. dust MAC n = 10 in the Dark Zone paper; ice-algae in vivo spectra; two cryoconite reflectance-derived spectra) | bounds of B come only from such measurements, never from the albedo data being fitted |
| Inferred per surface (as in A) | log10 algal cells, log10 mineral, log10 dark-particle loading, ln bubble radius | bounded as in config A |
| **New, inferred per surface (B)** | **at most one scalar per absorber class that has >= 2 independent measured end-members**: mixing weight w in [0, 1] between the two measured end-member spectra (`k(lambda) = (1-w) k1 + w k2`); and for dust a scalar mass-absorption multiplier m in the measured range | chemistry enters only as a convex combination of measured spectra or a bounded scalar of a measured spectrum; no per-wavelength degrees of freedom, no spectral offsets, no smoothing splines |
| Transferable (carried to unseen glaciers unchanged) | the end-member spectra, the bounds of w and m, and the error-scale rule | learned or set from development glaciers/measurements only |
| Not transferable | per-surface loadings, w and m, radius | re-inferred for every spectrum |

## Ambiguity and unresolved states
* w and algae/dark-particle loading can trade off (the Dark Zone authors report cryoconite versus dense-crust similarity as an unresolved confusion). For every spectrum report the **profile likelihood over w** and the width of the 95 % set; if the set covers more than 50 % of the prior range, w is **reported unresolved** and excluded from any claim that uses it.
* Estimates at a bound are flagged, not clipped quietly.
* Saturated, shadowed or liquid-water-dominated spectra are labelled as an out-of-scope state (A's quality rules plus a residual-pattern flag); they are not rescued by extra freedom.
* Report separately: (i) concentration-like outputs, which have no independent truth here, and (ii) reconstruction and prediction of held-out spectral regions.

## Data requirement before B is frozen
B needs **two measured end-member spectra** for at least one absorber, with provenance. At the time of writing the archive contains one ice-algae spectrum, one dust mixture and one cryoconite spectrum (the cryoconite derived at S21 from two reflectance measurements). **Two independent measured end-members are not demonstrated to be in hand**; locating them (e.g. in the 2023 ice-algae paper's per-sample data) is the first action of the next milestone. Until then B cannot be specified numerically and is not run.

## No chemistry calculation
No TD-DFT or MD is justified. The missing observable is a measured end-member spectrum, which is an experimental, not computational, gap; a calculation would also have no independent validation target here.


## Update after the matched-data check
See `extension_decision.md`: the only supported chemistry variable is one bounded scalar on per-cell algal absorption (bounds from the Williamson per-sample phenol per cell); it is degenerate with cell count except where pigment per cell and count are measured on the same sample, which they are not. The two-end-member mixing weight w is not supported by any in vivo, sample-paired data and should be treated as unspecified.
