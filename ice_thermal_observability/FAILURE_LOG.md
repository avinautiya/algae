# Numerical issue and versioned correction

The original frozen study `study_v1` completed all 120 smooth-profile comparisons, 80 independent vertex checks, refinements and scale cases. In one four-station thermal-age case (beta=0.11986456615013454), the separately solved lower/upper extrema were 0.4797997527779432 and 0.4797997527779305 W. Their 1.27e-14 W reversal produced a negative displayed radius. The result is numerically unresolved; it does not mean a negative uncertainty or exact continuum identifiability.

The correction preserves raw extremal values, orders near-equal values for presentation, and explicitly flags widths below the solver comparison scale. A larger ordering failure raises an error. A regression test covers this case. The primary forecast, dataset, seed, comparison methods, tolerances and protocol are unchanged.

`study_v1` is preserved as the original output. `study_v1_1` repeats the same protocol with the precision-reporting correction. Figures use the corrected version. This is not a new favorable test set or a change of the scientific decision rule.

The `study_v1_1` computation completed successfully, but the surrounding macOS `time -l` resource wrapper exited 1 because `sysctl kern.clockrate` is denied in the sandbox. Its partial output reports 1.90 s wall / 2.62 s user / 0.21 s system, with no peak-memory measurement. The shell therefore correctly stopped before the audit/plot commands. Those downstream commands were executed separately using the completed study; the computational campaign was not rerun for this wrapper failure.
