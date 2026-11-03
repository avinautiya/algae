# Independent pre-score review of Arolla interval models

No evaluation predictions or errors were computed by this reviewer. The source loader's5tests and model's9tests pass. The review checked the actual source workbook, metadata and code, not only a status summary.

* Source-provided MELT[mmw.e.] and density900 are correctly used. H1D is correctly primary; the other3stakes are spatially related secondary cases with unknown local debris thickness. No whole-glacier or cliff validation follows.
* Each fitted coefficient uses only the first interval. The held-out target enters scoring only. The target values were exposed in the data audit; the protocol correctly describes the design as retrospective rather than blind.
* The shared middle visit hour is identical in calibration and evaluation. Every declared scenario has finite hourly forcing across both intervals and all hour-label assumptions. END/CENTER/START conventions are hypotheses, not source-certified bin semantics.
* No native snow value is silently filled. These simple baselines omit snow dynamics and can fail because of them. An event mask based on evaluation errors would be inappropriate.
* For a dose ratio f and reading errors e0/e1/e2, residual uncertainty is `(f+1)e1 - f e0 - e2`. The conditional per-reading absolute bound2(f+1)delta and per-interval bound(f+1)delta are algebraically correct. The source0.1m is not a stated sigma; the code does not invent a confidence interval.
* The81 combinations are declared engineering scenarios, not81 independent tests. A robust-win flag over those scenarios is not a guarantee over every continuous unknown visit time or over all forcing/model errors. Dose-ratio extrema for fixed nonnegative drivers occur at visit-hour bounds, but pairwise model rankings need not be globally certified by a finite grid.
* Constant rate, PDD and net-shortwave dose are established empirical baselines, not novel physics. A favorable holdout would support a limited observational comparison, not Nature-level originality or universal superiority.

## Metadata correction before publication/scoring

Cell-level DOCX inspection showed the explicit Miage borrowing note belongs to roughness z0. The conductivity k=0.96 row has an empty provenance cell. The loader's previous association of that note with k was corrected. Original parsed data remain preserved; every forcing value and target is unchanged. See `PROVENANCE_AMENDMENT.json` and `AROLLA_PROPERTY_TABLE_ROWS.json`. The models do not use k, so this correction does not alter their equations or selected cases. The canonical input hash must nevertheless be updated in the model freeze before publication.

## Disposition

No remaining must-fix was identified for the bounded retrospective, weather-conditioned ice-ablation comparison, provided the published canonical freeze uses the corrected input hash. The exact bin timing, uncertainty semantics, snow omission, one shared evaluation weather event and spatial weather-transfer assumptions remain explicit limitations. This is not certification that all errors are absent or that predictive accuracy will be good.
