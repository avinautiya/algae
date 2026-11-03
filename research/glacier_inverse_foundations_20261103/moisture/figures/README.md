These four figures reproduce the released nine-method Tower2 comparison. Run
`python3 make_publication_figures.py` from `moisture/` to regenerate PNG and
editable SVG outputs. The script checks prediction/result hashes, recomputes all
nine period MAEs from the frozen period records, and checks the five existing
development MAEs. It performs no parameter search or new model selection.

1. `figure1_all_nine_mae`: all nine methods, sorted for display. The original
   primary remains red and is never replaced by the best-scoring baseline.
   Primary MAE is13.353mm water equivalent; air-only7.265 and calendar-adjusted
   air3.435. There are four serial periods from one stake and one glacier.
2. `figure2_period_predictions`: primary, air and calendar predictions against
   the same four observed ice-ablation increments, with signed-error panel.
   All intervals have source-local03:00 endpoints. Connecting lines aid reading;
   they do not introduce observations between periods.
3. `figure3_development_transfer`: all nine development and transfer MAEs. Each
   stake has its own original TRAIN amplitude calibration; global parameters
   come only from Tower1. Development errors are not fresh validation. The
   legend is outside the plotting area and obscures no bars.
4. `figure4_moisture_transfer`: descriptive daily SM3 anomalies and the frozen
   exponential-versus-air prediction multiplier. The same beta5.4 reduces
   Tower1 predictions but increases Tower2 predictions. This algebraic effect
   does not establish a causal moisture effect. SM3 probes sample different
   nominal depths (Tower1 approximately0.50m, Tower2 approximately0.39m), and
   operating/calibration flags remain. Boundary-day daily means use available
   samples. Tower1's30December boundary is00:00; Tower2's is03:00.

All ablation axes use mm water equivalent. These are empirical conditional
hindcasts using observed covariates, not future-weather forecasts or an
energy-conserving glacier simulation. Camera precision is unknown; no error bars
or confidence intervals are fabricated. The source campaign is described in
[Melo-Velasco et al.2025](https://doi.org/10.1029/2025JF008360), with
[data/code](https://doi.org/10.5281/zenodo.11581073).

`FIGURE_NUMBERS.json` contains plotted aggregate numbers and fixed-parameter
development predictions. `FIGURE_MANIFEST.json` binds source/script/output hashes
and records visual inspection of the final four PNGs. SVG text remains editable.
