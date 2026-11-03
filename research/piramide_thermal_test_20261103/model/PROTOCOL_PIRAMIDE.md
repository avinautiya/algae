# Preregistered conditional sub-debris melt prediction campaign

This is a **new held-out chronological field-ablation comparison**, not an
extension of Arolla's inspected target or proof of novelty. Process data are the
independently decoded 2022–23 Pirámide campaign, DOI10.5281/zenodo.11581073;
measurement geometry/error is from the final paper DOI10.1029/2025JF008360.
Native camera outcomes must remain sealed until this protocol, executable
sources and outcome operators are published. Inspection of process metadata and
min/max temperatures is acknowledged; this is not blind in that broad sense.

## Fixed calendar and information

Primary tower is Tower1. Initialize from its measured thermal profile at
**2022-12-02 03:00**, more than seven days after deployment. Thermal training
ends at **2022-12-16 02:30**; the final half-hour state is propagated forward.
Camera gain training is Dec2 03:00→Dec16 03:00. Evaluation is
**Dec16 03:00→Dec30 03:00**, with no reinitialization or calibration against its
outcomes. These are equal14-day blocks chosen from feature/clock coverage before
camera values, not selected by an outcome or fit. Missing exact camera endpoints
refuse scoring; do not shift dates toward favorable observations.

The native source clock has no timezone. Keep naive labels from the same source
and explicitly assume clock agreement. Never assign UTC without authority.
Thermistors are native30-minute samples; no hourly aggregation or gap filling.
Left-sample hold is primary; right hold is an explicit reconstruction sensitivity.
The two interpretations are numerical assumptions, not a claim the source
observations are hourly/half-hour means. Both are retained.

## Physical candidate and geometry ambiguity

Solve the constant-diffusivity, signed sensible heat equation between the buried
upper thermistor and an assumed temperate ice interface. Dirichlet upper
temperature is observed DT1; assumed bottom temperature is0C. DT5 is an interior
thermistor, **never identified with the ice interface**. The source nominal
upper0cm sensor was covered by about2cm gravel. Absolute sensor depths and debris
thickness carry source-reported ±2cm errors. The following explicit static
geometry scenarios all lie within those marginal bounds; they are assumptions,
not independent measured exact geometry, a joint confidence region or proof
that every intermediate geometry is covered:

| Scenario | Upper m | DT2/3/4 m | DT5 m | Ice m | Ice-to-DT5 gap |
|---|---:|---|---:|---:|---:|
| small_gap_5mm |.020|.125/.250/.375|.495|.500|.005m|
| middle_gap_20mm **primary** |.020|.125/.250/.375|.500|.520|.020m|
| large_gap_40mm |.020|.125/.250/.375|.480|.520|.040m|

Initial temperature is piecewise-linear interpolation of all five initial
measurements plus the assumed0C ice boundary. Do not clip negative measurements.
The new signed-temperature solver preserves them; it does not include pore-water
freezing latent heat. Its physical use is conditional on negligible phase-change
effects and constant effective geometry. Count all negative basal-flux substeps;
if any occur, that scenario is not eligible for a no-refreezing melt conversion.
No unfavorable interval or scenario is deleted. Report phase/geometry scope
beside field errors; source-bound geometry does not independently verify it.

The solver is cell-centered finite volume, backward Euler,32cells, maximum300s
steps split exactly at native30-minute boundaries, with half-cell Dirichlet
resistances and a normalized discrete storage/flux ledger. Record refinement
separately before any final numerical-accuracy claim; no target-guided change to
these settings is allowed. It computes heat normalized by bulk heat capacity:
units K·m. This is **not** joules, mm water equivalent or a measured basal flux.

## Training-only diffusivity and established comparators

For each geometry/reconstruction scenario, search the fixed25-point logarithmic
grid `1e-8 … 3e-6 m²/s`. Fit one diffusivity using only DT2–DT5 on the training
rows, minimizing squared distance outside ±.4C source sensor-error intervals.
This is an engineering error interval, not a Gaussian likelihood. Exact ties
choose the smallest diffusivity. Report boundary optima and flat compatible
sets; do not extend the grid because evaluation errors look poor.

The source-independent derivative comparator regresses `T_t=αT_zz` on training
interior time/depth nodes, using signed temperatures and actual geometry. It is
explicitly **zero-intercept ordinary least squares**, not the paper's normalized
total-least-squares CRh/CRi method. A nonpositive estimate is recorded as refusal,
never clipped into the grid. It cannot support superiority over the source's
stronger estimators or professional glacier models. The campaign must retain a
temperature-only no-memory/equilibrium and simple melt baselines in the separate
analysis, with their information/calibration budgets declared.

## Camera operator and absolute melt calibration

The camera loader independently gates exact source timestamps and raw exposed
stake `distance` in metres of ice. The approved operator is
`(distance_end-distance_start)*917` mm water equivalent. Provided `ablation`
column units are conflicting and are prohibited. Signed differences and resets
are flagged, not clipped or silently repaired. Camera precision is not supplied
as a statistical standard deviation; no confidence interval is invented.

The thermal model first freezes normalized heat traces and training/evaluation
totals **without reading any camera outcome**. A separately published sealed
routine may then fit **one positive effective melt gain** to the TRAIN camera
interval: `gain = measured_train_mmwe / normalized_train_heat_K_m`. Freeze
`predicted_eval_mmwe = gain * normalized_eval_heat_K_m` before inspecting later
camera truth. The gain is empirical, not an independently measured bulk heat
capacity; its implied `C = gain*Lf` is a physical diagnostic, not justification
to select a favorable geometry. All geometry/hold results must be reported.
Constant-rate/PDD comparators may fit one gain to exactly the same training
camera information. No evaluation melt rescales conductivity, depth or gain.

The primary model is fixed **middle_gap_20mm + left hold + training forward fit**.
All other geometries, holds and derivative estimates are sensitivities/comparators.
No winner may be chosen after held-out melt is read. Primary metric is absolute
error of the single14-day evaluation ice-ablation increment. This is one shared
weather event, not independent multi-glacier evidence. Subsequent daily interval
metrics are secondary and cannot become the primary claim by postselection.

## Execute without outcome access

```
python3 -m unittest discover -v
python3 piramide_experiment.py freeze PIRAMIDE_EXPERIMENT_SPEC.json --out NEW_FREEZE_DIR
# Publish frozen artifacts; preserve their SHA and commit before proceeding.
python3 piramide_experiment.py predict NEW_FREEZE_DIR --commit PUBLISHED_40_CHARACTER_SHA
```

The output contains no camera targets, no evaluation-temperature score and no
absolute melt prediction yet. Parent analysis owns the separate training gain,
prediction publication and later target unlock. Failed phase/clock/geometry or
positive-gain gates are results, not permission to alter the frozen experiment.

This campaign may establish conditional field skill or falsify the predictor.
It does not assume a novel heat equation, a new environmental mechanism,
professional-model superiority, coupled cliff validation or publication success.
