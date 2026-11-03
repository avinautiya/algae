# Fixed-model future batch after the first reported field test

The first version3 field event is known: the primary did not beat the
air-temperature baseline. That failed all-baseline comparison is retained.
This new batch tests **unchanged** models and their already fitted training gains
on later camera intervals whose values have not been inspected. Prior test
knowledge is disclosed; the new future target values remain sealed.

## Calendar and source-only selection

Predeclared nominal14-day anchors were Dec30 03:00, Jan13 03:00, Jan27 03:00,
Feb10 03:00 and Feb24 03:00. Use the separately audited rule: nearest finite
native camera endpoint within±12h, tie earlier, no interpolation and no outcome
performance selection. The known first start is the already published
version3 amendment, **Dec30 00:00**, because the nominal03:00 target was missing.
The four later anchors are exactly native and finite, from a missingness-only
audit that exposed no future target magnitudes.

Before opening future outcomes, require source forcing coverage for every
predeclared comparator. The air-temperature source ends Feb21 21:40. Consequently
the entire Feb10→Feb24 block is **excluded before scoring**: filling its final
53.33h would exceed the established18h gap limit. No February24 camera target
will be opened for this campaign. Thermal forcing alone covers it but cannot
justify dropping the strongest air baseline.

The three retained periods are:

1. 2022-12-30 00:00→2023-01-13 03:00;
2. 2023-01-13 03:00→2023-01-27 03:00;
3. 2023-01-27 03:00→2023-02-10 03:00.

Keep source-naive local labels; do not invent UTC. Endpoints lie exactly on the
native30-minute temperature grid. The first window is14days3h; models integrate
its actual duration rather than silently assuming all periods last14days.

## Immutable models, gains and state

All36 physical members retain their published version3 diffusivity, capacity,
geometry, upper-temperature reconstruction, ice properties and initial state.
The primary remains`forward__middle_gap_20mm__left__proxy`. No new training,
geometry selection, alpha inference, wetness choice or gain calibration is
allowed. The separate analysis reuses **all39** already frozen model gains:
36physical members plus constant-rate, buried-temperature and air-temperature
comparators. It must verify the exact version3 gain artifact hash.

The thermal runner starts each column at the original **Dec2 03:00** measured
debris profile and0C ice initial state and rolls it continuously to Feb10. No
state reset occurs at Dec30 or any later test boundary. It preserves the frozen
32debris/64ice cells,300s maximum step,2m ice domain, source-half-hour forcing,
fixed-geometry complementarity scheme and cold-content bookkeeping. The known
version3 TRAIN and first-evaluation raw-melt prefixes must reproduce before new
blocks are accepted. This checks lineage/state continuity without using a
future outcome.

For every retained period, output raw positive melt in mmwe, latent energy and
`normalized_dose = latent_energy/C_scene` in K·m. No new empirical gain is
estimated. Record cumulative sensible-plus-latent ledgers, complementarity,
cold-interface counts and raw ice-equivalent retreat. These physical ledgers
belong only to the unscaled raw component; multiplying by a retained empirical
gain does not create a new energy-conserving physical solve.

The future batch runner never imports camera values or camera/calibration
loaders. It hashes the fixed gain artifact as immutable lineage only. All future
point predictions, including the three statistical comparators, must be sealed
and published before their future camera endpoints are released.

## Comparison and interpretation

Primary metric is **mean absolute interval ice-ablation error across the three
new periods**, avoiding cancellation of positive and negative cumulative errors.
Compare the unchanged primary against each unchanged declared baseline on exactly
those intervals. All other family members remain secondary; no winner selection.
Also report a clearly separated pooled summary including the first failed
version3 event. Do not erase or relabel that previous result.

These are three adjacent periods at one tower, not three independent glacier
replicates. They share training calibration and an underlying continuous thermal
state. Observed future surface temperatures mean conditional hindcasts, not
weather forecasts. The same fixed-geometry, omitted ice-advection and wet-debris
pore-phase limitations remain, potentially more consequential over the extended
horizon. Neither successful numerical ledgers nor an improved score automatically
establish professional-model superiority, novel physics or broad field accuracy.

## Execute

```
# From model/: create and publish the NEW batch freeze before computing it.
python3 batch_phase_experiment.py freeze FUTURE_BATCH_SPEC.json --out frozen_piramide_future_batch
python3 batch_phase_experiment.py predict frozen_piramide_future_batch --commit PUBLISHED_SHA40
```

The parent analysis must separately apply the existing gains, seal predictions,
publish that seal and only then unlock the three future target periods. The
excluded fourth period remains unopened and unavailable in the final report.
